import numpy as np
import iminuit
import scipy.special as spec

# plot tau_s 
from scipy.constants import hbar, electron_volt, Boltzmann, e, pi
import math
from scipy.integrate import quad

print(f"hbar = {hbar:.5e} J·s")
print(f"e = {e:.5e} C")  # Coulombs
print(f"k_B = {Boltzmann:.5e} J/K")

#############################
# iMinuit fitting functions #
#############################

######################
# Mattis-Bardeen fit #
######################

def J_GR_eabs(f, n_qp_0, tau_r, V_ind, Delta):
    """
    Compute GR-absorbed energy spectral density J_Eabs_GR(f).

    Parameters:
    - f : float or ndarray
        Frequency [Hz]
    - n_qp_0 : float
        Quasiparticle density [1/m^3]
    - tau_r : float
        Recombination time [s]
    - V_ind : float
        Inductor volume [m^3]
    - Delta : float
        Superconducting gap energy [J]

    Returns:
    - JEabs : float or ndarray
        Energy spectral density [J^2 / Hz / m^3]
    """
    numerator = 2 * n_qp_0 * tau_r * Delta**2 * V_ind 
    denominator = (1 + (2 * np.pi * f * tau_r)**2)
    j_eabs_gr = numerator / denominator
    return j_eabs_gr

def convert_psd_eabs_to_dNqp(J_eabs, V_ind, Delta_0):
    """
    Convert energy spectral density J_Eabs(f) to quasiparticle number fluctuation PSD δn_qp(f).

    Parameters:
    - J_eabs : float or ndarray
        Energy spectral density [J² / Hz / m³]
    - V_ind : float
        Inductor volume [m³]
    - Delta_0 : float
        Superconducting gap energy [J]

    Returns:
    - dn_qp_psd : float or ndarray
        Quasiparticle number fluctuation PSD [1 / Hz / m³]
    """
    # scale_factor = 1 / (V_ind * Delta_0)
    scale_factor = 1 / (Delta_0)
    return J_eabs * scale_factor**2

def psd_dNqp_to_df_over_f(J_dNqp, alpha, gamma, kappa2, V_ind):
    """
    Convert PSD from delta n_qp to delta f_r / f_r,0.

    Parameters:
    - J_dnqp : array-like
        PSD of quasiparticle density fluctuations [1/Hz/μm³]
    - alpha : float
        Kinetic inductance fraction (unitless)
    - gamma : float
        Geometry factor (unitless)
    - kappa : float
        Material response parameter (unitless)
    - V_ind : float
        Inductor volume [m³]

    Returns:
    - J_df_over_f : array-like
        PSD of fractional frequency fluctuations [1/Hz]
    """
    factor = (alpha * gamma * kappa2 / (2 * V_ind))**2
    return J_dNqp * factor

def psd_dNqp_to_d1_over_Qi(J_dNqp, alpha, gamma, kappa1, V_ind):
    """
    Convert PSD from δn_qp to δ(1/Qi) using the specified parameters.

    Parameters:
    - J_dnqp : array-like
        PSD of quasiparticle density fluctuations [1/Hz/μm³]
    - alpha : float
        Kinetic inductance fraction (unitless)
    - gamma : float
        Geometry factor (unitless)
    - kappa1 : float
        Loss coupling coefficient (unitless)
    - V_ind : float
        Inductor volume [m³]

    Returns:
    - J_d1_over_Qi : array-like
        PSD of 1/Qi fluctuations [1/Hz]
    """
    factor = (alpha * gamma * kappa1 / V_ind)**2
    return J_dNqp * factor

def psd_df_over_f_to_ImS21(J_df_over_f, Q_r, Q_c):
    """
    Convert PSD from δf_r / f_r,0 to Im[δS21].

    Parameters:
    - J_df_over_f : array-like
        PSD of fractional frequency fluctuations [1/Hz]
    - Q_r : float
        Total quality factor (unitless)
    - Q_c : float
        Coupling quality factor (unitless)

    Returns:
    - J_ImS21 : array-like
        PSD of imaginary part of S21 fluctuations [1/Hz]
    """
    factor = (2 * Q_r**2 / Q_c)**2
    return J_df_over_f * factor

def psd_d1qi_to_ReS21(J_d1_over_qi, Q_r, Q_c):
    """
    Convert PSD from δf_r / f_r,0 to Im[δS21].

    Parameters:
    - J_df_over_f : array-like
        PSD of fractional frequency fluctuations [1/Hz]
    - Q_r : float
        Total quality factor (unitless)
    - Q_c : float
        Coupling quality factor (unitless)

    Returns:
    - J_ImS21 : array-like
        PSD of imaginary part of S21 fluctuations [1/Hz]
    """
    factor = (Q_r**2 / Q_c)**2
    return J_d1_over_qi * factor

def power_to_dbm(P_watts):
    """
    Convert power from Watts to dBm.

    Parameters:
    - P_watts: Power in Watts

    Returns:
    - Power in dBm
    """
    return 10 * np.log10(P_watts * 1e3)

def compute_p_feed(f_r, L, Q_i, N_0, Delta_0, rho_n, w_ind, t_ind, debug=False):
    if debug: 
        # Debug print
        print("Inputs:")
        print(f"f_r      = {f_r:.3e} Hz")
        print(f"L        = {L:.3e} H")
        print(f"Q_i      = {Q_i:.3e}")
        print(f"N_0      = {N_0:.3e} 1/(J·m³)")
        print(f"Delta_0  = {Delta_0:.3e} J")
        print(f"rho_n    = {rho_n:.3e} Ohm·m")
        print(f"w_ind    = {w_ind:.3e} m")
        print(f"t_ind    = {t_ind:.3e} m")
    Delta_0_J = Delta_0 * e  # Convert eV to J
    N_0_J = N_0 / e

    sqrt_term = np.sqrt(
        (pi * N_0_J * Delta_0_J**3) / (hbar * rho_n) ) * w_ind * t_ind

    P_feed = pi * f_r * L / Q_i * (0.42 * sqrt_term)**2
    P_dBm = power_to_dbm(P_feed)
    print(f"Readout power: {P_dBm:.2f} dBm")
    print(f"Readout power: {P_feed:.2g} J")

    return P_feed

def amp_psd(T_N, P_feed):
    """
    Calculate amplifier noise PSD contribution to Im[S21] or Re[S21].

    Parameters:
    - T_N      : Amplifier noise temperature [K]
    - P_feed   : Feedline power [W]

    Returns:
    - J_amp    : Amplifier noise spectral density [1/Hz]
    """
    return Boltzmann * T_N / (4 * P_feed) 

def calculate_inductance(C, f_r):
    """
    Calculate inductance L from capacitance C and resonant frequency f_r.

    Parameters:
    - C (float): Capacitance in Farads (F)
    - f_r (float): Resonant frequency in Hertz (Hz)

    Returns:
    - L (float): Inductance in Henrys (H)
    """
    L = 1 / ((2 * np.pi * f_r)**2 * C)

    return L 

def amp_psd_all(j_amp_ds21, Q_r, Q_c, V_ind, alpha, gamma, kappa_1, kappa_2, Delta_0):
    """
    Convert AMP PSD to dn_qp, df/f, d1/Qi, and J_eabs.

    Parameters:
    - T_N : Noise temperature [K]
    - P_feed : Feedline power [W]
    - V_ind : Inductor volume [m^3]
    - alpha, gamma, kappa_2 : resonator parameters
    - Delta_0 : Superconducting gap [J]

    Returns:
    - Dictionary of PSDs: {J_dn_qp, J_df/f, J_d1_Qi, J_eabs}
    """
    # Base AMP PSD

    # Directly passed through
    factor = (Q_r**2 / Q_c)**2
    J_d1_Qi = j_amp_ds21 / factor

    factor = (2*Q_r**2 / Q_c)**2
    J_df_f = j_amp_ds21 / factor

    # Convert to quasiparticle number fluctuation PSD
    J_dN_qp_diss = J_d1_Qi * ( V_ind / (alpha * gamma * kappa_1) )**2
    J_dN_qp_freq = J_df_f * ( 2 * V_ind / (alpha * gamma * kappa_2) )**2

    # Convert to energy PSD
    J_eabs_diss = J_dN_qp_diss * (Delta_0)**2
    J_eabs_freq = J_dN_qp_freq * (Delta_0)**2

    return {
        "J_d1/Qi": J_d1_Qi,
        "J_df/f": J_df_f,
        "J_dN_qp_diss": J_dN_qp_diss,
        "J_dN_qp_freq": J_dN_qp_freq,
        "J_eabs_diss": J_eabs_diss,
        "J_eabs_freq": J_eabs_freq,
    }

def update_tls_psd(T, T0, V_C, V_C0, E_C, E_C0, j_tls_0, beta, debug=False):
    """
    Update the TLS PSD J^{δfr/fr0}_TLS(1kHz) using scaling laws.

    Parameters:
    - T : float
        Target temperature [K]
    - T0 : float
        Reference temperature [K]
    - V_C : float
        Capacitor volume [m³]
    - V_C0 : float
        Reference capacitor volume [m³]
    - E_C : float
        Electric field [V/m]
    - E_C0 : float
        Reference electric field [V/m]
    - j_tls_0 : float
        Reference TLS PSD [e.g., 1/Hz]
    - beta : float
        Temperature scaling exponent

    Returns:
    - j_tls : float
        Scaled TLS PSD at 1 kHz
    """
    temp_ratio = (T / T0)
    volume_ratio = (V_C / V_C0)
    field_ratio = (E_C / E_C0)

    j_tls = j_tls_0 * temp_ratio**(-beta) * volume_ratio**(-1) * field_ratio**(-1)
    if debug: 
        print(f"T / T0**-beta = {temp_ratio**(-beta):.3e}")
        print(f"V_C / V_C0**-1 = {volume_ratio**(-1):.3e}")
        print(f"E_C / E_C0**-1 = {field_ratio**(-1):.3e}")
        print(f"Updated j_tls/j_tls_0 = {j_tls/j_tls_0:.3e}")
    return j_tls

def compute_W_res(L, N_0, Delta_0, rho_n, w_ind, t_ind):
    """
    Compute stored energy W_res in the resonator.

    Parameters:
    - L : float
        Inductance [H]
    - N_0 : float
        Density of states [1/J/m^3]
    - Delta_0 : float
        Superconducting gap [J]
    - rho_n : float
        Normal-state resistivity [Ohm·m]
    - w_ind : float
        Inductor width [m]
    - t_ind : float
        Inductor thickness [m]

    Returns:
    - W_res : float
        Stored energy [J]
    """
    pre_factor = 0.42
    Delta_0_J = Delta_0 * e  # Convert eV to J
    N_0_J = N_0 / e

    sqrt_term = np.sqrt((pi * N_0_J * Delta_0_J**3) / (hbar * rho_n))
    W_res = 0.25 * L * (pre_factor * sqrt_term * w_ind * t_ind)**2
    return W_res

def compute_E_field(W_res, C, t_aS):
    """
    Compute electric field amplitude E_C from stored energy.

    Parameters:
    - W_res : float
        Stored energy in the resonator [J]
    - C : float
        Capacitance [F]
    - t_aS : float
        Dielectric thickness (e.g., a-Si) [m]

    Returns:
    - E_C : float
        Electric field amplitude [V/m]
    """
    E_C = np.sqrt(4 * W_res / C) / t_aS
    return E_C

def full_psd_tls(f, J_1kHz, f_rolloff, n):
    """
    Compute TLS frequency noise PSD as a function of frequency.

    Parameters:
    - f : float or ndarray
        Frequency array [Hz]
    - J_1kHz : float
        TLS noise PSD at 1 kHz [unit consistent with output]
    - f_rolloff : float
        Roll-off frequency [Hz]
    - n : float
        TLS noise exponent

    Returns:
    - J_tls : ndarray
        TLS frequency noise PSD at frequency f
    """
    f = np.asarray(f)
    J_tls = J_1kHz * (1 / (1 + (f / f_rolloff)**2)) * (np.abs(f / 1e3))**(-n)
    return J_tls

def compute_f_rolloff(f_r, Q_r):
    """
    Compute TLS roll-off frequency from resonant frequency and loaded quality factor.

    Parameters:
    - f_r : float
        Resonant frequency [Hz]
    - Q_r : float
        Loaded quality factor (Q_r = (1/Q_i + 1/Q_c)^-1)

    Returns:
    - f_rolloff : float
        Roll-off frequency [Hz]
    """
    f_rolloff = f_r / (2 * Q_r)
    print(f"Roll-off frequency: {f_rolloff:.2e} Hz")
    return f_rolloff

def convert_dff_tls_psd_to_all(j_dff_tls, V_ind, alpha, gamma, kappa_2, Delta_0, Q_r, Q_c):
    """
    Convert TLS frequency noise PSD to various other observable PSDs.

    Parameters:
    - f_range : array
        Frequency array [Hz]
    - j_dff_tls : array
        TLS frequency noise PSD J_df/f [unitless/Hz]
    - V_ind : float
        Inductor volume [m^3]
    - alpha : float
        Fractional kinetic inductance contribution
    - kappa_2 : float
        Responsivity coefficient (unitless)
    - Delta_0 : float
        Superconducting gap [J]
    - Q_r : float
        Loaded quality factor
    - Q_c : float
        Coupling quality factor

    Returns:
    - dict of np.ndarray:
        Keys: 'J_dn_qp', 'J_eabs', 'J_Re(S21)', 'J_Im(S21)'
    """
    # dn_qp
    J_dN_qp = (2 * V_ind / (alpha * gamma * kappa_2))**2 * j_dff_tls
    
    # e_abs
    J_eabs = J_dN_qp * Delta_0**2
    
    # s21 terms
    s21_factor = (2 * Q_r**2 / Q_c)**2
    J_s21_imag = j_dff_tls * s21_factor

    return {
        "J_dN_qp": J_dN_qp,
        "J_eabs": J_eabs,
        "J_Im(S21)": J_s21_imag
    }

def s_exponential(t, tau):
    """
    Generate a causal exponential decay signal:
        s(t) = 0 for t < 0
        s(t) = exp(-t / tau) for t >= 0

    Parameters:
    - t : float or ndarray
        Time (can be scalar or array)
    - tau : float
        Time constant of decay [s]

    Returns:
    - s : float or ndarray
        Signal evaluated at time t
    """
    t = np.asarray(t)  # ensures compatibility with both scalars and arrays
    s = np.zeros_like(t)
    s[t >= 0] = np.exp(-t[t >= 0] / tau)
    return s

def compute_gr_resolution(n_qp_0, V_ind, Delta):
    """
    Compute the GR energy resolution sigma_Eabs.

    Parameters:
    - n_qp_0 : float
        Equilibrium quasiparticle density [1/m^3]
    - V_ind : float
        Inductor volume [m^3]
    - Delta : float
        Superconducting gap energy [J]

    Returns:
    - sigma_Eabs : float
        GR energy resolution [J]
    """
    return np.sqrt(2 * np.pi * n_qp_0 * V_ind) * Delta

def compute_amp_resolution(tau_qp, T_N, P_feed, debug=False):
    """
    Compute amplifier-limited resolution (std dev) for Re[δS21] and Im[δS21].

    Parameters:
    - tau_qp : float
        Quasiparticle lifetime [s]
    - T_N : float
        Noise temperature [K]
    - P_feed : float
        Feedline power [W]
    - debug : bool, optional
        If True, prints intermediate debug information.

    Returns:
    - sigma_amp : float
        Amplifier noise resolution for both Re and Im parts of δS21
    """
    numerator = Boltzmann * T_N
    denominator = 2 * tau_qp * P_feed
    sigma_amp = np.sqrt(numerator / denominator)

    if debug:
        print("=== Debug: Amplifier Resolution Calculation ===")
        print(f"tau_qp   = {tau_qp:.3e} s")
        print(f"T_N      = {T_N:.3e} K")
        print(f"P_feed   = {P_feed:.3e} W")
        print(f"k_B      = {Boltzmann:.3e} J/K")
        print(f"Numerator   (k_B * T_N)     = {numerator:.3e} J")
        print(f"Denominator (2 * tau_qp * P_feed) = {denominator:.3e}")
        print(f"Result σ_amp = {sigma_amp:.3e}")

    return sigma_amp

def convert_amp_res_to_eabs_res(sigma_ds21,
                                V_ind, Delta_0, alpha, gamma,
                                kappa_1, kappa_2, Q_c, Q_r, debug=False):
    """
    Convert amplifier resolution in δS21 to energy resolution σ_Eabs via frequency and dissipation channels.

    Parameters:
    - sigma_ds21_real : float
        Amplifier-limited resolution in Re[δS21]
    - sigma_ds21_imag : float
        Amplifier-limited resolution in Im[δS21]
    - V_ind : float
        Inductor volume [m³]
    - Delta_0 : float
        Superconducting gap [J]
    - alpha : float
        Kinetic inductance fraction
    - gamma : float
        Geometry factor
    - kappa_1 : float
        Responsivity for dissipation [m³]
    - kappa_2 : float
        Responsivity for frequency [m³]
    - Q_c : float
        Coupling quality factor
    - Q_r : float
        Loaded quality factor

    Returns:
    - sigma_eabs_diss : float
        Energy resolution via dissipation readout [J]
    - sigma_eabs_freq : float
        Energy resolution via frequency readout [J]
    """
    prefactor_diss = V_ind * Delta_0 / (alpha * abs(gamma) * kappa_1) * (Q_c / Q_r**2)
    prefactor_freq = V_ind * Delta_0 / (alpha * abs(gamma) * kappa_2) * (Q_c / Q_r**2)

    sigma_eabs_diss = prefactor_diss * sigma_ds21
    sigma_eabs_freq = prefactor_freq * sigma_ds21

    if debug:
        print(f"Prefactor (diss): {prefactor_diss:.3e}")
        print(f"Prefactor (freq): {prefactor_freq:.3e}")
        print(f"σ[Im(δS21)] or σ[Re(δS21)] = {sigma_ds21:.3e}")
        print(f"Q_c            : {Q_c:.3e}")
        print(f"Q_r            : {Q_r:.3e}")
        print(f"Q_c / Q_r**2      : {Q_c/Q_r**2:.3e}")

    return sigma_eabs_diss, sigma_eabs_freq

def tls_variance(tau_r, J_tls_1khz, f_roll, deltaf):
    """
    Compute (σ^{δf_r/f_r0}_TLS)^2 using numerical integration.
    
    Parameters:
    - tau_r : float
        Recombination time [s]
    - J_tls_1khz : float
        J_TLS at 1 kHz [1/Hz]
    - f_roll : float
        TLS roll-off frequency [Hz]

    Returns:
    - sigma_tls_sq : float
        Variance of fractional frequency TLS noise
    """
    # Define integrand
    def integrand(f):
        return (abs(f)**0.5 * (1 + f/f_roll)**2) / (1 + (2*np.pi*f*tau_r)**2)

    # Integrate from 0 to ∞, double it for symmetry (|f|)
    integral_val, _ = quad(integrand, 0, deltaf, limit=500)
    print("TLS Integral =", integral_val)

    prefactor = tau_r**2 / (J_tls_1khz * np.sqrt(1e3))
    sigma_tls_sq = (prefactor * 2 * integral_val)**-1
    sigma_tls_dff = np.sqrt(sigma_tls_sq)

    return sigma_tls_dff

def convert_tls_res_to_eabs_res(sigma_dff,
                                V_ind, Delta_0, alpha, gamma, kappa_2, debug=False):
    """
    Convert amplifier resolution in δS21 to energy resolution σ_Eabs via frequency and dissipation channels.

    Parameters:
    - sigma_ds21_real : float
        Amplifier-limited resolution in Re[δS21]
    - sigma_ds21_imag : float
        Amplifier-limited resolution in Im[δS21]
    - V_ind : float
        Inductor volume [m³]
    - Delta_0 : float
        Superconducting gap [J]
    - alpha : float
        Kinetic inductance fraction
    - gamma : float
        Geometry factor
    - kappa_1 : float
        Responsivity for dissipation [m³]
    - kappa_2 : float
        Responsivity for frequency [m³]
    - Q_c : float
        Coupling quality factor
    - Q_r : float
        Loaded quality factor

    Returns:
    - sigma_eabs_diss : float
        Energy resolution via dissipation readout [J]
    - sigma_eabs_freq : float
        Energy resolution via frequency readout [J]
    """
    prefactor_freq = 2 * V_ind * Delta_0 / (alpha * abs(gamma) * kappa_2) 

    sigma_eabs_freq = prefactor_freq * sigma_dff

    if debug:
        print(f"Prefactor (freq): {prefactor_freq:.3e}")
        print(f"σ_input = {sigma_dff:.3e}")

    return sigma_eabs_freq

def convert_tls_res_to_dnqp_res(sigma_dff, alpha, gamma, kappa_2, debug=False):
    """
    Convert amplifier resolution in δS21 to energy resolution σ_Eabs via frequency and dissipation channels.

    Parameters:
    - sigma_ds21_real : float
        Amplifier-limited resolution in Re[δS21]
    - sigma_ds21_imag : float
        Amplifier-limited resolution in Im[δS21]
    - V_ind : float
        Inductor volume [m³]
    - Delta_0 : float
        Superconducting gap [J]
    - alpha : float
        Kinetic inductance fraction
    - gamma : float
        Geometry factor
    - kappa_1 : float
        Responsivity for dissipation [m³]
    - kappa_2 : float
        Responsivity for frequency [m³]
    - Q_c : float
        Coupling quality factor
    - Q_r : float
        Loaded quality factor

    Returns:
    - sigma_eabs_diss : float
        Energy resolution via dissipation readout [J]
    - sigma_eabs_freq : float
        Energy resolution via frequency readout [J]
    """
    prefactor_freq = 2 / (alpha * abs(gamma) * kappa_2) 

    sigma_eabs_freq = prefactor_freq * sigma_dff

    if debug:
        print(f"Prefactor (freq): {prefactor_freq:.3e}")
        print(f"σ_input = {sigma_dff:.3e}")

    return sigma_eabs_freq

def compute_total_resolution(gr, amp_freq, amp_diss, tls_freq):
    """
    Compute total frequency and dissipation energy resolutions by quadrature sum.
    
    Parameters:
    - gr : float
        Generation-recombination resolution [J]
    - amp_freq : float
        Amplifier-limited resolution (frequency channel) [J]
    - amp_diss : float
        Amplifier-limited resolution (dissipation channel) [J]
    - tls_freq : float
        TLS resolution (frequency channel) [J]
    - tls_diss : float, optional
        TLS resolution (dissipation channel) [J] (default=0)
    
    Returns:
    - total_freq : float
        Total resolution for frequency readout [J]
    - total_diss : float
        Total resolution for dissipation readout [J]
    """
    total_freq = np.sqrt(gr**2 + amp_freq**2 + tls_freq**2)
    total_diss = np.sqrt(gr**2 + amp_diss**2)

    return total_freq, total_diss

def compute_total_resolution_list(res_list):
    """
    Compute total resolution from a list of contributions.
    Handles both scalars and array-like terms.
    
    Parameters
    ----------
    res_list : list
        List of resolution contributions.
        Example: [array_like, scalar, scalar]
    
    Returns
    -------
    total_res : ndarray or float
        Quadrature sum of contributions. Shape matches the largest array input.
    """
    # broadcast everything to arrays
    arrays = [np.atleast_1d(item) for item in res_list]
    
    # find target shape = shape of the largest array
    target_shape = np.broadcast_shapes(*[arr.shape for arr in arrays])
    
    # broadcast each to target_shape
    arrays_broadcasted = [np.broadcast_to(arr, target_shape) for arr in arrays]
    
    # quadrature sum element-wise
    squared_sum = sum(arr**2 for arr in arrays_broadcasted)
    total_res = np.sqrt(squared_sum)
    
    # return scalar if all inputs were scalar
    if total_res.size == 1:
        return total_res.item()
    return total_res

def compute_delta_f(tau_qp):
    """
    Compute frequency bandwidth (Δf) from quasiparticle lifetime τ_qp.

    Parameters:
    - tau_qp : float
        Quasiparticle lifetime [s]

    Returns:
    - delta_f : float
        Bandwidth Δf [Hz]
    """
    delta_f = 1 / (2 * np.pi * tau_qp)
    return delta_f

def get_phase_at_freq(f_query, f_paa, theta_deg, method="interp"):
    """
    Get phase theta (in degrees) at given query frequency/frequencies.

    Parameters
    ----------
    f_query : float or array-like
        Frequency or list of frequencies to query.
    f_paa : array-like
        Frequency sweep array (must be sorted).
    theta_deg : array-like
        Phase values corresponding to f_paa (in degrees).
    method : str, optional
        "nearest" -> pick the closest frequency in f_paa
        "interp"  -> linear interpolation between frequencies

    Returns
    -------
    phase_out : float or np.ndarray
        Phase at the requested frequency/frequencies (in degrees).
    """

    f_query = np.atleast_1d(f_query)  # ensure array for processing

    if method == "nearest":
        idx = np.abs(f_paa[:, None] - f_query).argmin(axis=0)
        phase_out = theta_deg[idx]
    elif method == "interp":
        phase_out = np.interp(f_query, f_paa, theta_deg)
    else:
        raise ValueError("method must be 'nearest' or 'interp'")

    return phase_out if len(phase_out) > 1 else phase_out.item()

def get_unwrapped_phase_deg(s21_shifted, direction="positive"):
    """
    Compute the unwrapped phase of shifted S21 in degrees, 
    ensuring it is monotonic either clockwise (negative) 
    or counter-clockwise (positive).

    Parameters
    ----------
    s21_shifted : array-like (complex)
        Complex S21 values shifted so that the resonance circle center is at (0,0).
    direction : str, optional
        "positive" -> output phases increase counter-clockwise (0 → 360 → 720 ...)
        "negative" -> output phases decrease clockwise (0 → -360 → -720 ...)

    Returns
    -------
    theta_deg : ndarray
        Unwrapped phase in degrees, all positive or all negative depending on `direction`.
    """
    theta = np.unwrap(np.angle(s21_shifted)) * 180 / np.pi

    if direction == "positive":
        # Force all values to be >= 0
        if np.mean(np.gradient(theta)) < 0:  # flipped direction
            theta = -theta
        theta = theta - np.floor(theta.min() / 360) * 360

    elif direction == "negative":
        # Force all values to be <= 0
        if np.mean(np.gradient(theta)) > 0:  # flipped direction
            theta = -theta
        theta = theta - np.ceil(theta.max() / 360) * 360

    else:
        raise ValueError("direction must be 'positive' or 'negative'")

    return theta

def angle_diff_from_im_shift(x, dy):
    """
    Compute the angle difference in degrees between (x, 0) 
    and (x, dy) on a circle centered at (0,0).

    Parameters
    ----------
    x : float
        Real part of the baseline point.
    dy : float
        Imaginary shift from the baseline.

    Returns
    -------
    delta_deg : float
        Angle difference in degrees.
    """
    theta1 = np.angle(x + 0j)
    theta2 = np.angle(x + 1j*dy)
    delta_deg = (theta2 - theta1) * 180 / np.pi
    return delta_deg

def intercept_energy_offset(e_abs_list, r_eabs, amp_res):
    """
    For each absorbed energy, find ΔE such that r(E±amp_res) intersects
    the nominal r(E) curve, handling monotonic increasing or decreasing cases.

    Parameters
    ----------
    e_abs_list : ndarray
        Absorbed energy values [meV].
    r_eabs : ndarray
        Corresponding radius values.
    amp_res : float
        Amplifier resolution [meV].

    Returns
    -------
    delta_Es : ndarray
        Energy offsets ΔE for each input E.
    """
    delta_E_ps = []
    delta_E_ms = []
    for i, E in enumerate(e_abs_list):
        r_cur = r_eabs[i]
        r_target_p = r_cur+amp_res
        r_target_m = r_cur-amp_res
        # ensure r_eabs is increasing for np.interp
        if r_eabs[0] < r_eabs[-1]:  # increasing
            E_r_p = np.interp(r_target_p, r_eabs, e_abs_list)
            E_r_m = np.interp(r_target_m, r_eabs, e_abs_list)
        else:  # decreasing
            E_r_p = np.interp(r_target_p, r_eabs[::-1], e_abs_list[::-1])
            E_r_m = np.interp(r_target_m, r_eabs[::-1], e_abs_list[::-1])

        delta_E_ps.append(E_r_p - E)
        delta_E_ms.append(E_r_m - E)

    return np.array(delta_E_ps), np.array(delta_E_ms)

