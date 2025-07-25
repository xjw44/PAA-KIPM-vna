import numpy as np
import iminuit
import scipy.special as spec

# plot tau_s 
from scipy.constants import hbar, electron_volt, Boltzmann, e, pi
import math

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
    numerator = 2 * n_qp_0 * tau_r * Delta**2
    denominator = V_ind * (1 + (2 * np.pi * f * tau_r)**2)
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

def update_tls_psd(T, T0, V_C, V_C0, E_C, E_C0, j_tls_0, beta):
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
    j_tls = j_tls_0 * (T / T0)**(-beta) * (V_C / V_C0)**(-1) * (E_C / E_C0)**(-1)
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