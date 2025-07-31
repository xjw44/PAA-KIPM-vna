import numpy as np
import iminuit
import scipy.special as spec

# plot tau_s 
from scipy.constants import hbar, electron_volt
import math
from scipy.constants import Boltzmann, e

#############################
# iMinuit fitting functions #
#############################

######################
# Mattis-Bardeen fit #
######################

Boltz_k = 8.6173303E-5 # eV/K
Planck_h = 4.135667662E-15 # eV*s
kB = 8.617333262145e-5  # Boltzmann constant in eV/K
hbar_ev = hbar / electron_volt #ev*s

def s21_z_to_mag(s21_z):
    """
    Convert complex S21 values to magnitude in decibels (dB).

    Parameters:
    - s21_z : complex or array-like
        Complex S21 value(s), e.g., from simulation or measurement.

    Returns:
    - s21_db : float or array-like
        Magnitude of S21 in dB.
    """
    return np.abs(s21_z)

def calculate_Qr(Qi, Qc):
    """
    Calculate loaded quality factor Qr from Qi and Qc.

    Parameters:
    - Qi : float or array-like
        Internal quality factor
    - Qc : float or array-like
        Coupling quality factor

    Returns:
    - Qr : float or array-like
        Loaded quality factor
    """
    return 1.0 / (1.0 / Qi + 1.0 / Qc)

def eabs_to_dnqp(e_abs, delta, vol):
    """
    Convert absorbed energy to quasiparticle density n_qp.

    Parameters:
    - e_abs : float or array
        Absorbed energy [J]
    - delta : float
        Superconducting gap energy [J]
    - vol : float
        Absorber volume [m^3]

    Returns:
    - n_qp : float or array
        Quasiparticle density [1/m^3]
    """
    return e_abs / (delta * vol)

def eabs_to_qi(e_abs, T, f0, Delta0, alpha_gamma, N_0, vol, qi0):
    """
    Convert absorbed energy to inverse Qi = 1/Qi

    Parameters:
    - e_abs : float or array
        Absorbed energy [J]
    - T : float
        Effective temperature [K]
    - f0 : float
        Resonance frequency [Hz]
    - Delta0 : float
        Superconducting gap energy [J]
    - alpha_Q : float
        Kinetic inductance fraction
    - N_0 : float
        Single-spin DoS at Fermi level [1/(eV·m³)]
    - Qi0 : float
        Baseline internal quality factor
    - vol : float
        Inductor volume [m³]

    Returns:
    - 1/Qi : float or array
    """
    dnqp = eabs_to_dnqp(e_abs, Delta0, vol)  # 1/m^3
    k1 = kappa_1(T, f0, Delta0, N_0)       # m^3
    return 1/ (alpha_gamma * k1 * dnqp + 1/qi0)

def eabs_to_fr(e_abs, T, f0, Delta0, alpha_gamma, N_0, vol):
    """
    Convert absorbed energy to shifted resonance frequency f_r.

    Parameters:
    - e_abs : float or array
        Absorbed energy [J]
    - T : float
        Effective temperature [K]
    - f0 : float
        Baseline resonance frequency [Hz]
    - Delta0 : float
        Superconducting gap energy [J]
    - alpha_gamma : float
        Combined kinetic inductance fraction and gamma
    - N_0 : float
        Single-spin DoS at Fermi level [1/(eV·m³)]
    - vol : float
        Inductor volume [m³]
    - kappa2_func : function
        Function to compute kappa_2(T, f0, Delta0, N_0)
    - base_nqp : float or None
        Baseline quasiparticle density at base temperature

    Returns:
    - fr : float or array
        Shifted resonance frequency [Hz]
    """
    dnqp = eabs_to_dnqp(e_abs, Delta0, vol)  # [1/m³]
    k2 = kappa_2(T, f0, Delta0, N_0)   # [m³]
    return f0 * (-1/2 * alpha_gamma * k2 * dnqp) + f0

def eabs_to_dS21(e_abs, T, f0, Delta0, alpha_gamma, N_0, vol, qi0, qc):
    """
    Compute change in S21 due to absorbed energy.

    Parameters:
    - e_abs : float or array
        Absorbed energy [J]
    - f : float or ndarray
        Probe frequency [Hz]
    - T : float
        Effective temperature [K]
    - f0 : float
        Baseline resonance frequency [Hz]
    - Delta0 : float
        Superconducting gap energy [J]
    - alpha_gamma : float
        Combined kinetic inductance fraction and gamma
    - N_0 : float
        DoS at Fermi level [1/(eV·m³)]
    - vol : float
        Inductor volume [m³]
    - qi0 : float
        Baseline Qi
    - qc : float
        Coupling quality factor

    Returns:
    - dS21 : complex
        Complex shift in transmission S21
    """
    # Compute dn_qp from absorbed energy
    dnqp = eabs_to_dnqp(e_abs, Delta0, vol)

    # Get kappa values at T
    k1 = kappa_1(T, f0, Delta0, N_0)
    k2 = kappa_2(T, f0, Delta0, N_0)

    # New Qi and fr
    qi = eabs_to_qi(e_abs, T, f0, Delta0, alpha_gamma, N_0, vol, qi0)
    qr = 1 / (1/qi + 1/qc)
    ds21 = qr**2/qc*alpha_gamma*(k1+ 1j*k2)*dnqp

    return ds21

def s21_ideal_eabs(f, e_abs, T, f0, Delta0, alpha_gamma, N_0, vol, qi0, qc):
    """
    Ideal S21 transmission function (complex-valued).
    
    Parameters:
    - f : ndarray or float
        Frequency [Hz]
    - fr : float
        Resonance frequency [Hz]
    - Qr : float
        Loaded quality factor
    - Qc : float
        Coupling quality factor

    Returns:
    - S21 : complex ndarray or float
        Complex transmission S21(f)
    """
    fr = eabs_to_fr(e_abs, T, f0, Delta0, alpha_gamma, N_0, vol)
    qi = eabs_to_qi(e_abs, T, f0, Delta0, alpha_gamma, N_0, vol, qi0)
    qr = 1 / (1 / qi + 1 / qc)
    x = (f - fr) / fr
    s21 = 1 - (qr / qc) / (1 + 2j * qr * x)
    return s21 

def delta_to_tc(delta_eV):
    """
    Convert energy gap Δ (eV) to critical temperature Tc (K)
    using BCS relation: Tc = Δ / ((3.528 / 2) * kB)
    """
    kB_eV = 8.617333262145e-5  # eV/K
    return delta_eV / ((3.528 / 2) * kB_eV)

def R_const(delta, n0, z1_zero=1.43, b_mat=317):
    """
    wilson 
    https://journals.aps.org/prb/pdf/10.1103/PhysRevB.69.094524 
    """
    tc = delta_to_tc(delta)
    tau_0 = z1_zero*hbar_ev/(2*math.pi*b_mat)/(kB*tc)**3
    r = 4*delta**2 / (kB*tc)**3 / (n0*tau_0)
    
    return r

def tc_to_delta(tc_K):
    """
    Convert critical temperature Tc (K) to energy gap Δ (eV)
    using BCS relation: Δ = (3.528 / 2) * kB * Tc
    """
    kB_eV = 8.617333262145e-5  # eV/K
    return (3.528 / 2) * kB_eV * tc_K

def signed_log10(x):
    return np.log10(np.abs(x)) * x/np.abs(x)

## Mazin thesis, equation 2.3: The density of thermally-excited quasiparticles
def n_qp(T, Delta0, N_0):
    # [K, eV]
    return 2.*N_0*np.sqrt(2.*np.pi*Boltz_k*T*Delta0)*np.exp(-1.*Delta0/(Boltz_k*T))

def tau_R(T, Delta, z1_zero=1.43, b_mat=317):
    """
    Compute the quasiparticle recombination time tau_R(T) [seconds].

    Parameters:
    - T     : Temperature in Kelvin
    - Delta : Superconducting gap energy in eV (Δ)
    - Tc    : Critical temperature in Kelvin (Tc)
    - tau0  : Material constant time in seconds

    Returns:
    - tau_R : Recombination time in seconds
    - eq: https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.106.167004 
    z1_zero = 1.43 # renormalization factor 
    # al_b = 0.317*1e-3 # mev**-2 
    al_b = 317 # ev**-2 
    """
    tc = delta_to_tc(Delta)
    tau_0 = z1_zero*hbar_ev/(2*math.pi*b_mat)/(kB*tc)**3

    prefactor = tau_0 / np.sqrt(np.pi)
    factor = (kB * tc / (2 * Delta))**2.5
    sqrt_term = np.sqrt(tc / T)
    expo = np.exp(Delta / (kB * T))
    return prefactor * factor * sqrt_term * expo

## Siegel thesis, equation 2.43
def kappa_1(T, f0, Delta0, N_0):
    """
    input: T (K), f0 (Hz), Delta0 (eV)
    output units: (m3) 
    """
    xi = 1./2.*(Planck_h*f0)/(Boltz_k*T)
    # return (1/(np.pi*N_0))*np.sqrt(2./(np.pi*Boltz_k*T*Delta0))*np.sinh(xi)*spec.k0(xi)
    return (1/(np.pi*Delta0*N_0))*np.sqrt((2.*Delta0)/(np.pi*Boltz_k*T))*np.sinh(xi)*spec.k0(xi)

## Siegel thesis, equation 2.44
def kappa_2(T, f0, Delta0, N_0):
    """
    input: T (K), f0 (Hz), Delta0 (eV)
    output units: (m3) 
    """
    xi = 1./2.*(Planck_h*f0)/(Boltz_k*T) #unitless
    return (1/(2.*Delta0*N_0))*(1.+np.sqrt((2.*Delta0)/(np.pi*Boltz_k*T))*np.exp(-1.*xi)*spec.i0(xi)) #m^3

## Siegel thesis, equation 2.59
def f_T(T, f0, Delta0, alpha_f, N_0, min_T=False):
    # [K, Hz, eV, _]
    # xi = 1./2.*(Planck_h*f0)/(Boltz_k*T)
    # return -1.*alpha_f/(4.*Delta0*N_0) * ( 1. + np.sqrt((2.*Delta0)/(np.pi*Boltz_k*T)) * np.exp(-1.*xi) * spec.i0(xi) ) * n_qp(T,Delta0) * f0 + f0
    if not min_T: 
        f_t = f0 * (-0.5* alpha_f * kappa_2(T, f0, Delta0, N_0) * (n_qp(T,Delta0,N_0) - n_qp(np.min(T),Delta0,N_0))) + f0
    else: 
        f_t = f0 * (-0.5* alpha_f * kappa_2(T, f0, Delta0, N_0) * (n_qp(T,Delta0,N_0) - n_qp(min_T,Delta0,N_0))) + f0
    return f_t

## Siegel thesis, equation 2.60
def Qi_T(T, f0, Delta0, alpha_Q, N_0, Qi0, min_T=False):
    # xi = 1./2.*(Planck_h*f0)/(Boltz_k*T)
    # return ( alpha_Q/(np.pi*N_0) * np.sqrt(2./(np.pi*Boltz_k*T*Delta0)) * np.sinh(xi) * spec.k0(xi) * n_qp(T,Delta0) + 1./Qi0 )**-1.
    Qi0_exp = 1/ (alpha_Q * kappa_1(np.min(T), f0, Delta0, N_0) * n_qp(np.min(T),Delta0,N_0))
    if not min_T: 
        Qi_t = 1./( alpha_Q * kappa_1(T, f0, Delta0, N_0) * (n_qp(T,Delta0,N_0) - n_qp(np.min(T),Delta0,N_0)) + 1./Qi0)
    else: 
        Qi_t = 1./( alpha_Q * kappa_1(T, f0, Delta0, N_0) * (n_qp(T,Delta0,N_0) - n_qp(min_T,Delta0,N_0)) + 1./Qi0)
    print(Qi0_exp)
    return Qi_t

def Qr_T(T, f0, Qi0, Delta0, alpha_Q):
    Qi = Qi_T(T, f0, Qi0, Delta0, alpha_Q)
    Qc = Qc_T(T, f0, Qi0, Delta0, alpha_Q)
    return (1./Qi) + (1./Qc)

## 2 * Delta * N0 * k1
def S_1(fr,T,Delta):
    # [Hz, K, eV]
    xi = 1./2.*(Planck_h*fr)/(Boltz_k*T)
    return (2/np.pi)*np.sqrt(2*Delta/(np.pi*Boltz_k*T))*np.sinh(xi)*spec.k0(xi) # unitless

## 2 * Delta * N0 * k2
def S_2(fr,T,Delta):
    # [Hz, K, eV]
    xi = 1./2.*(Planck_h*fr)/(Boltz_k*T)
    return 1+np.sqrt(2*Delta/(np.pi*Boltz_k*T))*np.exp(-1*xi)*spec.i0(xi) # unitless

## Fits to Qi, all parameters free
def MB_fitter(T_fit, Qi_fit, f_fit, fixed_alpha=False, fixed_delta=False):

    ## Define the chi-squared expression
    def chisq(f0, Delta0, alpha, Qi0):
        ## Variances of input (f,Qi) vs T values to fit
        var_Qi = np.var(Qi_fit)
        var_f  = np.var(f_fit)

        ## First term in x^2 expression
        x2_t1 = (Qi_T(T_fit, f0, Qi0, Delta0, alpha) - Qi_fit)**2./var_Qi

        ## Second term in x^2 expression
        x2_t2 = (f_T(T_fit, f0, Delta0, alpha) - f_fit)**2./var_f

        return sum( x2_t1 +  x2_t2 )
        #return sum((f_T(T_fit, f0, Delta0, alpha_f) - f_fit)**2./var_f )

    ## Initialize parameters with a guess
    f0_in     = f_fit[0]  ## Hz
    Delta0_in = 0.17e-3   ## eV
    alpha_in  = 0.03801   ## frac
    Qi0_in    = Qi_fit[0]

    ## Do the minimization problem for 500 iterations
    for j in range(500):
        minimizer = iminuit.Minuit(chisq, 
            f0=f0_in, Delta0=Delta0_in, alpha=alpha_in, Qi0=Qi0_in) 
            #limit_f0     = (f_fit[0]/1.1,f_fit[0]*1.1), 
            #limit_Delta0 = (Delta0_in,Delta0_in) if fixed_delta else (1.2e-4,2.2e-4), 
            #limit_alpha  = (alpha_in ,alpha_in ) if fixed_delta else (0.002,0.05), 
            #limit_Qi0    = (1.e2,1.e7), 
            #pedantic=False, print_level=-1)
        minimizer.limits['f0']=(f_fit[0]/1.1,f_fit[0]*1.1)
        minimizer.limits['Delta0']=(Delta0_in,Delta0_in) if fixed_delta else (1.2e-4,2.2e-4)
        minimizer.limits['alpha']=(alpha_in ,alpha_in ) if fixed_delta else (0.002,0.05)
        minimizer.limits['Qi0']=(1.e2,1.e7)
        f0_in     = minimizer.values["f0"]
        Delta0_in = minimizer.values["Delta0"]
        alpha_in  = minimizer.values["alpha"]
        Qi0_in    = minimizer.values["Qi0"]

        minimizer.migrad()

    ## Extract the final values from the minimization problem
    f0     = minimizer.values["f0"]
    Delta0 = minimizer.values["Delta0"]
    alpha  = minimizer.values["alpha"]
    Qi0    = minimizer.values["Qi0"]

    ## Get the degrees of freedom and reduced chisq
    ndof   = 4.0
    if (fixed_alpha):
        ndof -= 1.0
    if (fixed_delta):
        ndof -= 1.0

    chi_sq_dof = chisq(f0, Delta0, alpha, Qi0)/ndof

    ## F(T=0) [GHz] ; Delta(T=0) [meV] ; alpha(T=0) [frac.] ; Qr(T=0) ; reduced x2
    return f0/1.e9, Delta0*1000., alpha, Qi0, chi_sq_dof

## Fits to Qr rather than Qi
def MB_fitter_Qr(Tvals_K, Qr_fit, Fr_fit):
    fit_result = []

    def chisq(f0, D0, a, Qr0):
        ## Variances of input (f,Qi) vs T values to fit
        var_Qr = np.var(Qr_fit)
        var_fr = np.var(Fr_fit)

        ## First term in the chisq expression
        x2_t1 = np.power(Qr_T(Tvals_K, f0, Qr0, D0, a) - Qr_fit, 2) / var_Qr

        ## Second term in the chisq expression
        x2_t2 = np.power(f_T( Tvals_K, f0,      D0, a) - Fr_fit, 2) / var_fr 

        ## Calculate and return chi-squared
        x2 = np.sum( x2_t1 + x2_t2 )
        return x2

    ## Initialize parameters with a guess
    f0_in  = Fr_fit[0]
    D0_in  = 0.17e-3 ## eV
    a_in   = 0.03801 ## frac
    Qr0_in = Qr_fit[0]

    ## Do the minimization problem for 500 iterations
    for j in range(500):
        minimizer = iminuit.Minuit(chisq, 
            f0=f0_in, D0=D0_in, a=alpha_in, Qr0=Qr0_in, 
            limit_f0  = (Fr_fit[0]/1.1,Fr_fit[0]*1.1), 
            limit_D0  = (1.2e-4,2.2e-4), 
            limit_a   = (0.002,0.05), 
            limit_Qr0 = (1.e2,1.e7), 
            pedantic=False, print_level=-1)

        f0_in  = minimizer.values["f0"]
        D0_in  = minimizer.values["D0"]
        a_in   = minimizer.values["a"]
        Qr0_in = minimizer.values["Qr0"]

        minimizer.migrad()

    ## Extract the final values from the minimization problem
    f0_out  = minimizer.values["f0"]
    D0_out  = minimizer.values["D0"]
    a_out   = minimizer.values["a"]
    Qr0_out = minimizer.values["Qr0"]
    red_x2  = chisq(f0_out, D0_out, a_out, Qr0_out)/4.

    ## F(T=0) [GHz] ; Delta(T=0) [meV] ; alpha(T=0) [frac.] ; Qr(T=0) ; reduced x2
    return f0_out/1.e9, D0_out*1000., a_out, Qr0_out, red_x2

#appends caltech data for f0 to higher temp vals. doesn't really work, too much of a jump
def MB_fitter2(T_fit, Qi_fit, f_fit,added_points=11):
    fit_result = []
    added_points=10
    def chisq_f(f0, Delta0, alpha):
        
        alpha_f = alpha

        var_f = np.var(f_fit)

        #return sum( (Qi_T(T_fit[0:(len(T_fit)-added_points)], f0, Qi0, Delta0, alpha_Q) - Qi_fit)**2./var_Qi + (f_T(T_fit, f0, Delta0, alpha_f) - f_fit)**2./var_f )
        return sum((f_T(T_fit, f0, Delta0, alpha) - f_fit)**2./var_f )

    def fit_chisq_test(T_fit, f_fit, Qi_fit, f0, Delta0, alpha, Qi0):
        var_Qi = np.var(Qi_fit)
        var_f = np.var(f_fit)

        return sum(( (Qi_T(T_fit[0:(len(T_fit)-added_points)], f0, Qi0, Delta0, alpha) - Qi_fit)**2./var_Qi) +sum((f_T(T_fit, f0, Delta0, alpha) - f_fit)**2./var_f))/4.

    f0_in = f_fit[0]
    Delta0_in = 4.e-4
    alpha_in = 0.03801
    #Qi0_in = Qi_fit[0]

    for j in range(100):
        minimizer = iminuit.Minuit(chisq_f, f0=f0_in, Delta0=Delta0_in, alpha=alpha_in, limit_f0=(f_fit[0]/1.1,f_fit[0]*1.1), limit_Delta0=(1.e-4,1.e-3), limit_alpha=(0.,0.5), pedantic=False, print_level=-1)

        f0_in = minimizer.values["f0"]
        Delta0_in = minimizer.values["Delta0"]
        alpha_in = minimizer.values["alpha"]
        #Qi0_in =minimizer.values["Qi0"]

        minimizer.migrad()

    f0 = minimizer.values["f0"]
    Delta0 = minimizer.values["Delta0"]
    alpha = minimizer.values["alpha"]
    #Qi0 =minimizer.values["Qi0"]

    def chisq_qi(Qi0):
        var_qi=np.var(Qi_fit)
        return sum((Qi_T(T_fit[0:(len(T_fit)-added_points)], f0, Qi0, Delta0, alpha) - Qi_fit)**2./var_qi )

    Qi0_in = Qi_fit[0]

    for j in range(100):
        minimizer = iminuit.Minuit(chisq_qi, Qi0=Qi0_in, limit_Qi0=(1e2,1e7), pedantic=False, print_level=-1)

        Qi0_in =minimizer.values["Qi0"]

        minimizer.migrad()

    Qi0=minimizer.values["Qi0"]
    chi_sq_dof = fit_chisq_test(T_fit, f_fit, Qi_fit, f0, Delta0, alpha, Qi0)

    fit_result.append([f0/1.e9,Delta0*1000,alpha,Qi0,chi_sq_dof])

    T_smooth = np.linspace(T_fit[0],T_fit[-1],10000)

    return f0/1.e9, Delta0*1000., alpha, Qi0, chi_sq_dof

def s21_ideal(f, temp, f0, delta0, alpha_gamma, n_0, qi0, qc, min_T):
    """
    Ideal S21 transmission function (complex-valued).
    
    Parameters:
    - f : ndarray or float
        Frequency [Hz]
    - fr : float
        Resonance frequency [Hz]
    - Qr : float
        Loaded quality factor
    - Qc : float
        Coupling quality factor

    Returns:
    - S21 : complex ndarray or float
        Complex transmission S21(f)
    """
    fr = f_T(temp, f0, delta0, alpha_gamma, n_0, min_T)
    qi = Qi_T(temp, f0, delta0, alpha_gamma, n_0, qi0, min_T)
    qr = 1 / (1 / qi + 1 / qc)
    x = (f - fr) / fr
    return 1 - (qr / qc) / (1 + 2j * qr * x)

def s21_circle_radius(e_abs, T, f0, Delta0, alpha_gamma, N_0, vol, qi0, qc):
    """
    Compute the radius of the S21 resonance circle in the IQ-plane.

    Parameters
    ----------
    e_abs : ndarray or float
        Absorbed energy [eV]
    T : float
        Effective temperature [K]
    f0 : float
        Nominal resonance frequency [Hz]
    Delta0 : float
        Superconducting gap [eV]
    alpha_gamma : float
        Kinetic inductance fraction * geometry factor
    N_0 : float
        Single-spin density of states at Fermi level [1/eV μm^3]
    vol : float
        Active volume [μm^3]
    qi0 : float
        Nominal internal quality factor (baseline)
    qc : float
        Coupling quality factor

    Returns
    -------
    radius : float or ndarray
        Radius of the resonance circle in the IQ-plane.
    """
    # energy-dependent resonance parameters
    qi = eabs_to_qi(e_abs, T, f0, Delta0, alpha_gamma, N_0, vol, qi0)

    # loaded Q
    qr = 1 / (1 / qi + 1 / qc)

    # circle radius formula
    radius = qr / (2 * qc)

    # circle radius formula
    xc = 1 - qr / (2 * qc)
    return radius, xc 

