import numpy as np
import iminuit
import scipy.special as spec

# plot tau_s 
from scipy.constants import hbar, electron_volt
import math
from scipy.constants import Boltzmann, e

def phonon_collection_efficiency_vol(V_al_act, V_al_tot, V_nb, xi_nb_al_ratio, r_e_loss_xi_al):
    """
    Calculate phonon collection efficiency η_ph.

    Parameters:
    - V_ind : float
        Inductor volume [m³]
    - V_Al : float
        Aluminum volume [m³]
    - V_Nb : float
        Niobium volume [m³]
    r_e_loss_xi_al m**3

    Returns:
    - eta_ph : float
        Phonon collection efficiency (dimensionless, between 0 and 1)
    """
    denominator = V_al_tot + xi_nb_al_ratio * V_nb + r_e_loss_xi_al
    eta_ph = V_al_act / denominator
    return eta_ph

def phonon_collection_efficiency_area(f_al_act, a_tot, t_al_act, v_nb, xi_nb_al_ratio, r_e_loss_xi_al):
    """
    Calculate phonon collection efficiency η_ph using active Al area fraction.

    Parameters:
    - f_al_act : float
        Fraction of Al area that is active (unitless)
    - a_tot : float
        Total device area [m²]
    - t_al_act : float
        Active Al thickness [m]
    - v_al : float
        Total Al volume [m³]
    - v_nb : float
        Nb volume [m³]
    - xi_nb_al_ratio : float
        ξ_Nb / ξ_Al (unitless)
    - r_e_loss_xi_al : float
        Effective energy loss volume in Al units [m³]

    Returns:
    - eta_ph : float
        Phonon collection efficiency (dimensionless, between 0 and 1)
    """
    numerator = f_al_act * a_tot * t_al_act
    denominator = numerator + xi_nb_al_ratio * v_nb + r_e_loss_xi_al
    eta_ph = numerator / denominator
    return eta_ph

def phonon_eff_lifetime(tau_life, tau_collect):
    """
    Compute phonon collection efficiency.

    Parameters:
    tau_life (float): Phonon lifetime [s]
    tau_collect (float): Phonon collection time constant [s]

    Returns:
    float: Collection efficiency (unitless, between 0 and 1)
    """
    return tau_life / (tau_collect + tau_life)

def tau_collect_eq(eta, f_abs, n_abs, c_s):
    """
    Compute phonon collection time constant.

    Parameters:
    eta (float): Collection efficiency (unitless)
    f_abs (float): Fraction of phonons absorbed (unitless)
    n_abs (float): Absorber number density [1/m³]
    c_s (float): Speed of sound in material [m/s]

    Returns:
    float: Collection time constant [s]
    """
    return (4 * eta) / (f_abs * n_abs * c_s)

def infer_tau_life(eff, tau_collect):
    """
    Infer phonon lifetime from measured efficiency and tau_collect.

    Parameters:
    eff (float): Phonon collection efficiency (between 0 and 1)
    tau_collect (float): Phonon collection time constant [s]

    Returns:
    float: Phonon lifetime [s]
    """
    if eff >= 1.0:
        raise ValueError("Efficiency must be < 1.0 to compute finite tau_life.")
    if eff <= 0:
        raise ValueError("Efficiency must be > 0.")
    return (eff / (1 - eff)) * tau_collect

