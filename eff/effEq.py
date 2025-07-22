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

