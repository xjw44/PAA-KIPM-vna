import numpy as np
import iminuit
import scipy.special as spec

# plot tau_s 
from scipy.constants import hbar, electron_volt, Boltzmann, e, pi
import math
from scipy.integrate import quad

import sys
# sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/mb')
from pathlib import Path

# Path to the directory where this script lives
here = Path(__file__).resolve().parent

# Append ../res and ../mb relative to config/
sys.path.append(str(here.parent / "mb"))

from mbEquations import *

# print(f"hbar = {hbar:.5e} J·s")
# print(f"e = {e:.5e} C")  # Coulombs
# print(f"k_B = {Boltzmann:.5e} J/K")

#############################
# iMinuit fitting functions #
#############################

######################
# Mattis-Bardeen fit #
######################

nqp0_scan_list = np.linspace(0.2*1e18, 200*1e18, 1000) # m**-3

def t_eff_nqp_scan(delta_0, n_0, nqp0_list=nqp0_scan_list):
    """
    Find T_eff (in K) for each target quasiparticle density n_qp0
    such that n_qp(T_eff) ≈ n_qp0 (nearest match).

    Parameters
    ----------
    nqp0_list : array-like
        Target quasiparticle densities [1/m^3].
    delta_0 : float
        Superconducting gap (passed to n_qp).
    n_0 : float
        Single-spin DOS at E_F (passed to n_qp).
    temp_min_mK, temp_max_mK : float
        Temperature scan range in mK (converted to K internally).
    npts : int
        Number of temperature samples for building the lookup.

    Returns
    -------
    temps_K : np.ndarray
        Temperatures [K] corresponding to each target.
    """
    # Temperature grid
    temp = np.linspace(10, 100, 1000) # mk
    T_K  = temp * 1e-3

    # Evaluate n_qp(T)
    nqp_T = n_qp(T_K, delta_0, n_0)

    targets = np.asarray(nqp0_list, dtype=float)
    temps_K = np.zeros_like(targets, dtype=float)

    for i, y in enumerate(targets):
        j = int(np.argmin(np.abs(nqp_T - y)))
        temps_K[i] = T_K[j]

    return temps_K

def area_coverage(covered_area, total_area, *, debug=False):
    """
    Calculate area coverage percentage.

    Parameters
    ----------
    covered_area : float
        Covered area [m^2].
    total_area : float
        Total area [m^2].
    debug : bool, default=False
        If True, print intermediate values.

    Returns
    -------
    coverage_percent : float
        Coverage as a percentage of total area.
    """
    if total_area <= 0:
        raise ValueError("Total area must be positive.")

    coverage_percent = (covered_area / total_area) * 100.0

    if debug:
        print("=== area_coverage DEBUG ===")
        print(f"Covered area = {covered_area:.4e} m^2")
        print(f"Total area   = {total_area:.4e} m^2")
        print(f"Coverage     = {coverage_percent:.2f} %")
        print("===========================")

    return coverage_percent