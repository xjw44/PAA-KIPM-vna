import math 
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator, LogFormatter

# integrate diffusion
from scipy import integrate
import numpy as np

# plot tau_s 
from scipy.constants import hbar, electron_volt
import math
from scipy.constants import Boltzmann, e

import sys
sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna/MB')
from mbEquations import *

N_0_al = 1.72E28 # 1/(m^3*eV), Single-spin density of states (aluminum, from Jiansong's Thesis)
N_0_hf = 3.6*1e28 # 1/(m^3*eV), Single-spin density of states (aluminum, from Jiansong's Thesis)
delta_0_hf = 38*1e-6 # ev
delta_0_al = 180*1e-6 # ev

def plotN_qp(plot_dir): 
    plt.figure(figsize=(8, 6))
    temp = np.linspace(30, 300, 1000) # mk
    nqp_t_hf = n_qp(temp*1e-3, delta_0_hf, N_0_hf)
    nqp_t_al = n_qp(temp*1e-3, delta_0_al, N_0_al)
    long_label_al = rf"Al @ $\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+\
                    rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"
    long_label_hf = rf"Hf @ $\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+\
                    rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"
    plt.plot(temp, nqp_t_hf*1e-18, label=long_label_hf)
    plt.plot(temp, nqp_t_al*1e-18, label=long_label_al)

    # Mark each tau_target with vertical & horizontal lines
    for nqp_target in [20]:
        # Find closest match for Hf
        idx_hf = np.argmin(np.abs(nqp_t_hf - nqp_target*1e18))
        temp_hf_match = temp[idx_hf]

        # Find closest match for Al
        idx_al = np.argmin(np.abs(nqp_t_al - nqp_target*1e18))
        temp_al_match = temp[idx_al]

        # Add markers to the plot
        plt.axhline(nqp_target, color='gray', linestyle='--',
            label=rf'$n_{{qp}}$ = {nqp_target:.0f} $\mathrm{{\mu m^{{-3}}}}$')
        plt.plot(temp_hf_match, nqp_t_hf[idx_hf]*1e-18, 'bo', label=rf'$T_{{Hf}}$={temp_hf_match:.1f} mK')
        plt.plot(temp_al_match, nqp_t_al[idx_al]*1e-18, 'go', label=rf'$T_{{Al}}$={temp_al_match:.1f} mK')

    plt.xlabel('temperature (mK)')
    plt.ylabel(r'$n_{qp}$ ($\mathrm{{\mu m^{-3}}}$)')
    plt.title(r'expected $n_{qp}$')
    plt.grid(True)
    plt.yscale('log')

    # Minor ticks at 2–9 in each decade
    plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    plt.grid(True)
    
    # Save figure
    save_dir = os.path.dirname(plot_dir)
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.legend(loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig(plot_dir+".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir+".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True
