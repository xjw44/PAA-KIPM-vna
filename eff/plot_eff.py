import math 
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator, LogFormatter
from matplotlib import rcParams
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# integrate diffusion
from scipy import integrate
import numpy as np
from decimal import Decimal

# plot tau_s 
from scipy.constants import hbar, electron_volt
import math
from scipy.constants import Boltzmann, e

import sys
sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/eff')
from effEq import *
from config.eff_config import *

rcParams.update({'font.size': 14})

# xi_nb_al_ratio=0.1292
# r_e_loss_xi_al=118000*1e-18
xi_nb_al_ratio=0.09531
r_e_loss_xi_al=413000*1e-18

def plot_ph_eff_vol(plot_dir):
    # Compute recombination constants
    # Loop through devices
    plt.figure(figsize=(8, 6))

    for idx, (dev, vals) in enumerate(vol_kid.items()):
        v_al = vals["v_al"]
        v_ind = vals["v_al_act"]
        v_al_list = np.linspace(0, v_al, 1000)  # hz
        v_nb = vals.get("v_nb", 0.0)
        eta_list = phonon_collection_efficiency_vol(v_al_list, v_al, v_nb, xi_nb_al_ratio, r_e_loss_xi_al)
        eta_pred = phonon_collection_efficiency_vol(v_ind, v_al, v_nb, xi_nb_al_ratio, r_e_loss_xi_al)
        long_label = (rf"$\xi_{{Nb}}/\xi_{{Al}} = {xi_nb_al_ratio}$"+'\n'+
                      rf'$E_{{mount}}/\xi_{{Al}} = {r_e_loss_xi_al*1e18:.3e}\mathrm{{\,\mu m^3}}$'+'\n'+
                      rf"{dev}@$V_{{Al,tot}} = {v_al*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$V_{{Nb}} = {v_nb*1e18} \mathrm{{\,\mu m^3}}$")
        short_label = (
                      rf"{dev}@$V_{{Al,tot}} = {v_al*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$V_{{Nb}} = {v_nb*1e18:.3e} \mathrm{{\,\mu m^3}}$")
        if idx ==0:
            plt.plot(v_al_list*1e18, eta_list*100*0.46, label=long_label)
        else:
            plt.plot(v_al_list*1e18, eta_list*100*0.46, label=short_label)

        short_label = (rf"{eta_pred*100*0.46:.1f}%@$V_{{ind}} = {v_ind*1e18:.3e} \mathrm{{\,\mu m^3}}$")
        plt.plot(v_ind*1e18, eta_pred*100*0.46, label=short_label, marker='o', markersize=8, linestyle='None')
        short_label = (rf"{eta_pred*100*0.46:.1f}%@$V_{{ind}} = {v_ind*1e18:.3e} \mathrm{{\,\mu m^3}}$")
        plt.plot(v_ind*1e18, {}*100*0.46, label=short_label, marker='o', markersize=8, linestyle='None')

    # Plot
    plt.xscale('log')
    plt.gca().xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    plt.yscale('log')
    plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    plt.xlabel(r'$V_{\mathrm{Al}}$ [$\mu$m$^3$]')
    plt.ylabel(r'$\eta_{\mathrm{ph}}$ (%)')
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scientific y-axis
    plt.title('Phonon Collection Efficiency vs Al Volume')
    
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()

    # Save figure
    save_dir = os.path.dirname(plot_dir)
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir+".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir+".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True

