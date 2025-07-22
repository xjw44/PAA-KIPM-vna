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
sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/mb')
from mbEquations import *
from config.mb_plot_config import *

rcParams.update({'font.size': 14})

N_0_al = 1.72E28 # 1/(m^3*eV), Single-spin density of states (aluminum, from Jiansong's Thesis)
N_0_hf = 3.6*1e28 # 1/(m^3*eV), Single-spin density of states (aluminum, from Jiansong's Thesis)
delta_0_hf = 38*1e-6 # ev
delta_0_al = 180*1e-6 # ev
tc_hf = delta_to_tc(delta_0_hf)
tc_al = delta_to_tc(delta_0_al)

z1_al = 1.43 # renormalization factor 
# al_b = 0.317*1e-3 # mev**-2 
al_b = 317 # ev**-2 

t_eff_al = 170*1e-3 # k 
t_eff_hf = 38*1e-3 # k 

alpha_nom = 0.038
gamma_nom = -1 # 
alpha_paa = 0.46
alpha_gamma_al = alpha_nom*abs(gamma_nom)
alpha_gamma_paa = alpha_paa*abs(gamma_nom)

Qi0_nom = 1*1e6
qc0_nom = 300*1e3
qr0_nom = calculate_Qr(Qi0_nom, qc0_nom)
target_real = 1 - qr0_nom / qc0_nom

vol_al = 11900*1e-18 # m**3
vol_hf = 21*1e-18 # m**3

nqp_target = 20*1e-18 # m**3

f_al = np.linspace(4*1e9*(1-0.001), 4*1e9*(1+0.001), 10000)  # hz
f_hf = np.linspace(4*1e9-0.005*1e9, 4*1e9+0.005*1e9, 10000)  # hz


def plot_R_const(plot_dir):
    # Compute recombination constants
    r_al = R_const(delta_0_al, N_0_al, z1_zero=z1_al, b_mat=al_b)
    r_hf = R_const(delta_0_hf, N_0_hf, z1_zero=z1_al, b_mat=al_b)
    #     print(r_al, r_hf)

    plt.figure(figsize=(8, 6))
    # Plot each bar separately with labels
    label_al = {rf'Al@R = {r_al*1e18:.2f} $\mathrm{{\mu m/s}}$'+'\n'+
               rf'$Z_1$(0) = {z1_al}'+'\n'+
               rf'b = {al_b} $\mathrm{{eV^{{-2}}}}$'+'\n'+
               rf"$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
               rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"}
    label_hf = {rf'Hf@R = {r_hf*1e18:.2f} $\mathrm{{\mu m/s}}$'+'\n'+
               rf'@$Z_1$(0) = {z1_al}'+'\n'+
               rf'b = {al_b} $\mathrm{{eV^{{-2}}}}$'+'\n'+
               rf"$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
               rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"}
    plt.bar('Al', r_al*1e18, color='gray', label=label_al)
    plt.bar('Hf', r_hf*1e18, color='blue', label=label_hf)
    plt.ylabel(r'$R$ ($\mathrm{{\mu m/s}}$)')
    plt.title('Quasiparticle Recombination Constant $R$')
    plt.grid(True)
    plt.yscale('log')

    # Minor ticks at 2–9 in each decade
    plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    # Save figure
    save_dir = os.path.dirname(plot_dir)
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.legend(loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig(plot_dir+".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir+".png", dpi=300, bbox_inches='tight')
    plt.close()

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

    plt.xlim(0, 300)
    plt.axvline(tc_hf*1e3, color='blue', linestyle=':', label=rf'$T_c^{{Hf}}$ = {tc_hf*1e3:.0f} mK')
    plt.axvline(tc_al*1e3, color='gray', linestyle=':', label=rf'$T_c^{{Al}}$ = {tc_al*1e3:.0f} mK')

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
    plt.ylabel(r'$n_{qp,\,0}$ ($\mathrm{{\mu m^{-3}}}$)')
    plt.title(r'expected $n_{qp,\,0}$')
    plt.grid(True)
    plt.yscale('log')

    # Minor ticks at 2–9 in each decade
    plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    
    # Save figure
    save_dir = os.path.dirname(plot_dir)
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.legend(loc='lower right', frameon=True)
    plt.tight_layout()
    plt.savefig(plot_dir+".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir+".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True

def plot_tau_r(plot_dir):
    plt.figure(figsize=(8, 6))
    temp = np.linspace(30, 300, 1000)  # mK

    # Compute recombination time
    tau_r_hf = tau_R(temp*1e-3, delta_0_hf, z1_zero=z1_al, b_mat=al_b)
    tau_r_al = tau_R(temp*1e-3, delta_0_al, z1_zero=z1_al, b_mat=al_b)

    long_label_al = (rf"Al@$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$" + '\n' + 
                     rf'$Z_1$(0) = {z1_al}'+'\n'+
                     rf'b = {al_b} $\mathrm{{eV^{{-2}}}}$'+'\n')
    long_label_hf = (rf"Hf@$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$" + '\n' + \
                     rf'$Z_1$(0) = {z1_al}'+'\n'+
                     rf'b = {al_b} $\mathrm{{eV^{{-2}}}}$'+'\n')

    plt.plot(temp, tau_r_hf * 1e3, label=long_label_hf)  # Convert to μs
    plt.plot(temp, tau_r_al * 1e3, label=long_label_al)

    # Find index of closest temp value
    idx_target_al = np.argmin(np.abs(temp - t_eff_al*1e3)) # mk 
    idx_target_hf = np.argmin(np.abs(temp - t_eff_hf*1e3))

    # Get tau_r at T_target for both materials (convert to ms)
    tau_hf_at_T = tau_r_hf[idx_target_hf] * 1e3
    tau_al_at_T = tau_r_al[idx_target_al] * 1e3

    # Plot vertical line at T_target
    label_al = (rf'$Al@\tau_r = {tau_al_at_T:0.2f}\,\mathrm{{ms}}$'+'\n'+
                rf'$T_{{eff}}^{{Al}}$ = {t_eff_al*1e3} mK')
    label_hf = (rf'$Hf@\tau_r = {tau_hf_at_T:0.2f}\,\mathrm{{ms}}$'+'\n'+
                rf'$T_{{eff}}^{{Hf}}$ = {t_eff_hf*1e3} mK')
    plt.axvline(t_eff_al*1e3, color='gray', linestyle='--')
    plt.axvline(t_eff_hf*1e3, color='gray', linestyle='--')
    plt.plot(t_eff_al*1e3, tau_al_at_T, 'bo', label=label_al)
    plt.plot(t_eff_hf*1e3, tau_hf_at_T, 'go', label=label_hf)

    plt.xlabel('Temperature (mK)')
    plt.ylabel(r'$\tau_r$ ($\mathrm{ms}$)')
    plt.title(r'Quasiparticle Recombination Time $\tau_r(T)$')
    plt.grid(True)
    plt.yscale('log')
    plt.xlim(0, 300)
    plt.axvline(tc_hf*1e3, color='blue', linestyle=':', label=rf'$T_c^{{Hf}}$ = {tc_hf*1e3:.0f} mK')
    plt.axvline(tc_al*1e3, color='gray', linestyle=':', label=rf'$T_c^{{Al}}$ = {tc_al*1e3:.0f} mK')

    # Minor ticks at 2–9 per decade
    plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    save_dir = os.path.dirname(plot_dir)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.legend(loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True

def plot_kappas(plot_dir):
    plt.figure(figsize=(8, 6))
    temp = np.linspace(5, 300, 1000)  # mK

    # Compute recombination time
    kappa_1_al = kappa_1(temp*1e-3, f_0_paa, delta_0_al, N_0_al)
    kappa_1_hf = kappa_1(temp*1e-3, f_0_paa, delta_0_hf, N_0_hf)

    kappa_2_al = kappa_2(temp*1e-3, f_0_paa, delta_0_al, N_0_al)
    kappa_2_hf = kappa_2(temp*1e-3, f_0_paa, delta_0_hf, N_0_hf)

    long_label_al = (rf"Al $\kappa_1$@$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n' + 
                     rf"$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                     rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$")
    long_label_hf = (rf"Hf $\kappa_1$@$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n' + 
                     rf"$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                     rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$")

    plt.plot(temp, kappa_1_al * 1e18, label=long_label_al)  # Convert to μs
    plt.plot(temp, kappa_2_al * 1e18, label=r'Al $\kappa_2$')  # Convert to μs
    plt.plot(temp, kappa_1_hf * 1e18, label=long_label_hf)  # Convert to μs
    plt.plot(temp, kappa_2_hf * 1e18, label=r'Hf $\kappa_2$')  # Convert to μs

    plt.xlim(0, 300)
    plt.axvline(tc_hf*1e3, color='blue', linestyle=':', label=rf'$T_c^{{Hf}}$ = {tc_hf*1e3:.0f} mK')
    plt.axvline(tc_al*1e3, color='gray', linestyle=':', label=rf'$T_c^{{Al}}$ = {tc_al*1e3:.0f} mK')

    # Define target temperature in mK
    T_target_mK = 10

    # Find closest index in the temperature array
    idx_target = np.argmin(np.abs(temp - T_target_mK))

    # Get corresponding kappa values (already in units of 1e18 * m³)
    k1_al_val = kappa_1_al[idx_target] * 1e18
    k2_al_val = kappa_2_al[idx_target] * 1e18
    k1_hf_val = kappa_1_hf[idx_target] * 1e18
    k2_hf_val = kappa_2_hf[idx_target] * 1e18

    # Plot vertical line at target temperature
    plt.axvline(T_target_mK, color='gray', linestyle='--', label=rf'$T = {T_target_mK}$ mK')

    # Optional: plot dots at the corresponding kappa values
    plt.plot(T_target_mK, k1_al_val, 'o', color='C0', label=rf'Al $\kappa_1$ = {k1_al_val:.1e} $\mathrm{{\mu m^3}}$')
    plt.plot(T_target_mK, k2_al_val, 's', color='C1', label=rf'Al $\kappa_2$ = {k2_al_val:.1e} $\mathrm{{\mu m^3}}$')
    plt.plot(T_target_mK, k1_hf_val, 'o', color='C2', label=rf'Hf $\kappa_1$ = {k1_hf_val:.1e} $\mathrm{{\mu m^3}}$')
    plt.plot(T_target_mK, k2_hf_val, 's', color='C3', label=rf'Hf $\kappa_2$ = {k2_hf_val:.1e} $\mathrm{{\mu m^3}}$')

    plt.xlabel('Temperature (mK)')
    plt.ylabel(r'$\kappa$ ($\mathrm{\mu m^3}$)')
    plt.title(r'$\kappa$ vs Temperature')
    plt.grid(True)
    plt.legend(loc='upper right', frameon=True)
    plt.tight_layout()

    save_dir = os.path.dirname(plot_dir)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True

def plot_dissipation(plot_dir):
    plt.figure(figsize=(8, 6))
    temp_al = np.linspace(t_eff_al*1e3, 500, 1000)  # mK
    temp_hf = np.linspace(t_eff_hf*1e3, 500, 1000)  # mK

    # Compute recombination time
    qi_al = Qi_T(temp_al*1e-3, f_0_paa, delta_0_al, alpha_gamma_al, N_0_al, Qi0_nom)
    qi_hf = Qi_T(temp_hf*1e-3, f_0_paa, delta_0_hf, alpha_gamma_paa, N_0_hf, Qi0_nom)

    long_label_al = (rf"Al@$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n' + 
                     rf"$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                     rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
                     rf"$\alpha = {alpha_nom}$"+'\n'+
                     rf"$\gamma = {gamma_nom}$")
    long_label_hf = (
        rf"Hf@$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$"+'\n'+ 
        rf"$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
        rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
        rf"$\alpha = {alpha_paa}$"+'\n'+
        rf"$\gamma = {gamma_nom}$"
    )

    plt.plot(temp_al, qi_al, label=long_label_al)  # Convert to μs
    plt.plot(temp_hf, qi_hf, label=long_label_hf)  # Convert to μs

    plt.xlim(0, 300)
    plt.axvline(tc_hf*1e3, color='blue', linestyle=':', label=rf'$T_c^{{Hf}}$ = {tc_hf*1e3:.0f} mK')
    plt.axvline(tc_al*1e3, color='gray', linestyle=':', label=rf'$T_c^{{Al}}$ = {tc_al*1e3:.0f} mK')

    # Plot vertical line at T_target
    label_al = (
        rf"$Q_{{i,0}} = {Qi0_nom:.2g}$"+'\n'+
        rf'$T_{{eff}}^{{Al}}$ = {t_eff_al*1e3} mK'
    )
    label_hf = (
        rf"$Q_{{i,0}} = {Qi0_nom:.2g}$"+'\n'+
        rf'$T_{{eff}}^{{Hf}}$ = {t_eff_hf*1e3} mK'
    )
    plt.axvline(t_eff_al*1e3, color='gray', linestyle='--')
    plt.axvline(t_eff_hf*1e3, color='gray', linestyle='--')
    plt.plot(t_eff_al*1e3, Qi0_nom, 'bo', label=label_al)
    plt.plot(t_eff_hf*1e3, Qi0_nom, 'go', label=label_hf)

    plt.xlabel('Temperature (mK)')
    plt.ylabel(r'$Q_i(T)$')
    plt.title(r'Internal Quality Factor $Q_i$ vs Temperature')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()

    plt.yscale('log')
    # Minor ticks at 2–9 in each decade
    plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    save_dir = os.path.dirname(plot_dir)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True

def plot_frequency(plot_dir):
    plt.figure(figsize=(8, 6))
    temp_al = np.linspace(t_eff_al*1e3, 500, 1000)  # mK
    temp_hf = np.linspace(t_eff_hf*1e3, 500, 1000)  # mK

    # Compute recombination time
    fr_al = f_T(temp_al*1e-3, f_0_paa, delta_0_al, alpha_gamma_al, N_0_al)
    fr_hf = f_T(temp_hf*1e-3, f_0_paa, delta_0_hf, alpha_gamma_paa, N_0_hf)

    plt.xlim(0, 300)
    plt.axvline(tc_hf*1e3, color='blue', linestyle=':', label=rf'$T_c^{{Hf}}$ = {tc_hf*1e3:.0f} mK')
    plt.axvline(tc_al*1e3, color='gray', linestyle=':', label=rf'$T_c^{{Al}}$ = {tc_al*1e3:.0f} mK')

    long_label_al = (rf"Al@$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                     rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
                     rf"$\alpha = {alpha_nom}$"+'\n'+
                     rf"$\gamma = {gamma_nom}$"+'\n')
    long_label_hf = (
        rf"Hf@$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
        rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
        rf"$\alpha = {alpha_paa}$"+'\n'+
        rf"$\gamma = {gamma_nom}$"+'\n'
    )

    plt.plot(temp_al, fr_al*1e-9, label=long_label_al)  # Convert to μs
    plt.plot(temp_hf, fr_hf*1e-9, label=long_label_hf)  # Convert to μs

    # Plot vertical line at T_target
    label_al = (
        rf"$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n' + 
        rf'$T_{{eff}}^{{Al}}$ = {t_eff_al*1e3} mK'
    )
    label_hf = (
        rf"$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$"+'\n'+ 
        rf'$T_{{eff}}^{{Hf}}$ = {t_eff_hf*1e3} mK'
    )
    plt.axvline(t_eff_al*1e3, color='gray', linestyle='--')
    plt.axvline(t_eff_hf*1e3, color='gray', linestyle='--')
    plt.plot(t_eff_al*1e3, f_0_paa*1e-9, 'bo', label=label_al)
    plt.plot(t_eff_hf*1e3, f_0_paa*1e-9, 'go', label=label_hf)

    plt.xlabel('Temperature (mK)')
    plt.ylabel(r'$f_r(T)\,(GHz)$')
    plt.title(r'Resonant frequency $f_r$ vs Temperature')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()

    # plt.yscale('log')
    # # Minor ticks at 2–9 in each decade
    # plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    save_dir = os.path.dirname(plot_dir)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True

def s21_colors(temp_al, temp_hf):
    num_temp_al = len(temp_al)+1
    num_temp_hf = len(temp_hf)+1

    cmap_al = plt.cm.inferno
    colors_al = cmap_al(np.linspace(0, 1, num_temp_al))

    cmap_hf = plt.cm.viridis
    colors_hf = cmap_hf(np.linspace(0, 1, num_temp_hf))

    # Al colorbar
    norm_al = mcolors.Normalize(vmin=min(temp_al), vmax=max(temp_al))
    sm_al = cm.ScalarMappable(cmap=cmap_al, norm=norm_al)
    sm_al.set_array([])

    # Hf colorbar
    norm_hf = mcolors.Normalize(vmin=min(temp_hf), vmax=max(temp_hf))
    sm_hf = cm.ScalarMappable(cmap=cmap_hf, norm=norm_hf)
    sm_hf.set_array([])

    return colors_al, colors_hf, sm_al, sm_hf

def plot_s21(plot_dir):
    fig, axs = plt.subplots(2, 2, figsize=(16, 12), sharex=False, constrained_layout=True)
    temp_al = np.linspace(t_eff_al, tc_al/3, 10)  # mK
    temp_hf = np.linspace(t_eff_hf, tc_hf/3, 10)  # mK

    colors_al, colors_hf, sm_al, sm_hf = s21_colors(temp_al, temp_hf)

    # Compute recombination time
    for ind, temp in enumerate(temp_al): 
        s21_al = s21_ideal(f_al, temp_al[ind], f_0_paa, delta_0_al, alpha_gamma_al, 
            N_0_al, Qi0_nom, qc0_nom, t_eff_al)
        s21_hf = s21_ideal(f_hf, temp_hf[ind], f_0_paa, delta_0_hf, alpha_gamma_paa, 
            N_0_hf, Qi0_nom, qc0_nom, t_eff_hf)

        if ind==0: 
            label_al = rf'$T_{{eff}}^{{Al}}$ = {t_eff_al*1e3} mK'
            axs[0, 0].plot(f_al*1e-9, 20 * np.log10(np.abs(s21_al)), color=colors_al[ind], label=label_al)
            label_hf = rf'$T_{{eff}}^{{Hf}}$ = {t_eff_hf*1e3} mK'
            axs[0, 1].plot(f_hf*1e-9, 20 * np.log10(np.abs(s21_hf)), color=colors_hf[ind], label=label_hf)
        elif ind==len(temp_al)-1: 
            label_al = (
                rf'$T_c^{{Al}}$ = {tc_al*1e3:.0f} mK'+'\n'+
                rf"$Q_{{i,0}} = {Qi0_nom:.2g}$"+'\n'+
                rf"$Q_{{c}} = {qc0_nom:.2g}$"+'\n'+
                rf"$f_{{r,0}} = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n' + 
                rf"$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
                rf"$\alpha = {alpha_nom}$"+'\n'+
                rf"$\gamma = {gamma_nom}$"
            )
            label_hf = (
                rf'$T_c^{{Hf}}$ = {tc_hf*1e3:.0f} mK'+'\n'+
                rf"$Q_{{i,0}} = {Qi0_nom:.2g}$"+'\n'+
                rf"$Q_{{c}} = {qc0_nom:.2g}$"+'\n'+
                rf"$f_{{r,0}} = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n' + 
                rf"$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
                rf"$\alpha = {alpha_paa}$"+'\n'+
                rf"$\gamma = {gamma_nom}$"
            )
            axs[0, 0].plot(f_al*1e-9, 20 * np.log10(np.abs(s21_al)), color=colors_al[ind], label=label_al)
            axs[0, 1].plot(f_hf*1e-9, 20 * np.log10(np.abs(s21_hf)), color=colors_hf[ind], label=label_hf)
        else: 
            axs[0, 0].plot(f_al*1e-9, 20 * np.log10(np.abs(s21_al)), color=colors_al[ind])
            axs[0, 1].plot(f_hf*1e-9, 20 * np.log10(np.abs(s21_hf)), color=colors_hf[ind])
        axs[1, 0].plot(np.real(s21_al), np.imag(s21_al), color=colors_al[ind])
        axs[1, 1].plot(np.real(s21_hf), np.imag(s21_hf), color=colors_hf[ind])

    cbar_al = fig.colorbar(sm_al, ax=[axs[0, 0], axs[0, 1]], orientation='vertical', fraction=0.015, pad=0.04)
    cbar_al.set_label('Al: temperature (mK)')

    cbar_hf = fig.colorbar(sm_hf, ax=[axs[1, 0], axs[1, 1]], orientation='vertical', fraction=0.015, pad=0.04)
    cbar_hf.set_label('Hf: temperature (mK)')

    axs[1, 0].set_title("Al: Re/Im S21")
    axs[1, 0].set_xlabel("Re S21")
    axs[1, 0].set_ylabel("Im S21")
    axs[1,0].axis('equal')
    axs[1,1].axis('equal')

    axs[1, 1].set_title("Hf: Re/Im S21")
    axs[1, 1].set_xlabel("Re S21")
    axs[1, 1].set_ylabel("Im S21")

    axs[0, 0].set_title("Al: S21 mag")
    axs[0, 0].set_xlabel("frequency (GHz)")
    axs[0, 0].set_ylabel("S21 (dB)")
    axs[0, 0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    axs[0, 1].set_title("Hf: S21 mag")
    axs[0, 1].set_xlabel("frequency (GHz)")
    axs[0, 1].set_ylabel("S21 (dB)")
    axs[0, 1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    # Top-left plot
    axs[0, 0].set_ylim(-1.5, 0)
    axs[0, 0].set_xlim(3.998, 4.001)
    # Top-right plot
    axs[0, 1].set_ylim(-1.5, 0)
    axs[0, 1].set_xlim(3.995, 4.001)

    for ax in axs.flat:
        ax.grid(True)

    # plt.yscale('log')
    # # Minor ticks at 2–9 in each decade
    # plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    save_dir = os.path.dirname(plot_dir)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True

def plot_energy_response(plot_dir, plot_log=False):
    n_x = 2
    n_y = 3
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x), constrained_layout=True)
    e_abs = np.linspace(1, 500, 1000)  # meV
    e_abs_al = np.linspace(0, 1000, 10)  # eV
    e_abs_hf = np.linspace(0, 500, 10)  # meV

    # Compute recombination time
    dnqp_al = eabs_to_dnqp(e_abs*1e-3, delta_0_al, vol_al)
    dnqp_hf = eabs_to_dnqp(e_abs*1e-3, delta_0_hf, vol_hf)

    qi_al = eabs_to_qi(e_abs*1e-3, t_eff_al, f_0_paa, delta_0_al, alpha_gamma_al, N_0_al, vol_al, Qi0_nom)
    qi_hf = eabs_to_qi(e_abs*1e-3, t_eff_hf, f_0_paa, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_hf, Qi0_nom)

    fr_al = eabs_to_fr(e_abs*1e-3, t_eff_al, f_0_paa, delta_0_al, alpha_gamma_al, N_0_al, vol_al)
    fr_hf = eabs_to_fr(e_abs*1e-3, t_eff_hf, f_0_paa, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_hf)

    s21_org = s21_ideal_eabs(f_0_paa, 0, t_eff_al, f_0_paa, delta_0_al, alpha_gamma_al, N_0_al, vol_al, Qi0_nom, qc0_nom)
    s21_al = s21_ideal_eabs(f_0_paa, e_abs*1e-3, t_eff_al, f_0_paa, delta_0_al, alpha_gamma_al, N_0_al, vol_al, Qi0_nom, qc0_nom)
    s21_hf = s21_ideal_eabs(f_0_paa, e_abs*1e-3, t_eff_hf, f_0_paa, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_hf, Qi0_nom, qc0_nom)
    s21_hf_volnom = s21_ideal_eabs(f_0_paa, e_abs*1e-3, t_eff_hf, f_0_paa, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_al, Qi0_nom, qc0_nom)
    s21_al_db = s21_z_to_mag(s21_al)
    s21_hf_db = s21_z_to_mag(s21_hf)
    s21_hf_db_volnom = s21_z_to_mag(s21_hf_volnom)

    # Subplot 1: n_qp vs E_abs
    indp = axs[0,0]
    label_al = (
                rf"Al@$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                rf"$V_{{ind}} = {vol_al*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"+'\n'+
                rf'$n_{{qp,0}}$ = {nqp_target*1e18:.0f} $\mathrm{{\mu m^{{-3}}}}$')
    label_hf = (
                rf"Hf@$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                rf"$V_{{ind}} = {vol_hf*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"+'\n')
    indp.plot(e_abs, dnqp_al*1e18, label=label_al)
    indp.plot(e_abs, dnqp_hf*1e18, label=label_hf)
    indp.set_xlabel(r'$E_{\mathrm{abs}}$ (meV)')
    indp.set_ylabel(r'$\delta n_{qp}$ ($\mu$m$^{-3}$)')
    indp.set_title(r'$\delta n_{qp}$ vs $E_{\mathrm{abs}}$')
    indp.legend(loc='upper right', frameon=True)
    indp.set_yscale('log')
    indp.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    # Subplot 1: n_qp vs E_abs
    indp = axs[0,1]
    label_al = (
                rf"Al@$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                rf"$V_{{ind}} = {vol_al*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"+'\n'+
                rf'$T_{{eff}}^{{Al}}$ = {t_eff_al*1e3} mK'+'\n'+
                rf"$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n'+
                rf"$Q_{{i,0}} = {Qi0_nom:.2g}$"+'\n'+
                rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
                rf"$\alpha = {alpha_nom}$"+'\n'+
                rf"$\gamma = {gamma_nom}$")
    label_hf = (
                rf"Hf@$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                rf"$V_{{ind}} = {vol_hf*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"+'\n'+
                rf'$T_{{eff}}^{{Hf}}$ = {t_eff_hf*1e3} mK'+'\n'+
                rf"$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n'+
                rf"$Q_{{i,0}} = {Qi0_nom:.2g}$"+'\n'+
                rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
                rf"$\alpha = {alpha_paa}$"+'\n'+
                rf"$\gamma = {gamma_nom}$")
    indp.plot(e_abs, qi_al, label=label_al)
    indp.plot(e_abs, qi_hf, label=label_hf)
    indp.set_xlabel(r'$E_{\mathrm{abs}}$ (meV)')
    indp.set_ylabel(r'$Q_i$')
    indp.set_title(r'$Q_i$ vs $E_{\mathrm{abs}}$')
    indp.legend(loc='lower right', frameon=True)
    indp.set_yscale('log')
    indp.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    # Subplot 1: n_qp vs E_abs
    indp = axs[0,2]
    label_al = (
                rf"Al@$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                rf"$V_{{ind}} = {vol_al*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"+'\n'+
                rf'$T_{{eff}}^{{Al}}$ = {t_eff_al*1e3} mK'+'\n'+
                rf"$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n'+
                rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
                rf"$\alpha = {alpha_nom}$"+'\n'+
                rf"$\gamma = {gamma_nom}$"+'\n')
    label_hf = (
                rf"Hf@$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                rf"$V_{{ind}} = {vol_hf*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"+'\n'+
                rf'$T_{{eff}}^{{Hf}}$ = {t_eff_hf*1e3} mK'+'\n'+
                rf"$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n'+
                rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
                rf"$\alpha = {alpha_paa}$"+'\n'+
                rf"$\gamma = {gamma_nom}$"+'\n')
    indp.plot(e_abs, fr_al*1e-9, label=label_al)
    indp.plot(e_abs, fr_hf*1e-9, label=label_hf)
    indp.set_xlabel(r'$E_{\mathrm{abs}}$ (meV)')
    indp.set_ylabel(r'$f_r$ (GHz)')
    indp.set_title(r'$f_r$ vs $E_{\mathrm{abs}}$')
    indp.legend(loc='lower right', frameon=True)
    # indp.set_yscale('log')
    # indp.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    # Subplot 1: n_qp vs E_abs
    colors_al, colors_hf, sm_al, sm_hf = s21_colors(e_abs_al, e_abs)

    indp = axs[1,0]
    label_al = (
                rf"Al@$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
                rf"$V_{{ind}} = {vol_al*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"+'\n'+
                rf'$T_{{eff}}^{{Al}}$ = {t_eff_al*1e3} mK'+'\n'+
                rf"$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n'+
                rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
                rf"$\alpha = {alpha_nom}$"+'\n'+
                rf"$\gamma = {gamma_nom}$"+'\n'+
                rf"$Q_{{i,0}} = {Qi0_nom:.2g}$"+'\n'+
                rf"$Q_{{c}} = {qc0_nom:.2g}$"+'\n'+
                rf"$Q_{{r,0}} = {qr0_nom:.2g}$"+'\n')
    for ind, temp in enumerate(e_abs_al): 
        s21_al = s21_ideal_eabs(f_al, temp, t_eff_al, f_0_paa, delta_0_al, alpha_gamma_al, N_0_al, vol_al, Qi0_nom, qc0_nom)
        s21_al_fr = s21_ideal_eabs(f_0_paa, temp, t_eff_al, f_0_paa, delta_0_al, alpha_gamma_al, N_0_al, vol_al, Qi0_nom, qc0_nom)
        indp.plot(np.real(s21_al), np.imag(s21_al), color=colors_al[ind])
        if ind==0: 
            indp.plot(np.real(s21_al_fr), np.imag(s21_al_fr), label=label_al, marker='o', markersize=8, linestyle='None')
        else: 
            indp.plot(np.real(s21_al_fr), np.imag(s21_al_fr), marker='o', markersize=8, linestyle='None')
    # indp.plot(np.real(ds21_hf), np.imag(ds21_hf), label=label_hf)
    indp.axvline(target_real, color='gray', linestyle='--', label=rf'$1 - Q_r / Q_c = {target_real:.3f}$')

    # Equal aspect ratio
    indp.set_aspect('equal', adjustable='box')
    # Axis labels and title
    indp.set_xlabel(r'Re[$S_{21}$]')
    indp.set_ylabel(r'Im[$S_{21}$]')
    indp.set_title(r'$S_{21}$ Complex Plane')
    indp.legend(loc='lower right', frameon=True)

    cbar_al = fig.colorbar(sm_al, ax=indp, orientation='vertical', fraction=0.015, pad=0.04)
    cbar_al.set_label(r'Al: $E_{abs}$ (eV)')

    # Subplot 1: n_qp vs E_abs
    indp = axs[1,1]

    label_hf = (
            rf"Hf@$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"+'\n'+
            rf"$V_{{ind}} = {vol_hf*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"+'\n'+
            rf'$T_{{eff}}^{{Hf}}$ = {t_eff_hf*1e3} mK'+'\n'+
            rf"$f_r = {f_0_paa*1e-9:.0f}\,\mathrm{{GHz}}$" + '\n'+
            rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"+'\n'+
            rf"$\alpha = {alpha_paa}$"+'\n'+
            rf"$\gamma = {gamma_nom}$"+'\n'+
            rf"$Q_{{i,0}} = {Qi0_nom:.2g}$"+'\n'+
            rf"$Q_{{c}} = {qc0_nom:.2g}$"+'\n'+
            rf"$Q_{{r,0}} = {qr0_nom:.2g}$"+'\n')

    for ind, temp in enumerate(e_abs_hf): 
        s21_hf = s21_ideal_eabs(f_hf, temp*1e-3, t_eff_hf, f_0_paa, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_hf, Qi0_nom, qc0_nom)
        s21_hf_fr = s21_ideal_eabs(f_0_paa, temp*1e-3, t_eff_hf, f_0_paa, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_hf, Qi0_nom, qc0_nom)
        indp.plot(np.real(s21_hf), np.imag(s21_hf), color=colors_hf[ind])
        if ind==0: 
            indp.plot(np.real(s21_hf_fr), np.imag(s21_hf_fr), label=label_hf, marker='o', markersize=8, linestyle='None')
        else: 
            indp.plot(np.real(s21_hf_fr), np.imag(s21_hf_fr), marker='o', markersize=8, linestyle='None')
    # indp.plot(np.real(ds21_hf), np.imag(ds21_hf), label=label_hf)
    indp.axvline(target_real, color='gray', linestyle='--', label=rf'$1 - Q_r / Q_c = {target_real:.3f}$')

    # Equal aspect ratio
    indp.set_aspect('equal', adjustable='box')
    # Axis labels and title
    indp.set_xlabel(r'Re[$S_{21}$]')
    indp.set_ylabel(r'Im[$S_{21}$]')
    indp.set_title(r'$S_{21}$ Complex Plane')
    indp.legend(loc='lower left', frameon=True)

    cbar_hf = fig.colorbar(sm_hf, ax=indp, orientation='vertical', fraction=0.015, pad=0.04)
    cbar_hf.set_label(r'Al: $E_{abs}$ (meV)')

    # Subplot 1: n_qp vs E_abs
    indp = axs[1,2]
    indp.plot(e_abs, s21_al_db, label=label_al)
    indp.plot(e_abs, s21_hf_db, label=label_hf)
    indp.plot(e_abs, s21_hf_db_volnom, label=rf'Hf with $V_{{ind}}={vol_al}\,\mathrm{{\mu m^{{3}}}}$')
    indp.set_xlabel(r'$E_{\mathrm{abs}}$ (meV)')
    indp.set_ylabel(r'$|S_{21}(f_{r,0})|$')
    indp.set_title(r'$\delta S_{21}$ vs $E_{\mathrm{abs}}$')
    indp.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    if plot_log:
        for ax in [axs[0,1], axs[0,2], axs[1,2]]: 
            ax.set_yscale('log')
            ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    for ax in axs.flat:
        ax.grid(True)

    # plt.yscale('log')
    # # Minor ticks at 2–9 in each decade
    # plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    save_dir = os.path.dirname(plot_dir)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True

