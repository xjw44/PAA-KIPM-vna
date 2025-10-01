import math 
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator, LogFormatter
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.collections import PolyCollection
from matplotlib.transforms import Bbox

# integrate diffusion
from scipy import integrate
import numpy as np
from decimal import Decimal

# plot tau_s 
from scipy.constants import hbar, electron_volt
import math
from scipy.constants import Boltzmann, e

import sys
# sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/eff')

from pathlib import Path

# Path to the directory where this script lives
here = Path(__file__).resolve().parent

# Append ../res and ../mb relative to the script’s parent folder
sys.path.append(str(here.parent / "eff"))

from effEq import *
from config.eff_config import *
from config.plt_helper import *

def plot_ph_eff_vol(plot_dir):
    # Compute recombination constants
    # Loop through devices
    # plt.figure(figsize=(8, 6))
    n_x = 1
    n_y = 3
    fig, axs = plt.subplots(n_x, n_y, figsize=(10*n_y, 6*n_x), constrained_layout=True)
    ax_vol, ax_area, ax_frac = axs

    for idx, (dev, vals) in enumerate(vol_kid.items()):
        v_al = vals["v_al"]
        v_ind = vals["v_al_act"]
        eff_measured = vals.get("eff_measured", 0)
        v_al_list = np.linspace(0, v_al, 1000)  # hz
        v_nb = vals.get("v_nb", 0.0)
        eta_list = phonon_collection_efficiency_vol(v_al_list, v_al, v_nb, xi_nb_al_ratio, r_e_loss_xi_al)
        eta_pred = phonon_collection_efficiency_vol(v_ind, v_al, v_nb, xi_nb_al_ratio, r_e_loss_xi_al)
        long_label = (rf"$\xi_{{Nb}}/\xi_{{Al}} = {xi_nb_al_ratio}$"+'\n'+
                      rf'$E_{{mount}}/\xi_{{Al}} = {r_e_loss_xi_al*1e18:.3e}\mathrm{{\,\mu m^3}}$'+'\n'+
                      rf'$\eta_{{pb}} = {eta_pb*100}\%$'+'\n'+
                      rf"{dev}@$V_{{Al,tot}} = {v_al*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$V_{{Nb}} = {v_nb*1e18} \mathrm{{\,\mu m^3}}$")
        short_label = (
                      rf"{dev}@$V_{{Al,tot}} = {v_al*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$V_{{Nb}} = {v_nb*1e18:.3e} \mathrm{{\,\mu m^3}}$")
        if idx ==0:
            ax_vol.plot(v_al_list*1e18, eta_list*100*eta_pb, label=long_label)
        else:
            ax_vol.plot(v_al_list*1e18, eta_list*100*eta_pb, label=short_label)

        short_label = (rf"{eta_pred*100*eta_pb:.1f}%@$V_{{Al,act}} = {v_ind*1e18:.3e} \mathrm{{\,\mu m^3}}$")
        ax_vol.plot(v_ind*1e18, eta_pred*100*eta_pb, label=short_label, marker='o', markersize=8, linestyle='None')
        short_label = (rf"{eff_measured*100:.1f}% measured")
        ax_vol.plot(v_ind*1e18, eff_measured*100, label=short_label, marker='*', markersize=8, linestyle='None')

        t_al = vals["t_al"]
        a_al_list = v_al_list/t_al  # hz
        a_ind = v_ind/t_al
        long_label = (rf"$\xi_{{Nb}}/\xi_{{Al}} = {xi_nb_al_ratio}$"+'\n'+
                      rf'$E_{{mount}}/\xi_{{Al}} = {r_e_loss_xi_al*1e18:.3e}\mathrm{{\,\mu m^3}}$'+'\n'+
                      rf"{dev}@$V_{{Al,tot}} = {v_al*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$V_{{Nb}} = {v_nb*1e18} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$t_{{Al}} = {t_al*1e9:.2e} \mathrm{{\,nm}}$")
        short_label = (
                      rf"{dev}@$V_{{Al,tot}} = {v_al*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$V_{{Nb}} = {v_nb*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$t_{{Al}} = {t_al*1e9:.2e} \mathrm{{\,nm}}$")
        if idx ==0:
            ax_area.plot(a_al_list*1e12, eta_list*100, label=long_label)
        else:
            ax_area.plot(a_al_list*1e12, eta_list*100, label=short_label)

        short_label = (rf"{eta_pred*100:.1f}%@$A_{{Al,act}} = {a_ind*1e6:.3e} \mathrm{{\,mm^2}}$")
        ax_area.plot(a_ind*1e12, eta_pred*100, label=short_label, marker='o', markersize=8, linestyle='None')

        a_tot = vals["a_tot"]
        a_frac_list = a_al_list/a_tot
        a_ind_frac = a_ind/a_tot

        long_label = (rf"$\xi_{{Nb}}/\xi_{{Al}} = {xi_nb_al_ratio}$"+'\n'+
                      rf'$E_{{mount}}/\xi_{{Al}} = {r_e_loss_xi_al*1e18:.3e}\mathrm{{\,\mu m^3}}$'+'\n'+
                      rf"{dev}@$V_{{Al,tot}} = {v_al*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$V_{{Nb}} = {v_nb*1e18} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$t_{{Al}} = {t_al*1e9:.2e} \mathrm{{\,nm}}$"+'\n'+
                      rf"$a_{{tot}} = {a_tot*1e6:.2e} \mathrm{{\,mm^2}}$")
        short_label = (
                      rf"{dev}@$V_{{Al,tot}} = {v_al*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$V_{{Nb}} = {v_nb*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
                      rf"$t_{{Al}} = {t_al*1e9:.2e} \mathrm{{\,nm}}$"+'\n'+
                      rf"$a_{{tot}} = {a_tot*1e6:.2e} \mathrm{{\,mm^2}}$")
        if idx ==0:
            ax_frac.plot(a_frac_list*100, eta_list*100, label=long_label)
        else:
            ax_frac.plot(a_frac_list*100, eta_list*100, label=short_label)

        short_label = (rf"{eta_pred*100:.1f}%@$f_{{A_{{Al,act}}}} = {a_ind_frac*100:.3f}\%$")
        ax_frac.plot(a_ind_frac*100, eta_pred*100, label=short_label, marker='o', markersize=8, linestyle='None')

    # Plot formatting: Volume
    ax_vol.set_xscale('log')
    ax_vol.set_yscale('log')
    ax_vol.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    ax_vol.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    ax_vol.set_xlabel(r'$V_{\mathrm{Al,act}}$ [$\mu$m$^3$]')
    ax_vol.set_ylabel(r'$\eta_{\mathrm{tot}}$ (%)')
    ax_vol.set_title('Total Efficiency vs Al Volume')
    ax_vol.grid(True)
    ax_vol.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    # Plot formatting: Area
    ax_area.set_xscale('log')
    ax_area.set_yscale('log')
    ax_area.set_xlabel(r'$A_{\mathrm{Al,act}}$ [$\mu$m$^2$]')
    ax_area.set_ylabel(r'$\eta_{\mathrm{ph}}$ (%)')
    ax_area.set_title('Phonon Collection Efficiency vs Device Area')
    ax_area.grid(True)
    ax_area.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    # Plot formatting: Area
    ax_frac.set_xscale('log')
    ax_frac.set_yscale('log')
    ax_frac.set_xlabel(r'$f_{A_{\mathrm{Al,act}}}$ (%)')
    ax_frac.set_ylabel(r'$\eta_{\mathrm{ph}}$ (%)')
    ax_frac.set_title('Phonon Collection Efficiency vs area coverage')
    ax_frac.grid(True)
    ax_frac.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    # fig.subplots_adjust(right=0.75)  # Add space for legend

    # Save figure
    save_dir = os.path.dirname(plot_dir)
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir+".pdf", dpi=300)
    plt.savefig(plot_dir+".png", dpi=300)
    plt.close()

    save_subplots(axs, plot_dir)

    return True

def df_kid_eff_area():
    # Build list of rows
    rows = []
    for name, vals in area_kid.items():
        v_al_act = vals["v_al_act"]
        v_al = vals["v_al"]
        t_al = vals["t_al"]
        a_tot = vals["a_tot"]
        v_nb = vals["v_nb"]
        eff_meas = vals.get("eff_measured", 0)
        eff_meas_ph = eff_meas/eta_pb
        
        f_al_act = v_al_act / (t_al * a_tot)  # area fraction
        tau_coll = tau_collect_eq(eta_sub_nom, f_al_act, n_bar_abs, c_s_ph)
        if eff_meas!=0: 
            tau_life = infer_tau_life(eff_meas_ph, tau_coll)
        else: 
            tau_life=0
        
        rows.append({
            "device": name,
            "v_al_act [m³]": v_al_act,
            "v_al [m³]": v_al,
            "v_nb [m³]": v_nb,
            "t_al [m]": t_al,
            "a_tot [m²]": a_tot,
            "f_al_act": f_al_act,
            "tau_life [s]": tau_life,
            "eff_measured_ph": eff_meas_ph,
        })

    # Create DataFrame
    df = pd.DataFrame(rows)
    # df.set_index("device", inplace=True)
    print(df)
    return df 

def plot_ph_eff_area(plot_dir, plot_log=False, large_kid=False):
    plt.figure(figsize=(8, 6))
    f_act = np.linspace(1e-5, 1, 1000)  # hz
    tau_collect_scdms = tau_collect_eq(eta_sub_scdms, f_act, n_bar_abs, c_s_ph)
    eff_tau_scdms = phonon_eff_lifetime(tau_life_ph, tau_collect_scdms)

    long_label_tau = (rf"$\eta = {eta_sub_scdms*1e3:.0f}\,\mathrm{{mm}}$"+'\n'+
                     rf"$\overline{{n}}_{{Si-Al}} = {n_bar_abs:.3f}$"+'\n'+
                     rf"$c_s = {c_s_ph}\,\mathrm{{m/s}}$"+'\n'+
                     rf"$\tau_{{ph}} = {tau_life_ph*1e3}\,\mathrm{{ms}}$")
    plt.plot(f_act*100, eff_tau_scdms*100, label=long_label_tau)  # Convert to μs

    hvev_f_act = area_ph["HVeV-NFC"]["a_frac"]
    hvev_eff = area_ph["HVeV-NFC"]["eff_measured"]
    long_label_tau = (rf"{hvev_eff*100}%@$f_{{Al,act}}={hvev_f_act*100}\%$")
    plt.plot(hvev_f_act*100, hvev_eff*100, label=long_label_tau, marker='o', markersize=8, linestyle='None')  # Convert to μs

    df_kid = df_kid_eff_area()

    # Get the second device by index
    kid = df_kid.index[0]
    row_kid = df_kid.loc[kid]

    tau_life = row_kid["tau_life [s]"]
    eff_measured_ph = row_kid["eff_measured_ph"]  # If not already in df
    f_al = row_kid["f_al_act"]
    
    # Model curve
    tau_collect_vals = tau_collect_eq(eta_sub_nom, f_act, n_bar_abs, c_s_ph)
    eff_model = phonon_eff_lifetime(tau_life, tau_collect_vals)

    # Plot model curve
    long_label_tau = (rf"{row_kid['device']}$@\eta = {eta_sub_nom*1e3:.0f}\,\mathrm{{mm}}$"+'\n'+
                 rf"$\overline{{n}}_{{Si-Al}} = {n_bar_abs:.3f}$"+'\n'+
                 rf"$c_s = {c_s_ph}\,\mathrm{{m/s}}$"+'\n'+
                 rf"$\tau_{{ph}} = {tau_life*1e6:.2f}\,\mathrm{{\mu s}}$"+'\n'+
                 rf"$\eta_{{pb}} = {eta_pb*100}\%$")
    plt.plot(f_act * 100, eff_model * 100, label=long_label_tau)
    short_label = (rf"{eff_measured_ph*100:.1f}%@$f_{{A_{{Al,act}}}} = {f_al*100:.3f}\%$")
    plt.plot(f_al * 100, eff_measured_ph * 100, marker='o', linestyle='None', markersize=8, label=short_label)

    # Get the second device by index
    paa = df_kid.index[1]
    row = df_kid.loc[paa]

    # Extract values
    f_al_paa = row["f_al_act"]
    tau_col_paa = tau_collect_eq(eta_sub_nom, f_al_paa, n_bar_abs, c_s_ph)
    eff_paa = phonon_eff_lifetime(tau_life, tau_col_paa)
    short_label = (
        rf"{row['device']}@{eff_paa*100:.1f}%"+'\n'+
        rf"$f_{{A_{{Al,act}}}} = {f_al_paa*100:.3f}\%$")
    plt.plot(f_al_paa * 100, eff_paa * 100, marker='o', linestyle='None', markersize=8, label=short_label)

    if large_kid: 
        tau_col_large_kid = tau_collect_eq(eta_sub_nom, f_al*4, n_bar_abs, c_s_ph)
        eff_large_kid = phonon_eff_lifetime(tau_life, tau_col_large_kid)
        short_label = (
            rf"{row_kid['device']}@{eff_large_kid*100:.1f}%"+'\n'+
            rf"$f_{{A_{{Al,act}}}} = {f_al*400:.3f}\%$")
        plt.plot(f_al* 400, eff_large_kid * 100, marker='*', linestyle='None', markersize=8, label=short_label)

    eff_vol_conv = phonon_collection_efficiency_area(f_act, row["a_tot [m²]"], row["t_al [m]"], 
        row["v_nb [m³]"], xi_nb_al_ratio, r_e_loss_xi_al)
    short_label = (
          rf"Volume@$\xi_{{Nb}}/\xi_{{Al}} = {xi_nb_al_ratio}$"
    #       rf'$E_{{mount}}/\xi_{{Al}} = {r_e_loss_xi_al*1e18:.3e}\mathrm{{\,\mu m^3}}$'+'\n'+
    #       # rf"$V_{{Al,tot}} = {row["v_al [m³]"]*1e18:.3e} \mathrm{{\,\mu m^3}}$"+'\n'+
    #       rf"$V_{{Nb}} = {row["v_nb [m³]"]*1e18} \mathrm{{\,\mu m^3}}$"+'\n'+
    #       rf"$t_{{Al}} = {row["t_al [m]"]*1e9:.2e} \mathrm{{\,nm}}$"+'\n'+
    #       rf"$a_{{tot}} = {row["a_tot [m²]"]*1e6:.2e} \mathrm{{\,mm^2}}$")
    # plt.plot(f_act * 100, eff_vol_conv * 100, '--', label=short_label
    )

    eff_vol_paa = phonon_collection_efficiency_area(f_al_paa, row["a_tot [m²]"], row["t_al [m]"], 
        row["v_nb [m³]"], xi_nb_al_ratio, r_e_loss_xi_al)
    short_label = (
        rf"{row['device']}@{eff_vol_paa*100:.1f}%"+'\n'+
        rf"$f_{{A_{{Al,act}}}} = {f_al_paa*100:.3f}\%$")
    plt.plot(f_al_paa * 100, eff_vol_paa * 100, marker='o', linestyle='None', markersize=8, label=short_label)

    plt.xlabel('Active Al Area Coverage (%)')
    plt.ylabel(r'Phonon Collection Efficiency (%)')
    plt.title('Phonon Collection Efficiency vs Absorber Area Coverage')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()

    if plot_log: 
        plt.yscale('log')
        plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
        plt.xscale('log')
        plt.gca().xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    save_dir = os.path.dirname(plot_dir)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

    return True

def plot_ph_eff_area_clean(plot_dir, plot_log=True):
    n_x = 2
    n_y = 2
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()  # flatten to 1D array for easier indexing

    x_labels = [
        'Active Al Area Coverage (%)',
        'Active Al Area Coverage (%)',
        'Active Al Area Coverage (%)',
        'Active Al Area Coverage (%)',
    ]

    y_labels = [
        r'Phonon Collection Efficiency (%)',
        r'Phonon Collection Efficiency (%)',
        r"$\tau_{col}\mathrm{\,[\mu s]}$",
        r"$\tau_{col}\mathrm{\,[\mu s]}$",
    ]

    titles = [
        r'SCDMS-eff',
        r'KIPM-eff',
        r'SCDMS-tau_col',
        r'KIPM-tau_col',
    ]

    f_act = np.linspace(1e-4, 1, 10000)  # hz
    tau_collect_scdms = tau_collect_eq(eta_sub_scdms, f_act, n_bar_abs, c_s_ph)
    eff_scdms = phonon_eff_lifetime(tau_life_ph_nom, tau_collect_scdms)

    # Model curve
    tau_collect_kid = tau_collect_eq(eta_sub_nom, f_act, n_bar_abs, c_s_ph)
    eff_kid = phonon_eff_lifetime(tau_life_ph_nom, tau_collect_kid)

    df_kid = df_kid_eff_area()

    f_act_qet = area_ph["HVeV-NFC"]["a_frac"]
    eff_qet = area_ph["HVeV-NFC"]["eff_measured"]
    qet_label = (
        rf"QET: {eff_qet*100:.1f}%"+ "\n"
        + rf"$f_{{\mathrm{{Al,act}}}} = {f_act_qet*100:.1f}\%$"
    )

    # Get the second device by index
    kid = df_kid.index[0]
    row_kid = df_kid.loc[kid]
    eff_obs_kid = row_kid["eff_measured_ph"]  # If not already in df
    f_act_kid = row_kid["f_al_act"]
    kid_label = (rf"KIPM: {eff_obs_kid*100:.1f}%"+ "\n"
        + rf"$f_{{\mathrm{{Al,act}}}} = {f_act_kid*100:.2f}\%$")

    # Get the second device by index
    paa = df_kid.index[1]
    row_paa = df_kid.loc[paa]
    f_act_paa = row_paa["f_al_act"]
    tau_col_paa = tau_collect_eq(eta_sub_nom, f_act_paa, n_bar_abs, c_s_ph)
    eff_paa = phonon_eff_lifetime(tau_life_ph_nom, tau_col_paa)
    paa_label = (rf"PAA-KIPM (exp): {eff_paa*100:.1f}%"+ "\n"
        + rf"$f_{{\mathrm{{Al,act}}}} = {f_act_paa*100:.1f}\%$")

    axs[2].plot(f_act*100, tau_collect_scdms*1e6, label='SCDMS')  # Convert to μs
    axs[2].plot(f_act*100, tau_collect_kid*1e6, label='KIPM')  # Convert to μs

    axs[0].plot(f_act*100, eff_scdms*100, label='SCDMS')  # Convert to μs
    axs[0].plot(f_act * 100, eff_kid*100, label='KIPM')

    axs[0].plot(f_act_qet*100, eff_qet*100, label=qet_label, 
        marker='o', markersize=8, linestyle='None')  # Convert to μs
    axs[0].plot(f_act_kid * 100, eff_obs_kid * 100, marker='o', 
        linestyle='None', markersize=8, label=kid_label)
    axs[0].plot(f_act_paa * 100, eff_paa * 100, marker='o', 
        linestyle='None', markersize=8, label=paa_label)

    for ind, ax in enumerate(axs):
        ax.set_xlabel(x_labels[ind])
        ax.set_ylabel(y_labels[ind])
        ax.set_title(titles[ind])
        ax.grid(True)
        if plot_log: 
            ax.set_xscale('log')
            ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
            ax.set_yscale('log')
            ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    for ax in [axs[0], axs[2]]:
        ax.legend()

    plt.tight_layout()

    save_dir = os.path.dirname(plot_dir)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

    save_each_axes(fig, axs, plot_dir)

    return True



