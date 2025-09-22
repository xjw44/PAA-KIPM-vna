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

rcParams.update({'font.size': 30})

# xi_nb_al_ratio=0.1292
# r_e_loss_xi_al=118000*1e-18
xi_nb_al_ratio=0.09531
r_e_loss_xi_al=413000*1e-18
eta_pb = 0.46

n_bar_abs = 0.91 # si to al 
c_s_ph = 5880 # m/s 
eta_sub = 1*1e-3 # mm
eta_sub_scdms = 4*1e-3 # mm
tau_life_ph = 2*1e-3 # ms

def save_subplots(axes, plot_dir, equalr=False):
    # Save individual subplots
    for i, ax_orig in enumerate(axes):
        # Create a new figure and axis
        fig_single, ax_single = plt.subplots(figsize=(8, 6))

        # Copy lines from the original axis
        for line in ax_orig.get_lines():
            ax_single.plot(*line.get_data(),
                           label=line.get_label(),
                           color=line.get_color(),
                           linestyle=line.get_linestyle(),
                           marker=line.get_marker(),
                           markersize=line.get_markersize())

        # --- Copy fill_between patches ---
        for coll in ax_orig.collections:
            if isinstance(coll, PolyCollection):
                # Extract vertices of the polygon(s)
                for path in coll.get_paths():
                    verts = path.vertices
                    x = verts[:, 0]
                    y = verts[:, 1]
                    # Because fill_between creates top and bottom polygons,
                    # you can just re-plot them
                    ax_single.fill_between(
                        x, y.min(), y.max(),
                        color=coll.get_facecolor()[0],
                        alpha=coll.get_alpha()
                    )

        # Copy axis formatting
        ax_single.set_xscale(ax_orig.get_xscale())
        ax_single.set_yscale(ax_orig.get_yscale())
        ax_single.set_xlabel(ax_orig.get_xlabel())
        ax_single.set_ylabel(ax_orig.get_ylabel())
        ax_single.set_title(ax_orig.get_title())
        if equalr: 
            ax_single.set_aspect('equal', adjustable='box')
        ax_single.grid(True)

        # Copy legend
        handles, labels = ax_orig.get_legend_handles_labels()
        if handles:
            # ax_single.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            # ax_single.legend(handles, labels, loc='lower left', frameon=True)
            # indp.legend(loc='lower left', frameon=True)
            # ax_single.legend(handles, labels)
            ax_single.legend(loc='best', frameon=True)

        # Save the individual figure
        # fig_single.tight_layout()
        fig_single.savefig(f"{plot_dir}_{i}.pdf", dpi=300, bbox_inches='tight')
        fig_single.savefig(f"{plot_dir}_{i}.png", dpi=300, bbox_inches='tight')
        plt.close(fig_single)

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
        tau_coll = tau_collect_eq(eta_sub, f_al_act, n_bar_abs, c_s_ph)
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
    tau_collect_vals = tau_collect_eq(eta_sub, f_act, n_bar_abs, c_s_ph)
    eff_model = phonon_eff_lifetime(tau_life, tau_collect_vals)

    # Plot model curve
    long_label_tau = (rf"{row_kid['device']}$@\eta = {eta_sub*1e3:.0f}\,\mathrm{{mm}}$"+'\n'+
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
    tau_col_paa = tau_collect_eq(eta_sub, f_al_paa, n_bar_abs, c_s_ph)
    eff_paa = phonon_eff_lifetime(tau_life, tau_col_paa)
    short_label = (
        rf"{row['device']}@{eff_paa*100:.1f}%"+'\n'+
        rf"$f_{{A_{{Al,act}}}} = {f_al_paa*100:.3f}\%$")
    plt.plot(f_al_paa * 100, eff_paa * 100, marker='o', linestyle='None', markersize=8, label=short_label)

    if large_kid: 
        tau_col_large_kid = tau_collect_eq(eta_sub, f_al*4, n_bar_abs, c_s_ph)
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


def save_each_axes(fig, axs, out_dir, prefix="subplot", dpi=300, pad_frac=0.02,
                   include_outside_legend=True):
    """
    Save each Axes from an existing figure by cropping the figure canvas
    to that Axes' (and legend's) tight bounding box.

    Works with lines, fill_between (PolyCollection), images, etc.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    axs : array-like of Axes
    out_dir : str
    prefix : str
    dpi : int
    pad_frac : float
        Extra fractional padding around the axes bbox.
    include_outside_legend : bool
        If True, expand bbox to include legend even if it's outside the axes.
    """
    os.makedirs(out_dir, exist_ok=True)
    axs = np.ravel(axs)

    # Ensure tightbbox values are up-to-date
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    for i, ax in enumerate(axs):
        bbox = ax.get_tightbbox(renderer)

        # Optionally include legend box if placed outside the axes
        leg = ax.get_legend()
        if include_outside_legend and leg is not None:
            leg_bbox = leg.get_window_extent(renderer)
            bbox = Bbox.union([bbox, leg_bbox])

        # Convert from pixels to inches and add a little padding
        bbox_in = bbox.transformed(fig.dpi_scale_trans.inverted())
        bbox_in = bbox_in.expanded(1 + pad_frac, 1 + pad_frac)

        out = os.path.join(out_dir, f"{prefix}_{i}.png")
        out_pdf = os.path.join(out_dir, f"{prefix}_{i}.pdf")
        fig.savefig(out, dpi=dpi, bbox_inches=bbox_in)
        fig.savefig(out_pdf, dpi=dpi, bbox_inches=bbox_in)
        print(f"Saved {out}")

