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
from matplotlib.patches import Circle

# integrate diffusion
from scipy import integrate
import numpy as np
from numpy import sign
from decimal import Decimal

# plot tau_s 
from scipy.constants import hbar, electron_volt
import math
from scipy.constants import Boltzmann, e

import sys
from pathlib import Path

# Path to the directory where this script lives
here = Path(__file__).resolve().parent

# Append ../res and ../mb relative to config/
sys.path.append(str(here.parent / "res"))
sys.path.append(str(here.parent / "eff"))

from plot_eff import save_subplots, save_each_axes
from config.const_config import *
from resEquations import *
from scaleEq import *

rcParams.update({'font.size': 20})

resolution_all = {
    "PAA":{
        "GR-diss": -gr_paa_diss,
        "GR-freq": gr_paa_freq,
        "AMP-freq": amp_ds21_res_paa_vol,
        "AMP-diss": -amp_ds21_res_paa_vol,
        "TLS-freq": tls_s21_freq_paa,
        "Total-freq": tot_freq_paa_s21,
        "Total-diss": -tot_diss_paa_s21,},
    "KID":{
        "GR": gr_kid,
        "AMP-freq": amp_eabs_res_kid_freq_vol,
        "AMP-diss": amp_eabs_res_kid_diss_vol,
        "TLS-freq": tls_eabs_kid,
        "Total-freq": tot_freq_kid,
        "Total-diss": tot_diss_kid,}
}

def find_intercept(x, y, y_target):
    """
    Find the x-value where y(x) crosses y_target.
    Uses simple linear interpolation between nearest points.
    """
    diff = y - y_target
    sign_change = sign(diff[:-1]) != sign(diff[1:])
    idx = np.where(sign_change)[0]
    if len(idx) == 0:
        return None  # no crossing found
    i = idx[0]
    # linear interpolation
    x0, x1 = x[i], x[i+1]
    y0, y1 = y[i], y[i+1]
    return x0 + (y_target - y0) * (x1 - x0) / (y1 - y0)

def plot_threshold(plot_dir, dev="PAA"):
    """
    Plot TLS-limited energy resolution vs different detector parameters in a 2x3 grid.
    Mirrors gr_res_vs_eabs.

    Parameters
    ----------
    plot_dir : str
        Directory where the plot will be saved.
    """
    n_x = 1
    n_y = 2
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()

    # Define variables
    e_abs_list = np.linspace(0, 10, 1000)  # meV

    x_labels = [
        r"$E$ (meV)",
        r"$E$ (meV)",
    ]
    y_labels = [
        r'$r(S_{21}(E_{abs}))-r_0$',
        r'$\theta(S_{21}(E_{abs}))-\theta_0$',
    ]

    r_eabs, xc_eabs = s21_circle_radius(e_abs_list*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
    s21_paa_fr = s21_ideal_eabs(f_0_nom, e_abs_list*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
    theta_eabs = get_unwrapped_phase_deg(s21_paa_fr-xc_eabs)

    # --- Panel 5/6: arc/phase direction ---
    axs[0].plot(e_abs_list, r_eabs-r_eabs[0], label=r"$E_{\mathrm{abs}}$")
    axs[0].plot(e_abs_list/total_eff_exp, r_eabs-r_eabs[0], label=r"$E_{\mathrm{dep}}$")
    axs[1].plot(e_abs_list, theta_eabs-theta_eabs[0], label=r"$E_{\mathrm{abs}}$")
    axs[1].plot(e_abs_list/total_eff_exp, theta_eabs-theta_eabs[0], label=r"$E_{\mathrm{dep}}$")
    x_dep = e_abs_list / total_eff_exp
    y_r = r_eabs-r_eabs[0]
    y_theta = theta_eabs - theta_eabs[0]

    # --- Overlay horizontal lines for PAA resolutions ---
    for key, val in resolution_all[dev].items():
        if "diss" in key.lower():
            x_cross = find_intercept(x_dep, y_r, val)
            axs[0].axhline(val, linestyle="--", label=f"{key}: {x_cross:.1f} meV")
        elif "freq" in key.lower():
            val_theta = angle_diff_from_im_shift(r_eabs[0], val)
            x_cross = find_intercept(x_dep, y_theta, val_theta)
            axs[1].axhline(val_theta, linestyle="--", label=f"{key}: {x_cross:.1f} meV")

    # --- Axis labels, grids, legends ---
    for i, label in enumerate(x_labels):
        axs[i].set_xlabel(label)
        axs[i].set_ylabel(y_labels[i])
        axs[i].grid(True)

    for ax in [axs[0], axs[1]]:
        ax.legend()

    plt.tight_layout()
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

def scale_nqp_plot(plot_dir):
    n_x = 2
    n_y = 2
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()

    x_labels = [
        r"$n_{qp,0}$ ($\mu m^{-3}$)",
        r"$n_{qp,0}$ ($\mu m^{-3}$)",
        r"$n_{qp,0}$ ($\mu m^{-3}$)",
        r"$n_{qp,0}$ ($\mu m^{-3}$)",
    ]
    y_labels = [
        r'$T_{eff}$ (mK)',
        r'$\kappa$ ($\mu m^3$)',
        r'$\sigma_{E_{abs}}$ (meV)',
        r'$T_{eff}$ (mK)',
    ]

    t_eff_list = t_eff_nqp_scan(delta_0_hf, N_0_hf)
    j_dff_tls_paa_1khz_list = update_tls_psd(t_eff_list, t_eff_music, v_c_paa, v_c_music, 
        eres_paa, eres_music, j_dff_tls_music_1khz, tls_beta)
    tls_dff_paa_list = tls_variance(tau_r_target, j_dff_tls_paa_1khz_list, froll_paa, deltaf_paa)
    tls_eabs_paa_list = convert_tls_res_to_eabs_res(tls_dff_paa_list, vol_paa, delta_0_hf, 
        alpha_paa, gamma_nom, k2_paa)
    kappa_1_list = kappa_1(t_eff_list, f_0_nom, delta_0_hf, N_0_hf)
    kappa_2_list = kappa_2(t_eff_list, f_0_nom, delta_0_hf, N_0_hf)
    gr_paa_list = compute_gr_resolution(nqp0_scan_list, vol_paa, delta_0_hf)
    tot_freq_paa = compute_total_resolution_list([gr_paa_list, 
        tls_eabs_paa_list, amp_eabs_res_paa_freq_vol])

    # --- Panel 5/6: arc/phase direction ---
    axs[0].plot(nqp0_scan_list*1e-18, t_eff_list*1e3)
    axs[1].plot(nqp0_scan_list*1e-18, kappa_1_list*1e18, label=r'$\kappa_1$')
    axs[1].plot(nqp0_scan_list*1e-18, kappa_2_list*1e18, label=r'$\kappa_2$')
    axs[2].plot(nqp0_scan_list*1e-18, tot_freq_paa*1e3, label=r'Total (freq)')
    axs[2].plot(nqp0_scan_list*1e-18, gr_paa_list*1e3, label=r'GR (freq)')
    axs[2].plot(nqp0_scan_list*1e-18, tls_eabs_paa_list*1e3, label=r'TLS (freq)')
    # axs[0].plot(e_abs_list/total_eff_exp, r_eabs-r_eabs[0])

    # --- Axis labels, grids, legends ---
    for i, label in enumerate(x_labels):
        axs[i].set_xlabel(label)
        axs[i].set_ylabel(y_labels[i])
        axs[i].grid(True)

    for ax in [axs[1], axs[2]]:
        ax.legend()

    # Save figure
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    fig.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    fig.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close(fig)

    save_each_axes(fig, axs, plot_dir)

def scale_taur_plot(plot_dir):
    n_x = 1
    n_y = 2
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()

    x_labels = [
        r"$\tau_{r}$ (ms)",
    ]
    y_labels = [
        r'$\sigma_{E_{abs}}$ (meV)',
    ]

    taur_list = np.linspace(2*1e-5, 2*1e-1, 1000) # s
    amp_ds21_res_paa_vol = compute_amp_resolution(taur_list, tn_nom, pfeed_paa_vol, debug=False)
    amp_eabs_res_paa_diss_vol, amp_eabs_res_paa_freq_vol = convert_amp_res_to_eabs_res(amp_ds21_res_paa_vol,
                vol_paa, delta_0_hf, alpha_paa, gamma_nom, k1_paa, k2_paa, qc0_nom, qr0_nom, debug=False)
    deltaf_paa_list = compute_delta_f(taur_list)

    tls_dff_paa_list = []
    for idx, taur in enumerate(taur_list):
        tls_dff_paa = tls_variance(taur, j_dff_tls_paa_1khz, froll_paa, deltaf_paa_list[idx])
        tls_dff_paa_list.append(tls_dff_paa)
    tls_dff_paa_list = np.asarray(tls_dff_paa_list)
    tls_eabs_paa = convert_tls_res_to_eabs_res(tls_dff_paa_list, vol_paa, delta_0_hf, 
        alpha_paa, gamma_nom, k2_paa)
    tot_freq_paa = compute_total_resolution_list([gr_paa, 
            tls_eabs_paa, amp_eabs_res_paa_freq_vol])

    # --- Panel 5/6: arc/phase direction ---
    axs[0].plot(taur_list*1e3, tot_freq_paa*1e3, label=r'Total (freq)')
    axs[0].plot(taur_list*1e3, tls_eabs_paa*1e3, label=r'TLS (freq)')
    axs[0].plot(taur_list*1e3, amp_eabs_res_paa_freq_vol*1e3, label=r'AMP (freq)')
    # axs[0].plot(e_abs_list/total_eff_exp, r_eabs-r_eabs[0])

    # --- Axis labels, grids, legends ---
    for i, label in enumerate(x_labels):
        axs[i].set_xlabel(label)
        axs[i].set_ylabel(y_labels[i])
        axs[i].grid(True)
        axs[i].set_xscale('log')
        # axs[i].set_yscale('log')
        axs[i].xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
        # axs[i].yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    for ax in [axs[0]]:
        ax.legend()

    # Save figure
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    fig.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    fig.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close(fig)

    save_each_axes(fig, axs, plot_dir)

def scale_qi_plot(plot_dir):
    n_x = 2
    n_y = 3
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()

    x_labels = [
        r'$\sigma_{E_{abs}}$ (meV)',
        r'$\sigma_{E_{abs}}$ (meV)',
        r"$Q_i$",
        r"$Q_i$",
    ]
    y_labels = [
        r"$Q_i$",
        r'$|S_{21}(f_{r,0})|$',
        r'$P_{feed}$ (dBm)',
        r'$\sigma_{E_{abs}}$ (meV)',
    ]

    # qi_list = np.linspace(1e5, 1e7, 10) # s
    qi0_list = [1e5, 200*1e3, 500*1e3, 800*1e3, 1e6, 5*1e6, 1e7]
    qi0_list_smooth = np.linspace(1e5, 1e7, 1000) # s
    e_abs = np.linspace(0.01, 500, 1000)  # meV
    for qi0 in qi0_list: 
        s21_qi = s21_ideal_eabs(f_0_nom, e_abs*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0, qc0_nom)
        qi_eabs = eabs_to_qi(e_abs*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0)
        s21_qi_db = s21_z_to_mag(s21_qi)
        axs[0].plot(e_abs, qi_eabs, label=rf'$Q_{{i,0}}$ = {qi0}')
        axs[1].plot(e_abs, s21_qi_db, label=rf'$Q_{{i,0}}$ = {qi0}')

    qr0_list = calculate_Qr(qi0_list_smooth, qc0_nom)
    pfeed_paa_vol = P_bif(N_0_hf, delta_0_hf, vol_paa, f_0_nom, qc0_nom, alpha_paa, qr0_list, debug=False)
    pfeed_paa_vol_dBm = power_to_dbm(pfeed_paa_vol, debug=debug)

    tls_dff_paa_list = []
    amp_eabs_res_paa_freq_vol_list = []
    amp_ds21_res_paa_vol = compute_amp_resolution(tau_r_target, tn_nom, pfeed_paa_vol, debug=False)

    froll_paa = compute_f_rolloff(f_0_nom, qr0_list)
    froll_music = compute_f_rolloff(f_0_music, qr0_list)

    j_dff_tls_paa_1khz = update_tls_psd(t_eff_hf, t_eff_music, v_c_paa, v_c_music, 
        eres_paa, eres_music, j_dff_tls_music_1khz, tls_beta)
    for idx, qr0 in enumerate(qr0_list):
        amp_eabs_res_paa_diss_vol, amp_eabs_res_paa_freq_vol = convert_amp_res_to_eabs_res(amp_ds21_res_paa_vol[idx],
                vol_paa, delta_0_hf, alpha_paa, gamma_nom, k1_paa, k2_paa, qc0_nom, qr0, debug=False)
        amp_eabs_res_paa_freq_vol_list.append(amp_eabs_res_paa_freq_vol)

        wres_paa_vol = W_er(qr0, qc0_nom, f_0_nom, pfeed_paa_vol[idx])
        eres_paa_local = compute_E_field(wres_paa_vol, c_paa, t_a_si_paa)
        j_dff_tls_paa_1khz = update_tls_psd(t_eff_hf, t_eff_music, v_c_paa, v_c_music, 
            eres_paa_local, eres_music, j_dff_tls_music_1khz, tls_beta)
        tls_dff_paa = tls_variance(tau_r_target, j_dff_tls_paa_1khz, froll_paa[idx], deltaf_paa)
        tls_dff_paa_list.append(tls_dff_paa)
    tls_dff_paa_list = np.asarray(tls_dff_paa_list)
    amp_eabs_res_paa_freq_vol_list = np.asarray(amp_eabs_res_paa_freq_vol_list)
    tls_eabs_paa = convert_tls_res_to_eabs_res(tls_dff_paa_list, vol_paa, delta_0_hf, 
        alpha_paa, gamma_nom, k2_paa)
    tot_freq_paa = compute_total_resolution_list([gr_paa, 
            tls_eabs_paa, amp_eabs_res_paa_freq_vol_list])

    # # --- Panel 5/6: arc/phase direction ---
    axs[2].plot(qi0_list_smooth, pfeed_paa_vol_dBm)
    axs[3].plot(qi0_list_smooth, tot_freq_paa*1e3, label=r'Total (freq)')
    axs[3].plot(qi0_list_smooth, tls_eabs_paa*1e3, label=r'TLS (freq)')
    axs[3].plot(qi0_list_smooth, amp_eabs_res_paa_freq_vol_list*1e3, label=r'AMP (freq)')
    # axs[0].plot(e_abs_list/total_eff_exp, r_eabs-r_eabs[0])

    # --- Axis labels, grids, legends ---
    for i, label in enumerate(x_labels):
        axs[i].set_xlabel(label)
        axs[i].set_ylabel(y_labels[i])
        axs[i].grid(True)

    for ax in [axs[1], axs[3]]:
        ax.legend()
        ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
        # ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    for ax in [axs[2]]:
        ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
        # ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    for ax in [axs[0]]:
        ax.legend()
        # ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # reserve 15% of width for legends

    # Save figure
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    fig.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    fig.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close(fig)

    save_each_axes(fig, axs, plot_dir)

def scale_nkids_plot(plot_dir):
    n_x = 1
    n_y = 2
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()

    x_labels = [
        r'$N_{KIDs}$',
    ]
    y_labels = [
        r"$\sigma_{E_{abs}}$ (meV)",
    ]

    # qi_list = np.linspace(1e5, 1e7, 10) # s
    nkid_list = list(range(1, 51))
    cover_list = area_coverage(l_abs**2*l_abs_num*np.asarray(nkid_list), l_sub**2*2)

    # # --- Panel 5/6: arc/phase direction ---
    axs[0].plot(nkid_list, tot_freq_paa*np.sqrt(nkid_list)*1e3, color="tab:blue", label="Resolution")
    ax0_twin = axs[0].twinx()
    line2, = ax0_twin.plot(nkid_list, cover_list, color="tab:red", label="Area coverage")
    ax0_twin.set_ylabel(r"Area coverage ($\%$)", color="tab:red")

    # --- Axis labels, grids, legends ---
    for i, label in enumerate(x_labels):
        axs[i].set_xlabel(label)
        axs[i].set_ylabel(y_labels[i])
        axs[i].grid(True)

    # Collect handles/labels from both axes
    handles1, labels1 = axs[0].get_legend_handles_labels()
    handles2, labels2 = ax0_twin.get_legend_handles_labels()

    # Put a single combined legend outside the plot (to the right)
    axs[0].legend(handles1 + handles2, labels1 + labels2, frameon=True)

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # reserve 15% of width for legends

    # Save figure
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    fig.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    fig.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close(fig)

    save_each_axes(fig, axs, plot_dir)

def scale_vol_plot(plot_dir):
    n_x = 3
    n_y = 3
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()

    x_labels = [
        r'Inductor width ($\mu m$)',
        r'Inductor width ($\mu m$)',
        r'Inductor width ($\mu m$)',
        r'Inductor width ($\mu m$)',
    ]
    y_labels = [
        r'Inductor length ($\mu m$)',
        r'Inductor length ($\mu m$)',
        r'Inductor length ($\mu m$)',
        r'Inductor length ($\mu m$)',
    ]

    w_ind_list = np.linspace(1e-6, 1e-5, 100) # mum
    l_ind_list = np.linspace(l_ind_paa, l_ind_paa*10, 100) # s
    # Create 2D grid
    W, L = np.meshgrid(w_ind_list, l_ind_list)
    vol = W * L * t_ind_paa   # element-wise product, in m^2
    kin_indhenry_tot = L/W * indhenry_sq
    alpha_list = kin_indhenry_tot/(indhenry_geo_paa+kin_indhenry_tot)
    fr_list = calculate_resonant_frequency(c_paa, kin_indhenry_tot+indhenry_geo_paa)

    gr_paa_list = compute_gr_resolution(nqp_target, vol, delta_0_hf)
    pfeed_paa_list = P_bif(N_0_hf, delta_0_hf, vol, f_0_nom, qc0_nom, alpha_list, qr0_nom, debug=False)
    pfeed_paa_list_dbm = power_to_dbm(pfeed_paa_list, debug=debug)

    amp_ds21_res_paa_list = compute_amp_resolution(tau_r_target, tn_nom, pfeed_paa_list, debug=debug)
    amp_eabs_res_paa_diss_list, amp_eabs_res_paa_freq_list = convert_amp_res_to_eabs_res(amp_ds21_res_paa_list,
                vol, delta_0_hf, alpha_list, gamma_nom, k1_paa, k2_paa, qc0_nom, qr0_nom, debug=debug)

    # # --- Panel 5/6: arc/phase direction ---
    c0 = axs[0].pcolormesh(W*1e6, L*1e6, vol*1e18, shading="auto", cmap="viridis")

    c1 = axs[1].pcolormesh(W*1e6, L*1e6, kin_indhenry_tot*1e9, shading="auto", cmap="viridis")
    # Overlay contour at 10 nH
    cont = axs[1].contour(W*1e6, L*1e6, kin_indhenry_tot*1e9, 
        levels=[3, 5, 10], colors=["brown", "orange", "red"], linewidths=2)
    # Add labels onto the contour lines
    axs[1].clabel(cont, fmt={3: "3 nH", 5: "5 nH", 10: "10 nH"}, 
        inline=True, colors=["brown", "orange", "red"])

    c2 = axs[2].pcolormesh(W*1e6, L*1e6, alpha_list, shading="auto", cmap="viridis")
    cont = axs[2].contour(W*1e6, L*1e6, alpha_list, 
        levels=[0.2, 0.4, 0.5], colors=["brown", "orange", "red"], linewidths=2)
    # Add labels onto the contour lines
    axs[2].clabel(cont, fmt={0.2: "0.2", 0.4: "0.4", 0.5: "0.5"}, 
        inline=True, colors=["brown", "orange", "red"])

    c3 = axs[3].pcolormesh(W*1e6, L*1e6, fr_list*1e-9, shading="auto", cmap="viridis")
    cont = axs[3].contour(W*1e6, L*1e6, fr_list*1e-9, 
        levels=[4, 4.5, 5], colors=["brown", "orange", "red"], linewidths=2)
    # Add labels onto the contour lines
    axs[3].clabel(cont, fmt={4: "4 GHz", 4.5: "4.5 GHz", 5: "5 GHz"}, 
        inline=True, colors=["brown", "orange", "red"])

    c4 = axs[4].pcolormesh(W*1e6, L*1e6, gr_paa_list*1e3, shading="auto", cmap="viridis")

    c5 = axs[5].pcolormesh(W*1e6, L*1e6, pfeed_paa_list_dbm, shading="auto", cmap="viridis")
    cont = axs[5].contour(W*1e6, L*1e6, pfeed_paa_list_dbm, 
        levels=[-110, -100, -90], colors=["brown", "orange", "red"], linewidths=2)
    # Add labels onto the contour lines
    axs[5].clabel(cont, fmt={-110: "-110 dBm", -100: "-100 dBm", -90: "-90 dBm"}, 
        inline=True, colors=["brown", "orange", "red"])

    c6 = axs[6].pcolormesh(W*1e6, L*1e6, amp_eabs_res_paa_freq_list*1e3, shading="auto", cmap="viridis")
    # cont = axs[5].contour(W*1e6, L*1e6, pfeed_paa_list_dbm, 
    #     levels=[-110, -100, -90], colors=["brown", "orange", "red"], linewidths=2)
    # # Add labels onto the contour lines
    # axs[5].clabel(cont, fmt={-110: "-110 dBm", -100: "-100 dBm", -90: "-90 dBm"}, 
    #     inline=True, colors=["brown", "orange", "red"])

    cb = fig.colorbar(c0, ax=axs[0], label=r"Inductor volume [$\mu$m$^3$]")
    cb = fig.colorbar(c1, ax=axs[1], label=r"Total kinetic inductance [nH]")
    cb = fig.colorbar(c2, ax=axs[2], label=r"$\alpha$")
    cb = fig.colorbar(c3, ax=axs[3], label=r"$f_{r,0}$ (GHz)")
    cb = fig.colorbar(c4, ax=axs[4], label=r"GR $\sigma_{E_{abs}}$ (meV)")
    cb = fig.colorbar(c5, ax=axs[5], label=r"$P_{feed}$ (dBm)")
    cb = fig.colorbar(c6, ax=axs[6], label=r"AMP $\sigma_{E_{abs}}$ (meV)")

    # --- Axis labels, grids, legends ---
    for i, label in enumerate(x_labels):
        axs[i].set_xlabel(label)
        axs[i].set_ylabel(y_labels[i])
        axs[i].grid(True)

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # reserve 15% of width for legends

    # Save figure
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    fig.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    fig.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close(fig)

    save_each_axes(fig, axs, plot_dir)
