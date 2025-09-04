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
from decimal import Decimal

# plot tau_s 
from scipy.constants import hbar, electron_volt
import math
from scipy.constants import Boltzmann, e

import sys
# sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/res')
from pathlib import Path

# Path to the directory where this script lives
here = Path(__file__).resolve().parent

# Append ../res and ../mb relative to config/
sys.path.append(str(here.parent / "res"))

from config.const_config import *
from resEquations import *

rcParams.update({'font.size': 20})

y_label_psd = {
    "J_eabs": r"$J_{E_{abs}}\,[\mathrm{eV^2/Hz}]$",
    "J_dn_qp": r"$J_{\delta N_{qp}}\,[\mathrm{1/Hz}]$",
    "J_df/f": r"$J_{\delta f/f_{r,0}}\,[\mathrm{1/Hz}]$",
    "J_d1/Qi": r"$J_{\delta 1/Q_i}\,[\mathrm{1/Hz}]$",
    "J_Re(S21)": r"$J_{Re[\delta s21]}\,[\mathrm{1/Hz}]$",
    "J_Im(S21)": r"$J_{Im[\delta s21]}\,[\mathrm{1/Hz}]$"
}

legend_psd = {
    "GR": {
    'paa': "\n".join(["paa@" + label_nqp_target, label_tau_r_target, 
        label_delta_0_hf, label_vol_paa, label_alpha_paa, label_gamma_nom, 
        label_qi0_nom, label_qc0_nom, label_qr0_nom, label_f_0_nom, label_t_eff_hf, 
        label_N_0_hf, label_k1_paa, label_k2_paa]),
    'kid': "\n".join(["kid@" + label_nqp_target, label_tau_r_target, 
        label_delta_0_al, label_vol_kid, label_alpha_kid, label_gamma_nom, 
        label_qi0_nom, label_qc0_nom, label_qr0_nom, label_f_0_nom, label_t_eff_al, 
        label_N_0_al, label_k1_kid, label_k2_kid]),
    },
    "AMP": {
    'paa': "\n".join(["paa_freq@" + label_l_paa, label_t_ind_paa, label_w_ind_paa, 
        label_tn_nom, label_rho_nom_al, label_pfeed_paa, 
        label_delta_0_hf, label_vol_paa, label_alpha_paa, label_gamma_nom, 
        label_qi0_nom, label_qc0_nom, label_qr0_nom, label_f_0_nom, label_t_eff_hf, 
        label_N_0_hf, label_k1_paa, label_k2_paa]),
    'kid': "\n".join(["kid_freq@" + label_l_kid, label_t_ind_kid, label_w_ind_kid, 
        label_tn_nom, label_rho_nom_al, label_pfeed_kid, 
        label_delta_0_al, label_vol_kid, label_alpha_kid, label_gamma_nom, 
        label_qi0_nom, label_qc0_nom, label_qr0_nom, label_f_0_nom, label_t_eff_al, 
        label_N_0_al, label_k1_kid, label_k2_kid]),
    },
    "TLS": {
    'paa': "\n".join(["paa_freq@" + label_tls_beta, label_tls_n, 
        label_l_paa, label_t_ind_paa, label_w_ind_paa, 
        label_tn_nom, label_rho_nom_al, label_pfeed_paa, 
        label_delta_0_hf, label_vol_paa, label_alpha_paa, label_gamma_nom, 
        label_qi0_nom, label_qc0_nom, label_qr0_nom, label_f_0_nom, label_t_eff_hf, 
        label_N_0_hf, label_k1_paa, label_k2_paa]),
    'kid': "\n".join(["kid@" + label_nqp_target, label_tau_r_target, label_delta_0_al, label_vol_kid]),
    'music': "\n".join(["music_freq@" + label_l_paa, label_t_ind_paa, label_w_ind_paa, 
        label_tn_nom, label_rho_nom_al, label_pfeed_music, 
        label_delta_0_al, label_gamma_nom]),
    },
}

f_range = np.linspace(1*1e-3, 1000, 10000)  # hz
f_range_tls = np.linspace(1*1e-3, 1e5, 10000)  # hz

def convert_dict_df(psd_dict):
    # Build a DataFrame: each cell contains a full PSD array across frequencies
    df_psd_compact = pd.DataFrame({
        source: {
            obs: psd_dict[source][obs]
            for obs in psd_dict[source]
        }
        for source in psd_dict
    }).T  # Transpose so sources are rows

    df_psd_compact.index.name = "Noise Source"
    return df_psd_compact

def df_psd_all():
    j_eabs_gr_paa = J_GR_eabs(f_range, nqp_target, tau_r_target, vol_paa, delta_0_hf, debug=False)
    j_eabs_gr_kid = J_GR_eabs(f_range, nqp_target, tau_r_target, vol_kid, delta_0_al, debug=False)
    j_dNqp_gr_paa = convert_psd_eabs_to_dNqp(j_eabs_gr_paa, vol_paa, delta_0_hf)
    j_dNqp_gr_kid = convert_psd_eabs_to_dNqp(j_eabs_gr_kid, vol_kid, delta_0_al)
    j_dff_gr_paa = psd_dNqp_to_df_over_f(j_dNqp_gr_paa, alpha_paa, gamma_nom, k2_paa, vol_paa)
    j_dff_gr_kid = psd_dNqp_to_df_over_f(j_dNqp_gr_kid, alpha_kid, gamma_nom, k2_kid, vol_kid)
    j_d1qi_gr_paa = psd_dNqp_to_d1_over_Qi(j_dNqp_gr_paa, alpha_paa, gamma_nom, k1_paa, vol_paa)
    j_d1qi_gr_kid = psd_dNqp_to_d1_over_Qi(j_dNqp_gr_kid, alpha_kid, gamma_nom, k1_kid, vol_kid)
    j_imds21_gr_paa = psd_df_over_f_to_ImS21(j_dff_gr_paa, qr0_nom, qc0_nom)
    j_imds21_gr_kid = psd_df_over_f_to_ImS21(j_dff_gr_kid, qr0_nom, qc0_nom)
    j_reds21_gr_paa = psd_d1qi_to_ReS21(j_d1qi_gr_paa, qr0_nom, qc0_nom)
    j_reds21_gr_kid = psd_d1qi_to_ReS21(j_d1qi_gr_kid, qr0_nom, qc0_nom)

    j_amp_ds21_paa = amp_psd(tn_nom, pfeed_paa)
    j_amp_ds21_kid = amp_psd(tn_nom, pfeed_kid)
    j_amp_rest_paa = amp_psd_all(j_amp_ds21_paa, qr0_nom, qc0_nom, vol_kid, alpha_paa, gamma_nom, k1_paa, k2_paa, delta_0_hf)
    j_amp_rest_kid = amp_psd_all(j_amp_ds21_kid, qr0_nom, qc0_nom, vol_kid, alpha_kid, gamma_nom, k1_kid, k2_kid, delta_0_al)

    j_dff_tls_music = full_psd_tls(f_range_tls, j_dff_tls_music_1khz, froll_music, tls_n)
    j_dff_tls_paa = full_psd_tls(f_range_tls, j_dff_tls_paa_1khz, froll_paa, tls_n)
    j_tls_rest_music = convert_dff_tls_psd_to_all(j_dff_tls_music, vol_music, alpha_music, gamma_nom, 
        k2_music, delta_0_al, qr0_music, qc0_music)
    j_tls_rest_paa = convert_dff_tls_psd_to_all(j_dff_tls_paa, vol_music, alpha_paa, gamma_nom, 
        k2_paa, delta_0_hf, qr0_nom, qc0_nom)

    psd_dict_paa = {
        "GR": {
            "J_eabs": j_eabs_gr_paa,
            "J_dn_qp": j_dNqp_gr_paa,
            "J_df/f": j_dff_gr_paa,
            "J_d1/Qi": j_d1qi_gr_paa,
            "J_Re(S21)": j_reds21_gr_paa,
            "J_Im(S21)": j_imds21_gr_paa,
        },
        "AMP": {
            "J_eabs_freq": j_amp_rest_paa["J_eabs_freq"]*np.ones_like(f_range),
            "J_eabs_diss": j_amp_rest_paa["J_eabs_diss"]*np.ones_like(f_range),
            "J_dn_qp_freq": j_amp_rest_paa["J_dN_qp_freq"]*np.ones_like(f_range),
            "J_dn_qp_diss": j_amp_rest_paa["J_dN_qp_diss"]*np.ones_like(f_range),
            "J_df/f": j_amp_rest_paa["J_df/f"]*np.ones_like(f_range),     # replaced
            "J_d1/Qi": j_amp_rest_paa["J_d1/Qi"]*np.ones_like(f_range),    # replaced
            "J_Re(S21)": j_amp_ds21_paa*np.ones_like(f_range),  # replaced
            "J_Im(S21)": j_amp_ds21_paa*np.ones_like(f_range),  # replaced
        },
        "TLS": {
            "J_eabs": j_tls_rest_paa["J_eabs"],     # replaced
            "J_dn_qp": j_tls_rest_paa["J_dN_qp"],    # replaced
            "J_df/f": j_dff_tls_paa, 
            "J_d1/Qi": np.ones_like(f_range),    # replaced
            "J_Re(S21)": np.ones_like(f_range),  # replaced
            "J_Im(S21)": j_tls_rest_paa["J_Im(S21)"],  # replaced
        }
    }

    psd_dict_kid = {
        "GR": {
            "J_eabs": j_eabs_gr_kid,
            "J_dn_qp": j_dNqp_gr_kid,
            "J_df/f": j_dff_gr_kid,
            "J_d1/Qi": j_d1qi_gr_kid,
            "J_Re(S21)": j_reds21_gr_kid,
            "J_Im(S21)": j_imds21_gr_kid,
        },
        "AMP": {
            "J_eabs_freq": j_amp_rest_kid["J_eabs_freq"]*np.ones_like(f_range),
            "J_eabs_diss": j_amp_rest_kid["J_eabs_diss"]*np.ones_like(f_range),
            "J_dn_qp_freq": j_amp_rest_kid["J_dN_qp_freq"]*np.ones_like(f_range),
            "J_dn_qp_diss": j_amp_rest_kid["J_dN_qp_diss"]*np.ones_like(f_range),
            "J_df/f": j_amp_rest_kid["J_df/f"]*np.ones_like(f_range),     # replaced
            "J_d1/Qi": j_amp_rest_kid["J_d1/Qi"]*np.ones_like(f_range),    # replaced
            "J_Re(S21)": j_amp_ds21_kid*np.ones_like(f_range),  # replaced
            "J_Im(S21)": j_amp_ds21_kid*np.ones_like(f_range),  # replaced
        },
        "TLS": {
            "J_eabs": np.ones_like(f_range),     # replaced
            "J_dn_qp": np.ones_like(f_range),    # replaced
            "J_df/f": np.ones_like(f_range),     # replaced
            "J_d1/Qi": np.ones_like(f_range),    # replaced
            "J_Re(S21)": np.ones_like(f_range),  # replaced
            "J_Im(S21)": np.ones_like(f_range),  # replaced
        }
    }

    psd_dict_music = {
        "GR": {
            "J_eabs": j_eabs_gr_kid,
            "J_dn_qp": j_dNqp_gr_kid,
            "J_df/f": j_dff_gr_kid,
            "J_d1/Qi": j_d1qi_gr_kid,
            "J_Re(S21)": j_reds21_gr_kid,
            "J_Im(S21)": j_imds21_gr_kid,
        },
        "AMP": {
            "J_eabs_freq": j_amp_rest_kid["J_eabs_freq"]*np.ones_like(f_range),
            "J_eabs_diss": j_amp_rest_kid["J_eabs_diss"]*np.ones_like(f_range),
            "J_dn_qp_freq": j_amp_rest_kid["J_dN_qp_freq"]*np.ones_like(f_range),
            "J_dn_qp_diss": j_amp_rest_kid["J_dN_qp_diss"]*np.ones_like(f_range),
            "J_df/f": j_amp_rest_kid["J_df/f"]*np.ones_like(f_range),     # replaced
            "J_d1/Qi": j_amp_rest_kid["J_d1/Qi"]*np.ones_like(f_range),    # replaced
            "J_Re(S21)": j_amp_ds21_kid*np.ones_like(f_range),  # replaced
            "J_Im(S21)": j_amp_ds21_kid*np.ones_like(f_range),  # replaced
        },
        "TLS": {
            "J_eabs": j_tls_rest_music["J_eabs"],     # replaced
            "J_dn_qp": j_tls_rest_music["J_dN_qp"],    # replaced
            "J_df/f": j_dff_tls_music, 
            "J_d1/Qi": np.ones_like(f_range),    # replaced
            "J_Re(S21)": np.ones_like(f_range),  # replaced
            "J_Im(S21)": j_tls_rest_music["J_Im(S21)"],  # replaced
        }
    }

    df_psd_kid = convert_dict_df(psd_dict_kid)
    df_psd_paa = convert_dict_df(psd_dict_paa)
    df_psd_music = convert_dict_df(psd_dict_music)

    return df_psd_kid, df_psd_paa, df_psd_music

def res_all(debug=False):
    resolution_all = {
        "PAA":{
        "GR": gr_paa,
        "AMP-freq": amp_eabs_res_paa_freq,
        "AMP-diss": amp_eabs_res_paa_diss,
        "TLS-freq": tls_eabs_paa,
        "Total-freq": tot_freq_paa,
        "Total-diss": tot_diss_paa,},
        "KID":{
        "GR": gr_kid,
        "AMP-freq": amp_eabs_res_kid_freq,
        "AMP-diss": amp_eabs_res_kid_diss,
        "TLS-freq": tls_eabs_kid,
        "Total-freq": tot_freq_kid,
        "Total-diss": tot_diss_kid,}
    }
    rows = []
    for device, resolution in resolution_all.items():
        for noise, res in resolution.items():
            key = f"{device}-{noise}"
            rows.append((key, res))

    df_res = pd.DataFrame(rows, columns=["Label", "Resolution"]).set_index("Label")
    if debug:
        print(df_res)

    return df_res

def plot_psd_all(plot_dir):
    # Compute recombination constants
    df_psd_kid, df_psd_paa, df_psd_music = df_psd_all()
    observables = list(df_psd_kid.columns)
    print(observables)

    # Loop over each noise source
    for noise_source in df_psd_kid.index:
        n_x = 2
        n_y = 3
        fig, axes = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
        axes = axes.flatten()

        for i, obs in enumerate(observables):
            if "_freq" in obs or "_diss" in obs:
                continue  # Skip this observable
            ax = axes[i]
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
            ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
            ax.set_xlabel("Frequency [Hz]")
            ax.grid(True, which='both', linestyle='--', alpha=0.5)
            ax.set_title(f"{obs}_{noise_source}")
            ax.set_ylabel(y_label_psd.get(obs, "PSD [arb. units]"))
            # Plot PAA and KID PSDs
            if (noise_source == "AMP") and (obs=="J_dn_qp" or obs=="J_eabs"):
                obs_freq = obs + "_freq"
                obs_diss = obs + "_diss"
                ax.plot(f_range, df_psd_paa.loc[noise_source, obs_freq], label=legend_psd[noise_source]['paa'])
                ax.plot(f_range, df_psd_kid.loc[noise_source, obs_freq], linestyle='--', label=legend_psd[noise_source]['kid'])
                ax.plot(f_range, df_psd_paa.loc[noise_source, obs_diss], label=f"paa_diss")
                ax.plot(f_range, df_psd_kid.loc[noise_source, obs_diss], linestyle='--', label=f"kid_diss")
            elif (noise_source=="TLS"):
                ax.plot(f_range_tls, df_psd_paa.loc[noise_source, obs], label=legend_psd[noise_source]['paa'])
                ax.plot(f_range_tls, df_psd_music.loc[noise_source, obs], label=legend_psd[noise_source]['music'], linestyle='--')
            else:
                ax.plot(f_range, df_psd_paa.loc[noise_source, obs], label=legend_psd[noise_source]['paa'])
                ax.plot(f_range, df_psd_kid.loc[noise_source, obs], label=legend_psd[noise_source]['kid'], linestyle='--')

        # Legend in last plot
        # Collect handles and labels from one axis
        handles, labels = axes[0].get_legend_handles_labels()

        fig.legend(
            handles, labels,
            loc='center left',
            bbox_to_anchor=(0.82, 0.5),  # Move legend closer to the plot
            ncol=1,
            frameon=True
        )
        fig.tight_layout(rect=[0, 0, 0.83, 1])  # Also adjust layout to leave just enough space

        # Save figure
        save_dir = os.path.dirname(f"{plot_dir}_{noise_source}")
        if save_dir:  # avoid error if save_path is just a filename
            os.makedirs(save_dir, exist_ok=True)
        fig.savefig(plot_dir+f"{noise_source}.pdf", dpi=300, bbox_inches='tight')
        fig.savefig(plot_dir+f"{noise_source}.png", dpi=300, bbox_inches='tight')
        plt.close(fig)

def plot_signal_time(plot_dir):
    t = np.linspace(-2*1e-3, 5*1e-3, 1000)  # time array from -2 to 5 s

    # Compute signal
    s_t = s_exponential(t, tau_r_target)

    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(t*1e3, s_t, label=rf"$\tau_r = {tau_r_target*1e3:.1f}\,\mathrm{{ms}}$", color='tab:blue')
    plt.title("Exponential Decay Signal")
    plt.xlabel("Time $t$ [ms]")
    plt.ylabel("$s(t)$")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper right', frameon=True)

        # Save figure
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir+".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir+".png", dpi=300, bbox_inches='tight')

    plt.close()

def compare_resolution(plot_dir):
    df_res = res_all(debug=True)
    # Compute signal
    res_values = df_res.iloc[:, 0].values  # assumes only 1 column

    x = np.arange(len(res_values))  # index array for bar positions
    zeros = np.zeros_like(x)
    labels = df_res.index.tolist()

    # Plot each error bar with its own legend entry
    for i, (lab, val) in enumerate(zip(labels, res_values)):
        plt.errorbar(
            x[i], 0,
            yerr=val*1000,  # Convert to meV
            capsize=8, elinewidth=2, marker='s', markersize=6,
            label=f"{lab}: {val*1000:.2f} meV"
        )
    plt.xticks(x, labels, rotation=45, ha='right')
    plt.ylabel(r"$\sigma_{E_{\mathrm{abs}}}$ [meV]")
    
    y_max = 5 # meV 
    plt.ylim(-y_max, y_max)

    plt.title("Energy Resolution")
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        # Save figure
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir+".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir+".png", dpi=300, bbox_inches='tight')

    plt.close()

def gr_res_vs_eabs(plot_dir):
    """
    Plot energy resolution vs different detector parameters in a 2x3 grid.

    Parameters
    ----------
    plot_dir : str
        Directory where the plot will be saved.
    e_abs : array-like
        Absorbed energy values.
    res : array-like
        Energy resolution values.
    dn_qp : array-like
        Change in quasiparticle number.
    dqi : array-like
        Change in internal quality factor.
    dfr : array-like
        Change in resonant frequency.
    """
    n_x = 2
    n_y = 3
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()  # flatten to 1D array for easier indexing

    # Define variables and labels for loop plotting
    e_abs = np.linspace(2, 500, 1000)  # meV
    e_abs_s21 = [0, 20, 200] # mev 
    zeros = np.zeros_like(e_abs)  # array of 1000 zeros
    ones  = np.ones_like(e_abs)   # array of 1000 ones

    x_labels = [
        r"$E_{\mathrm{abs}}$ (meV)",
        r"$E_{\mathrm{abs}}$ (meV)",
        r"$E_{\mathrm{abs}}$ (meV)",
        r"$E_{\mathrm{abs}}$ (meV)",
        r'Re[$S_{21}$]',
    ]
    y_labels = [
        r"$\sigma_{E_{\mathrm{abs}}}$ (meV)",
        r"$\sigma_{\delta n_{qp}}/\delta n_{qp}(E_{abs})$",
        r"$\sigma_{Q_i}/Q_i(E_{abs})$",
        r"$\sigma_{f_r}/f_r(E_{abs})$",
        r'Im[$S_{21}$]',
    ]

    axs[0].axhline(0, color="k", lw=1.5, linestyle="--", label="Baseline (0)")
    axs[0].fill_between(
        e_abs,
        -gr_res_paa * np.ones_like(e_abs) *1e3,
        +gr_res_paa * np.ones_like(e_abs) *1e3,
        color="tab:orange",
        alpha=0.3,
        label=fr"$\pm {gr_res_paa:.3f}$ band"
    )

    dnqp_paa = eabs_to_dnqp(e_abs*1e-3, delta_0_hf, vol_paa)
    dnqp_paa_err = eabs_to_dnqp(gr_res_paa, delta_0_hf, vol_paa)
    # Fill error band
    axs[1].fill_between(
        e_abs,
        (dnqp_paa - dnqp_paa_err)/dnqp_paa -1,
        (dnqp_paa + dnqp_paa_err)/dnqp_paa -1,
        color="tab:blue",
        alpha=0.3,
        label=rf"$\pm$ {dnqp_paa_err*1e-18:.2f} $\,\mathrm{{(\mu m^{{-3}})}}$"
    )

    qi_paa = eabs_to_qi(e_abs*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom)
    qi_paa_high = eabs_to_qi(e_abs*1e-3+gr_res_paa, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom)
    qi_paa_low = eabs_to_qi(e_abs*1e-3-gr_res_paa, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom)
    axs[2].fill_between(
        e_abs,
        qi_paa_high/qi_paa -1,
        qi_paa_low/qi_paa -1,
        color="tab:blue",
        alpha=0.3
    )

    fr_paa = eabs_to_fr(e_abs*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa)
    fr_paa_high = eabs_to_fr(e_abs*1e-3+gr_res_paa, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa)
    fr_paa_low = eabs_to_fr(e_abs*1e-3-gr_res_paa, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa)
    val_list = fr_paa_high/fr_paa-1
    val_label = val_list[0]
    axs[3].fill_between(
        e_abs,
        fr_paa_high/fr_paa-1,
        fr_paa_low/fr_paa-1,
        color="tab:blue",
        alpha=0.3, label=rf"$\pm{val_label:.3e}$"
    )

    for e_abs in e_abs_s21: 
        s21_paa = s21_ideal_eabs(f_paa, e_abs*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
        s21_paa_high = s21_ideal_eabs(f_paa, e_abs*1e-3+gr_res_paa, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
        s21_paa_low = s21_ideal_eabs(f_paa, e_abs*1e-3-gr_res_paa, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
        s21_paa_fr = s21_ideal_eabs(f_0_nom, e_abs*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
        s21_paa_high_fr = s21_ideal_eabs(f_0_nom, e_abs*1e-3+gr_res_paa, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
        s21_paa_low_fr = s21_ideal_eabs(f_0_nom, e_abs*1e-3-gr_res_paa, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)

                # Extract real and imag parts
        x_nom, y_nom = np.real(s21_paa), np.imag(s21_paa)
        x_low, y_low = np.real(s21_paa_low), np.imag(s21_paa_low)
        x_high, y_high = np.real(s21_paa_high), np.imag(s21_paa_high)

        # Build polygon path: high forward, low reversed
        x_band = np.concatenate([x_high, x_low[::-1]])
        y_band = np.concatenate([y_high, y_low[::-1]])

        axs[4].plot(np.real(s21_paa), np.imag(s21_paa), label=f'{e_abs} meV')
        axs[4].plot(np.real(s21_paa_fr), np.imag(s21_paa_fr), marker='o', markersize=8, linestyle='None')
        axs[4].plot(np.real(s21_paa_high_fr), np.imag(s21_paa_high_fr), marker='^', markersize=8, linestyle='None')
        axs[4].plot(np.real(s21_paa_low_fr), np.imag(s21_paa_low_fr), marker='v', markersize=8, linestyle='None')
        axs[4].fill(x_band, y_band, color="tab:blue", alpha=0.3)
    axs[4].set_aspect('equal', adjustable='box')

    for i, label in enumerate(x_labels):
        axs[i].set_xlabel(label)
        axs[i].set_ylabel(y_labels[i])
        axs[i].grid(True)

    # for ax in [axs[2]]: 
    #     ax.set_yscale('log')
    #     ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    for ax in [axs[0], axs[1], axs[3], axs[4]]: 
        ax.legend()

    plt.tight_layout()
    # Save figure
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir+".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir+".png", dpi=300, bbox_inches='tight')

    plt.close()

def tls_res_vs_eabs(plot_dir):
    """
    Plot TLS-limited energy resolution vs different detector parameters in a 2x3 grid.
    Mirrors gr_res_vs_eabs.

    Parameters
    ----------
    plot_dir : str
        Directory where the plot will be saved.
    """
    n_x = 2
    n_y = 4
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()

    # Define variables
    e_abs_list = np.linspace(0.1, 500, 1000)  # meV
    e_abs_s21 = [0, 20, 200]           # sample points for S21 circles
    colors = ["tab:blue", "tab:orange", "tab:green"]
    dnqp_paa = eabs_to_dnqp(e_abs_list*1e-3, delta_0_hf, vol_paa)

    x_labels = [
        r"$E_{\mathrm{abs}}$ (meV)",
        r"$E_{\mathrm{abs}}$ (meV)",
        r"$E_{\mathrm{abs}}$ (meV)",
        r'Re[$S_{21}$]',
        r'Re[$S_{21}$]',
        r'f (GHz)',
        r"$E_{\mathrm{abs}}$ (meV)",
        r"$E_{\mathrm{abs}}$ (meV)",
    ]
    y_labels = [
        r"$\sigma^{\mathrm{TLS}}_{E_{\mathrm{abs}}}$ (meV)",
        r"$\sigma^{\mathrm{TLS}}_{\delta n_{qp}}/\delta n_{qp}(E_{abs})$",
        r"$\sigma^{\mathrm{TLS}}_{f_r}/f_{r,0}(E_{abs})$",
        r'Im[$S_{21}$]',
        r'Im[$S_{21}$]',
        r'$\theta$ (deg)',
        r'radius',
        r'$\delta \theta/ \theta$',
    ]

    fr_res_dff = tls_variance(tau_r_target, j_dff_tls_paa_1khz, froll_paa, deltaf_paa)
    fr_res_eabs = convert_tls_res_to_eabs_res(fr_res_dff, vol_paa, delta_0_hf, 
        alpha_paa, gamma_nom, k2_paa)
    fr_res_dnqp = convert_tls_res_to_dnqp_res(fr_res_dff, alpha_paa, gamma_nom, k2_paa)
    r_eabs, xc_eabs = s21_circle_radius(e_abs_list*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
    s21_paa_fr = s21_ideal_eabs(f_0_nom, e_abs_list*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
    s21_paa_high_fr = s21_ideal_eabs(f_0_nom, e_abs_list*1e-3, t_eff_hf, f_0_nom*(1+fr_res_dff), delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
    s21_paa_low_fr = s21_ideal_eabs(f_0_nom, e_abs_list*1e-3, t_eff_hf, f_0_nom*(1-fr_res_dff), delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
    theta_eabs = get_unwrapped_phase_deg(s21_paa_fr-xc_eabs)
    theta_eabs_high = get_unwrapped_phase_deg(s21_paa_high_fr-xc_eabs)
    theta_eabs_low = get_unwrapped_phase_deg(s21_paa_low_fr-xc_eabs)

    # --- Panel 0: baseline ± TLS resolution ---
    axs[0].axhline(0, color="k", lw=1.5, linestyle="--")
    axs[0].fill_between(e_abs_list, fr_res_eabs * 1e3, -fr_res_eabs * 1e3, color="tab:orange",
        alpha=0.3, label=rf"$\pm{fr_res_eabs*1e3:.3f}\,\mathrm{{meV}}$")

    # --- Panel 1: dnqp with TLS error band ---
    axs[1].fill_between(e_abs_list,
        -(fr_res_dnqp/dnqp_paa),
        fr_res_dnqp/dnqp_paa,
        color="tab:blue", alpha=0.3,
        label=rf"$\pm$ {fr_res_dnqp*1e-18:.2f} $\,\mathrm{{(\mu m^{{-3}})}}$"
    )

    # --- Panel 3: f_r with TLS error band ---
    axs[2].fill_between(e_abs_list, fr_res_dff, -fr_res_dff, color="tab:blue", 
        alpha=0.3, label=rf"$\pm{fr_res_dff:.3e}$")

    # --- Panel 4: S21 resonance circle(s) with TLS band ---
    plot_ind = axs[3]
    for i, e_abs in enumerate(e_abs_s21):
        s21_tls = s21_ideal_eabs(f_paa, e_abs*1e-3, t_eff_hf, f_0_nom,
                                 delta_0_hf, alpha_gamma_paa, N_0_hf,
                                 vol_paa, qi0_nom, qc0_nom)
        s21_tls_high = s21_ideal_eabs(f_paa, e_abs*1e-3, t_eff_hf,
                                      f_0_nom*(1+fr_res_dff), delta_0_hf, alpha_gamma_paa,
                                      N_0_hf, vol_paa, qi0_nom, qc0_nom)
        s21_tls_low = s21_ideal_eabs(f_paa, e_abs*1e-3, t_eff_hf,
                                     f_0_nom*(1-fr_res_dff), delta_0_hf, alpha_gamma_paa,
                                     N_0_hf, vol_paa, qi0_nom, qc0_nom)
        s21_paa_fr = s21_ideal_eabs(f_0_nom, e_abs*1e-3, t_eff_hf, f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
        s21_paa_high_fr = s21_ideal_eabs(f_0_nom, e_abs*1e-3, t_eff_hf, f_0_nom*(1+fr_res_dff), delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
        s21_paa_low_fr = s21_ideal_eabs(f_0_nom, e_abs*1e-3, t_eff_hf, f_0_nom*(1-fr_res_dff), delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)

        x_nom, y_nom = np.real(s21_tls), np.imag(s21_tls)
        x_low, y_low = np.real(s21_tls_low), np.imag(s21_tls_low)
        x_high, y_high = np.real(s21_tls_high), np.imag(s21_tls_high)

        x_band = np.concatenate([x_high, x_low[::-1]])
        y_band = np.concatenate([y_high, y_low[::-1]])

        plot_ind.plot(np.real(s21_tls), np.imag(s21_tls), '--', label=f'{e_abs} meV', color=colors[i])
        plot_ind.plot(np.real(s21_paa_fr), np.imag(s21_paa_fr), marker='o', markersize=8, linestyle='None')
        plot_ind.plot(np.real(s21_paa_high_fr), np.imag(s21_paa_high_fr), marker='^', markersize=8, linestyle='None')
        plot_ind.plot(np.real(s21_paa_low_fr), np.imag(s21_paa_low_fr), marker='v', markersize=8, linestyle='None')
        plot_ind.fill(x_band, y_band, color="tab:blue", alpha=0.3)

        x = np.real(s21_tls)
        y = np.imag(s21_tls)
        A = np.c_[2*x, 2*y, np.ones_like(x)]
        b = x**2 + y**2
        c, residuals, _, _ = np.linalg.lstsq(A, b, rcond=None)
        x_c, y_c, c0 = c
        radius = np.sqrt(c0 + x_c**2 + y_c**2)
        circle = Circle((x_c, y_c), radius, color='r', fill=False, lw=2, 
            label=fr"$x_c={x_c:.3f},\; y_c={y_c:.3f},\; r={radius:.3f}$")
        axs[3].add_patch(circle)

        # Shift S21 data so that the circle center is at (0,0)
        s21_shifted = (x - x_c) + 1j*(y - y_c)

        # Polar coordinates
        r = np.abs(s21_shifted)
        theta = np.unwrap(np.angle(s21_shifted)) * 180 / np.pi  # unwrap to avoid 2π jumps

        axs[4].plot(np.real(s21_shifted), np.imag(s21_shifted), '--', label=f'{e_abs} meV', color=colors[i])
        axs[4].plot(np.real(s21_paa_fr)-x_c, np.imag(s21_paa_fr)-y_c, marker='o', markersize=8, linestyle='None')
        axs[4].plot(np.real(s21_paa_high_fr)-x_c, np.imag(s21_paa_high_fr)-y_c, marker='^', markersize=8, linestyle='None')
        axs[4].plot(np.real(s21_paa_low_fr)-x_c, np.imag(s21_paa_low_fr)-y_c, marker='v', markersize=8, linestyle='None')

        f_0_nom_theta = get_phase_at_freq(f_0_nom, f_paa, theta)
        f_high_theta = get_phase_at_freq(f_0_nom*(1+fr_res_dff), f_paa, theta)
        f_low_theta = get_phase_at_freq(f_0_nom*(1-fr_res_dff), f_paa, theta)
        axs[5].plot(f_paa*1e-9, theta, color=colors[i])
        axs[5].plot(f_0_nom*(1+fr_res_dff)*1e-9, f_high_theta, marker='^', markersize=8, linestyle='None')
        axs[5].plot(f_0_nom*(1-fr_res_dff)*1e-9, f_low_theta, marker='v', markersize=8, linestyle='None')
        axs[5].plot(f_0_nom*1e-9, f_0_nom_theta, marker='o', markersize=8, linestyle='None')

    # --- Panel 5/6: arc/phase direction ---
    axs[6].plot(e_abs_list, r_eabs)
    # axs[7].plot(e_abs_list, theta_eabs)
    axs[7].fill_between(e_abs_list, 
        (theta_eabs_low-theta_eabs)/theta_eabs, 
        (theta_eabs_high-theta_eabs)/theta_eabs,
        color="tab:blue", alpha=0.3, label=r"$\theta(f_{r,0} \pm \delta f_r)$")

    # --- Axis labels, grids, legends ---
    for i, label in enumerate(x_labels):
        axs[i].set_xlabel(label)
        axs[i].set_ylabel(y_labels[i])
        axs[i].grid(True)

    for ax in [axs[0], axs[1], axs[2], axs[3], axs[4], axs[7]]:
        ax.legend()

    for ax in [axs[3], axs[4]]:
        ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + "_TLS.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + "_TLS.png", dpi=300, bbox_inches='tight')
    plt.close()

def amp_res_vs_eabs(plot_dir, debug=False):
    """
    Plot TLS-limited energy resolution vs different detector parameters in a 2x3 grid.
    Mirrors gr_res_vs_eabs.

    Parameters
    ----------
    plot_dir : str
        Directory where the plot will be saved.
    """
    n_x = 2
    n_y = 2
    fig, axs = plt.subplots(n_x, n_y, figsize=(8*n_y, 6*n_x))
    axs = axs.flatten()

    # Define variables
    e_abs_list = np.linspace(0, 500, 1000)  # meV
    e_abs_s21 = [0, 20, 200]           # sample points for S21 circles
    colors = ["tab:blue", "tab:orange", "tab:green"]

    x_labels = [
        r"$E_{\mathrm{abs}}$ (meV)",
        r"$E_{\mathrm{abs}}$ (meV)",
        r"$E_{\mathrm{abs}}$ (meV)",
        r"$E_{\mathrm{abs}}$ (meV)",
    ]
    y_labels = [
        r"$\delta r$",
        r"$\delta \theta$ (deg)",
        r"$\sigma_{E_{abs},amp}$ (meV)",
        r"$\sigma_{E_{abs},tot}$ (meV)",
    ]

    amp_ds21_res_paa = compute_amp_resolution(tau_r_target, tn_nom, pfeed_paa)
    r_eabs, xc_eabs = s21_circle_radius(e_abs_list*1e-3, t_eff_hf, 
        f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
    s21_paa_fr = s21_ideal_eabs(f_0_nom, e_abs_list*1e-3, t_eff_hf, 
        f_0_nom, delta_0_hf, alpha_gamma_paa, N_0_hf, vol_paa, qi0_nom, qc0_nom)
    theta_eabs = get_unwrapped_phase_deg(s21_paa_fr-xc_eabs)
    d_theta = angle_diff_from_im_shift(r_eabs[0], amp_ds21_res_paa)

    if debug:
        print("[DEBUG] Parameters used for d_theta calculation:")
        print(f"  tau_r_target        = {tau_r_target}")
        print(f"  tn_nom              = {tn_nom}")
        print(f"  pfeed_paa           = {pfeed_paa}")
        print(f"  e_abs_list[0] (mJ)     = {e_abs_list[0]}")
        print(f"  t_eff_hf            = {t_eff_hf}")
        print(f"  f_0_nom             = {f_0_nom}")
        print(f"  delta_0_hf          = {delta_0_hf}")
        print(f"  alpha_gamma_paa     = {alpha_gamma_paa}")
        print(f"  N_0_hf              = {N_0_hf}")
        print(f"  vol_paa             = {vol_paa}")
        print(f"  qi0_nom             = {qi0_nom}")
        print(f"  qc0_nom             = {qc0_nom}")
        print(f"  r_eabs[0]           = {r_eabs[0]}")
        print(f"  amp_ds21_res_paa    = {amp_ds21_res_paa}")
        print(f"  xc_eabs             = {xc_eabs}")
        print(f"  s21_paa_fr          = {s21_paa_fr}")
        print(f"  theta_eabs          = {theta_eabs}")

    d_eabs_p_diss, d_eabs_m_diss = intercept_energy_offset(e_abs_list, r_eabs, amp_ds21_res_paa)
    d_eabs_p_freq, d_eabs_m_freq = intercept_energy_offset(e_abs_list, theta_eabs, d_theta)
    d_eabs_diss_tot = compute_total_resolution_list([d_eabs_p_diss, gr_paa*1e3])
    d_eabs_freq_tot = compute_total_resolution_list([d_eabs_p_freq, tls_eabs_paa*1e3, gr_paa*1e3])

    # --- Panel 4: S21 resonance circle(s) with TLS band ---
    # axs[0].plot(e_abs_list, r_eabs)
    axs[0].fill_between(e_abs_list, - amp_ds21_res_paa, + amp_ds21_res_paa, 
        color="tab:blue", alpha=0.3, label=rf"$\pm{amp_ds21_res_paa:.3e}$")
    # axs[1].plot(e_abs_list, theta_eabs)
    axs[1].fill_between(e_abs_list, - d_theta, + d_theta, 
        color="tab:blue", alpha=0.3, label=rf"$\pm{d_theta:.3e}$")

    # --- Panel 0: baseline ± TLS resolution ---
    axs[2].axhline(0, color="k", lw=1.5, linestyle="--")
    max_idx_diss = np.argmax(d_eabs_m_diss)
    axs[2].fill_between(e_abs_list, d_eabs_p_diss, d_eabs_m_diss, color="tab:orange",
        alpha=0.3, label=rf'diss-max@{d_eabs_m_diss[max_idx_diss]:.2e} meV')
    max_idx_freq = np.argmax(d_eabs_p_freq)
    axs[2].fill_between(e_abs_list, d_eabs_p_freq, d_eabs_m_freq, color="tab:blue",
    alpha=0.3, label=rf'freq-max@{d_eabs_p_freq[max_idx_freq]:.2e} meV')

    axs[3].axhline(0, color="k", lw=1.5, linestyle="--")
    max_idx_diss = np.argmax(d_eabs_diss_tot)
    axs[3].fill_between(e_abs_list, d_eabs_diss_tot, -d_eabs_diss_tot, color="tab:orange",
        alpha=0.3, label=rf'diss-max@{d_eabs_diss_tot[max_idx_diss]:.2e} meV')
    max_idx_freq = np.argmax(d_eabs_freq_tot)
    axs[3].fill_between(e_abs_list, d_eabs_freq_tot, -d_eabs_freq_tot, color="tab:blue",
    alpha=0.3, label=rf'freq-max@{d_eabs_freq_tot[max_idx_freq]:.2e} meV')

    # --- Axis labels, grids, legends ---
    for i, label in enumerate(x_labels):
        axs[i].set_xlabel(label)
        axs[i].set_ylabel(y_labels[i])
        axs[i].grid(True)

    for ax in [axs[0], axs[1], axs[2], axs[3]]:
        ax.legend()

    plt.tight_layout()
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(plot_dir + ".pdf", dpi=300, bbox_inches='tight')
    plt.savefig(plot_dir + ".png", dpi=300, bbox_inches='tight')
    plt.close()

