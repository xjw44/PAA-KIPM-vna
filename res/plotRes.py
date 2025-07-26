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
sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/res')
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
    'music': "\n".join(["music_freq@" + + label_l_paa, label_t_ind_paa, label_w_ind_paa, 
        label_tn_nom, label_rho_nom_al, label_pfeed_music, 
        label_delta_0_al, label_vol_music, label_alpha_music, label_gamma_nom, 
        label_qi0_music, label_qc0_music, label_qr0_music, label_f_0_nom, label_t_eff_hf, 
        label_N_0_hf, label_k1_paa, label_k2_paa]),
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
    j_eabs_gr_paa = J_GR_eabs(f_range, nqp_target, tau_r_target, vol_paa, delta_0_hf)
    j_eabs_gr_kid = J_GR_eabs(f_range, nqp_target, tau_r_target, vol_kid, delta_0_al)
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
    j_amp_rest_paa = amp_psd_all(j_amp_ds21_paa, qr0_nom, qc0_nom, vol_paa, alpha_paa, gamma_nom, k1_paa, k2_paa, delta_0_hf)
    j_amp_rest_kid = amp_psd_all(j_amp_ds21_kid, qr0_nom, qc0_nom, vol_kid, alpha_kid, gamma_nom, k1_kid, k2_kid, delta_0_al)

    j_dff_tls_music = full_psd_tls(f_range_tls, j_dff_tls_music_1khz, froll_music, tls_n)
    j_dff_tls_paa_1khz = update_tls_psd(t_eff_hf, t_eff_music, v_c_paa, v_c_music, 
        eres_paa, eres_music, j_dff_tls_music_1khz, tls_beta, debug=True)
    j_dff_tls_paa = full_psd_tls(f_range_tls, j_dff_tls_paa_1khz, froll_paa, tls_n)
    j_tls_rest_music = convert_dff_tls_psd_to_all(j_dff_tls_music, vol_music, alpha_music, gamma_nom, 
        k2_music, delta_0_al, qr0_music, qc0_music)
    j_tls_rest_paa = convert_dff_tls_psd_to_all(j_dff_tls_paa, vol_paa, alpha_paa, gamma_nom, 
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
    t = np.linspace(-2, 5, 1000)  # time array from -2 to 5 s

    # Compute signal
    s_t = causal_exponential(t, tau_r_target)

    # Plot
    plt.figure(figsize=(8, 4))
    plt.plot(t, s_t, label=rf"$s(t) = e^{{-t/\tau}},\ \tau = {tau_r_target*1e3:.1f} ms$", color='tab:blue')
    plt.title("Exponential Decay Signal")
    plt.xlabel("Time $t$ [s]")
    plt.ylabel("$s(t)$")
    plt.grid(True, linestyle='--', alpha=0.6)

        # Save figure
    save_dir = os.path.dirname(f"{plot_dir}")
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    fig.savefig(plot_dir+".pdf", dpi=300, bbox_inches='tight')
    fig.savefig(plot_dir+".png", dpi=300, bbox_inches='tight')
    plt.close(fig)





