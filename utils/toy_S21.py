import numpy as np
import matplotlib.pyplot as plt
import os

def s21_response(f, fr, Qc, Qi):
    Qr = 1 / (1/Qc + 1/Qi)
    x = (f - fr) / fr
    s21 = 1 - (Qr / Qc) / (1 + 2j * Qr * x)
    s21_db = 20 * np.log10(np.sqrt(s21.real**2+s21.imag**2))
    return s21_db

def plot_s21_magnitude(save_path, fr=90e6, Qc=10, Qi_ratios=[1, 10], f_span=0.3):
    """
    Plot |S21(f)| for varying Qi = ratio * Qc, same Qc and fr

    Parameters:
    - fr: Resonance frequency [Hz]
    - Qc: Coupling quality factor
    - Qi_ratios: List of Qi multipliers (Qi = ratio * Qc)
    - f_span: Frequency span as a fraction of fr (e.g. 0.02 = ±1%)
    """
    f = np.linspace(fr * (1 - f_span), fr * (1 + f_span), 1000)

    plt.figure(figsize=(7, 4))
    for ratio in Qi_ratios:
        Qi = Qc * ratio
        s21 = s21_response(f, fr, Qc, Qi)
        label = rf'$Q_i = {ratio:.0f} \times Q_c$'
        plt.plot(f / 1e6, s21, label=label)

    plt.xlabel('Frequency (MHz)')
    plt.ylabel(r'$S_{21}(f) (dB)$')
    plt.title(rf'$Q_c = {Qc:.0f},\ f_r = {fr/1e6:.2f}\,\mathrm{{MHz}}$')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    save_dir = os.path.dirname(save_path)
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path+'.pdf', dpi=300)
    plt.savefig(save_path+'.png', dpi=300)
    plt.show()
    plt.close()

def plot_s21_magnitude(save_path, fr=340e6, Qi=1e4, Qc_ratios=[1, 5, 10], f_span=0.003):
    """
    Plot |S21(f)| for varying Qc = ratio * Qi, fixed Qi and fr

    Parameters:
    - fr: Resonance frequency [Hz]
    - Qi: Internal quality factor
    - Qc_ratios: List of Qc multipliers (Qc = ratio * Qi)
    - f_span: Frequency span as a fraction of fr (e.g. 0.02 = ±1%)
    """
    f = np.linspace(fr * (1 - f_span), fr * (1 + f_span), 1000)

    plt.figure(figsize=(7, 4))
    for ratio in Qc_ratios:
        Qc = Qi * ratio
        s21 = s21_response(f, fr, Qc, Qi)
        label = rf'$Q_c = {ratio:.0f} \times Q_i$'
        plt.plot(f / 1e6, s21, label=label)

    plt.xlabel('Frequency (MHz)')
    plt.ylabel(r'$S_{21}(f) \, (\mathrm{dB})$')
    # plt.xlim(325, 350) # mhz
    plt.title(rf'$Q_i = {Qi:.0f},\ f_r = {fr/1e6:.2f}\,\mathrm{{MHz}}$')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    save_dir = os.path.dirname(save_path)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path + '.pdf', dpi=300)
    plt.savefig(save_path + '.png', dpi=300)
    plt.show()
    plt.close()


