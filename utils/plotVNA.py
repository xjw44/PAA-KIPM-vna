import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import scipy.signal as sig

def load_noise_data(filename):
    freqs = []
    s21_values = []

    with open(filename, 'r') as f:
        for line in f:
            # Skip comment or metadata lines
            if line.strip().startswith('!'):
                continue

            # Split the line into parts
            parts = line.strip().split()
            if len(parts) >= 2:
                freq = float(parts[0].split(",")[0])
                s21 = float(parts[1])
                freqs.append(freq)
                s21_values.append(s21)
                # except ValueError:
                #     continue  # Skip lines with non-numeric data

    return np.array(freqs), np.array(s21_values)

def overlay_traces(file_list, file_leg, plt_title, scan_range, scan_range_y, save_path):
    plt.figure(figsize=(8, 4))

    freq_list = []
    s21_list = []
    for i, file in enumerate(file_list):
        freqs, s21 = load_noise_data(file)
        freqs_mhz = freqs / 1e6
        plt.plot(freqs_mhz, s21, '-', label=file_leg[i])
        freq_list.append(freqs)
        s21_list.append(s21)

    plt.xlabel('Frequency (MHz)')
    plt.ylabel('S21 (dB)')
    plt.title(plt_title)
    plt.grid(True)
    plt.legend()
    plt.xlim(scan_range[0], scan_range[1])
    if scan_range_y: 
        plt.ylim(scan_range_y[0], scan_range_y[1])  # Adjust if needed
    plt.tight_layout()
    save_dir = os.path.dirname(save_path)
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path+'.pdf', dpi=300)
    plt.savefig(save_path+'.png', dpi=300)
    plt.show()
    plt.close()
    return freq_list, s21_list

def scan_peaks(file, peak_list):
    freqs, s21 = load_noise_data(file)  # assuming arrays # hz 
    # Convert to DataFrame
    df = pd.DataFrame({'frequency': freqs, 's21': s21})

    range_x_list = []
    range_y_list = []
    for i, peak in enumerate(peak_list):
        range_x_min = (freqs[peak] - 1e6)
        range_x_max = (freqs[peak] + 1e6)

        # Filter within frequency range
        df_range = df[(df['frequency'] >= range_x_min) & (df['frequency'] <= range_x_max)]

        # Get min and max s21 in that range
        s21_min = df_range['s21'].min()
        s21_max = df_range['s21'].max()
        range_x_list.append([range_x_min*1e-6, range_x_max*1e-6]) # mhz
        range_y_list.append([s21_min, s21_max])

    return range_x_list, range_y_list


def scan_window(file, save_path, title):
    freqs, s21 = load_noise_data(file)  # assuming arrays # hz 
    # Convert to DataFrame
    df = pd.DataFrame({'frequency': freqs, 's21': s21})

    step_size = 20 
    window_size = 100 

    for start in range(0, len(df) - window_size + 1, step_size):
        end = start + window_size
        freqs_win = freqs[start:end]*1e-6 # mhz 
        s21_win = s21[start:end]

        # Polynomial fit
        z = np.polyfit(freqs_win, s21_win, 3)
        p = np.poly1d(z)

        y_fit = p(freqs_win)
        residuals = s21_win - y_fit

        # Detect peaks in original trace and in residual
        peaklist, _ = sig.find_peaks(freqs_win)
        peaklist_res, properties_res = sig.find_peaks(residuals, width=[10, 200], prominence=[0, 0.02])

        # Create subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 6), gridspec_kw={'height_ratios': [1, 1]})

        # --- Top: trace + fits + peaks ---
        ax1.plot(freqs_win, s21_win, label='Filtered magnitude (mfz)')
        ax1.plot(freqs_win, y_fit, 'r--', label='100th-Degree Poly Fit')

        for i in peaklist:
            ax1.axvline(freqs_win[i], color='red', linestyle='--', alpha=0.3)

        ax1.set_ylabel('|S21| (dB)')
        ax1.set_title(f'Window {start}-{end}')
        ax1.grid(True)
        ax1.legend()

        # --- Bottom: residuals + residual peaks ---
        ax2.plot(freqs_win, residuals, label='Residual (100 deg)')
        ax2.axhline(0, color='k', linestyle='--', linewidth=1)
        ax2.set_ylim(-0.02, 0.02)

        for i in peaklist_res:
            ax2.axvline(freqs_win[i], color='red', linestyle='--', alpha=0.5)

        ax2.set_xlabel('Frequency (MHz)')
        ax2.set_ylabel('Residual')
        ax2.grid(True)
        ax2.legend()

        plt.tight_layout()

        # Save
        save_dir = os.path.dirname(save_path)
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)

        title_suf = f"{title}_start{start}_end{end}"
        plt.savefig(save_path + title_suf + '.pdf', dpi=300)
        plt.savefig(save_path + title_suf + '.png', dpi=300)
        plt.close()

        print(f"Saved window {start}-{end} with {len(peaklist_res)} residual peaks.")


    