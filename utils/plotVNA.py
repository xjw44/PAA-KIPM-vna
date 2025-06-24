import numpy as np
import matplotlib.pyplot as plt
import os

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

def overlay_traces(file_list, file_leg, plt_title, scan_range, save_path):
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
    # plt.ylim(-1.005, -0.995)  # Adjust if needed
    plt.tight_layout()
    save_dir = os.path.dirname(save_path)
    if save_dir:  # avoid error if save_path is just a filename
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path+'.pdf', dpi=300)
    plt.savefig(save_path+'.png', dpi=300)
    plt.show()
    return freq_list, s21_list