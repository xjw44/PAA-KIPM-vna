import numpy as np
import os
import pandas as pd

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)  # For long string cells

def read_s21_son(filename, ghz=True):
    # Define column names manually
    if ghz: 
        freq_unit = "GHz"
    col_names = [
        f"Frequency ({freq_unit})",
        "RE[S11]", "IM[S11]",
        "RE[S12]", "IM[S12]",
        "RE[S21]", "IM[S21]",
        "RE[S22]", "IM[S22]"
    ]

    # Read CSV with correct header and column names
    df = pd.read_csv(filename, skiprows=3, names=col_names, engine='python')
    # Convert all columns to float safely
    df = df.apply(pd.to_numeric, errors='coerce')  # invalid entries become NaN

    # Drop rows with any NaNs (i.e., corrupted rows like "R 50.00000")
    df = df.dropna()

    # Sort in descending order of frequency
    df = df.sort_values(by=f"Frequency ({freq_unit})", ascending=True).reset_index(drop=True)

    # Extract frequency, real(S21), imag(S21)
    freq = df[f'Frequency ({freq_unit})'].values
    s21_real = df['RE[S21]'].values
    s21_imag = df['IM[S21]'].values
    s21_mag = 20*np.log10(np.sqrt(s21_real**2 + s21_imag**2))
    s21_z = s21_real+1j*s21_imag

    return freq, s21_real, s21_imag, s21_mag, s21_z

def mask_frequency_range_list(freq_list, s21_list, fmin, fmax):
    """
    Masks multiple frequency and s21 traces based on a frequency window.

    Parameters:
    - freq_list: list of numpy arrays (frequencies for each trace)
    - s21_list: list of numpy arrays (S21 values for each trace)
    - fmin: minimum frequency (inclusive)
    - fmax: maximum frequency (inclusive)

    Returns:
    - masked_freqs: list of masked frequency arrays
    - masked_s21s: list of masked s21 arrays
    """
    masked_freqs = []
    masked_s21s = []
    
    for freqs, s21s in zip(freq_list, s21_list):
        freqs = np.array(freqs)
        s21s = np.array(s21s)
        mask = (freqs >= fmin) & (freqs <= fmax)
        masked_freqs.append(freqs[mask])
        masked_s21s.append(s21s[mask])
    
    return masked_freqs, masked_s21s