from utils.readFiles import *

# # ######################################### 7/10/2025 find the qc, fr, mb
folder_son = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-pcb-fin/data/2025-7-14-fin-cond-son/'
output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-pcb-fin/output/2025-7-14-pcb-cond-fin/'

# son_files_raw = {
# "lossless_ideal": folder_son+'2025-7-2-paa-lossless_idealcap_rep_s21'+'.csv',
# "lossless": folder_son+'2025-7-2-paa-lossless_s21'+'.csv',
# }

son_files_raw = {
# "lossless_ideal": folder_son+'2025-7-2-paa-lossless_idealcap_rep_s21'+'.csv',
# "lossless": folder_son+'2025-7-2-paa-lossless_s21'+'.csv',
"*1000cond": folder_son+'2025-7-2-paa-1000cond_s21'+'.csv',
"*100cond": folder_son+'2025-7-2-paa-100cond_s21'+'.csv',
# "*1cond": folder_son+'2025-6-26-paa-fr4-bridge_s21_primary'+'.csv',
# "*0.1cond": folder_son+'2025-7-2-paa-0.1cond_s21'+'.csv',
}

res_files = {
"lossless_ideal": [536.4*1e6, 571.4*1e6, 766.8*1e6],
"lossless": [340.1*1e6, 481.1*1e6, 522.6*1e6, 568.5*1e6, 650.2*1e6, 735.1*1e6],
}

s21mag_list = []
s21z_list = []
freq_list = []
file_list_leg = []
for label, filepath in son_files_raw.items():
    print(f"Processing '{label}' from file: {filepath}")
    freq, s21_real, s21_imag, s21_mag, s21_z = read_s21_son(filepath)
    s21mag_list.append(s21_mag)
    s21z_list.append(s21_real+1j*s21_imag)
    freq_list.append(freq*1e6) # hz
    file_list_leg.append(label)

# title = 'ideal_real_capacitor'
title = 'vary_cond'
# scan_range = [200*1e6, 800*1e6] # hz 
scan_range = [339.8*1e6, 340.3*1e6] # hz 
scan_range_y = [-2, 0] # db
fr_est = 340*1e6 # hz

masked_freq_list, masked_s21z_list = mask_frequency_range_list(freq_list, s21z_list, scan_range[0], scan_range[1])

son_files = {"s21mag_list": s21mag_list, 
"freq_list": freq_list,
"s21z_list": s21z_list,
"masked_s21z_list": masked_s21z_list,
"masked_freq_list": masked_freq_list,
"res_list": list(res_files.values()),
"file_list_leg": file_list_leg, 
"title": title,
"output_path": output_path,
"scan_range": scan_range,
"scan_range_y": scan_range_y,
"fr_est": fr_est,
}