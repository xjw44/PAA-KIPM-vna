from utils.readFiles import *

# # ######################################### 7/10/2025 find the qc, fr, mb
folder_son = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-7-10-son-paa-hf/'
output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-7-10-hf-m-b/'

son_repstrange = folder_son+'2025-7-10-PAAKID_feedline-copied-from-Shilin-s21'+'.csv'
freq, s21_real, s21_imag, s21_mag, s21_z = read_s21_son(son_repstrange)

s21_list = [s21_mag]
freq_list = [freq]
z_list = [s21_real+1j*s21_imag]
file_list_leg = ['sonnet_sim'] 
title = 'al-hf-son'
scan_range = [4.0382, 4.0384] # ghz
fr_est = 4.0383 # ghz 

masked_freq_list, masked_z_list = mask_frequency_range_list(freq_list, z_list, scan_range[0], scan_range[1])

son_files = {"s21_list": s21_list, 
"freq_list": freq_list,
"z_list": z_list,
"masked_z_list": masked_z_list,
"masked_freq_list": masked_freq_list,
"file_list_leg": file_list_leg, 
"title": title,
"scan_range": scan_range,
"output_path": output_path,
"fr_est": fr_est
}