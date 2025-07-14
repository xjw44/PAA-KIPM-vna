import sys
sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna')
from utils.plotVNA import *
from utils.findPeak import *
from utils.fitres import *
from utils.toy_S21 import *

# # resonance found 
# folder = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-23-pcb-cold/'
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-6-24-pcb-feedline/'

# folder_feedline = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-24-pcb-testmore/'
# # b1_cold = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_resfound60mhz'+'.csv'
# # b1_warm = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_resfound60mhz'+'.csv'
# # b2_warm = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B2_ROOM_resmaybe_ZOOM'+'.csv'
# # b2_cold = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B2_cold_resmaybe_ZOOM'+'.csv'
# # b3_cold = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B3_cold_resmaybe_ZOOM'+'.csv'

# b1_warm = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_room-test2'+'.csv'
# b1_cold = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_cold-test2'+'.csv'

# b1_feed = folder_feedline+'2025_6_24_xjw_noise_trace_ifbw_1000_hz_3dbm_feedline_res'+'.csv'
# # b1_warm_short = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_room_shortcable-nores'+'.csv'
# # b1_warm_mid = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_room_midcable-nores'+'.csv'
# # b1_warm_long_unstable = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_room_longcable-unstableres'+'.csv'
# # b2_warm = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B2_room-test2'+'.csv'
# # b2_cold = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B2_cold-test2'+'.csv'
# # b3_warm = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B3_room-test2'+'.csv'
# # b3_cold = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B3_cold-test2'+'.csv'

# b1_cold_long = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_4_cold_b1'+'.csv'
# b1_cold_real = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_10_cold_b1_real'+'.csv'
# b1_cold_imag = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_10_cold_b1_imag'+'.csv'

# file_list = [b1_warm, b1_cold, b1_feed] 
# file_list_leg = ['b1', 'b1_cold', 'b1_feed'] 
# title = 'feedline_res'
# scan_range = [50, 80] # mhz

# freq_list, s21_list = overlay_traces(file_list, file_list_leg, title, scan_range, False, output_path+title)

# file_list = [b1_cold_long]
# file_list_leg = ['b1_cold_long'] 
# title = 'finding_non_feedline_resonance'
# scan_range = [10, 300] # mhz

# freq_list, s21_list = overlay_traces(file_list, file_list_leg, title, scan_range, False, output_path+title)

###########################################
###########################################
###########################################

# ######################################### 6/25/2025 noise trace why noise so large?
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-6-25-pcb-peakfinding/'
# b1_room_noise = folder_feedline+'2025_6_24_xjw_noise_trace_ifbw_100_hz_3dbm_NOISE_0.001db_ave_10_room'+'.csv'
# file_list = [b1_room_noise]
# file_list_leg = ['b1_room_noise'] 
# title = 'noise_level_room'
# scan_range = [0.3, 0.31] # mhz

# freq_list, s21_list = overlay_traces(file_list, file_list_leg, title, scan_range, False, output_path+title)

# file_list = [b1_cold_long]
# file_list_leg = ['b1_cold_long'] 
# title = 'noise_level_cold'
# scan_range = [10, 10.002] # mhz
# scan_range_y = [-1.65175, -1.6515] # mhz

# freq_list, s21_list = overlay_traces(file_list, file_list_leg, title, scan_range, scan_range_y, output_path+title)

# file_list_real = [b1_cold_real]
# file_list_leg_re = ['b1_cold_real'] 
# title = 'finding_non_feedline_resonance_real'
# freq_list, s21_real = overlay_traces(file_list_real, file_list_leg_re, title, scan_range, False, output_path+title)
# file_list_imag = [b1_cold_imag]
# file_list_leg_im = ['b1_cold_imag'] 
# title = 'finding_non_feedline_resonance_imag'
# freq_list, s21_imag = overlay_traces(file_list_imag, file_list_leg_im, title, scan_range, False, output_path+title)

# # print(s21_real, s21_imag, s21_list)

# # title = "filter_trace_width_l10"
# # for index, pcb_freq in enumerate(freq_list): 
# #     print('index: ', index)
# #     freqs_ghz = pcb_freq/1e9
# #     z = s21_real[index] + 1j*s21_imag[index]
# #     plot_range = [10, 300]
# #     plot_range_y = False
# #     peak_res = plot_filtered_trace_with_peaks(freqs_ghz, s21_list[index], plot_range, plot_range_y, output_path+title, fwindow=500*1e-6, start_f=60*1e-3, stop_f=None, nsig=0.05)

# # range_x_list, range_y_list = scan_peaks(b1_cold_long, peak_res)
# # print(range_x_list)

# # for index, range_x in enumerate(range_x_list):
# #     overlay_traces([b1_cold_long], file_list_leg, title+f'_{index}', range_x, range_y_list[index], output_path+title+f'_{index}')

# title = "filter_trace_width_allscans"
# scan_window(b1_cold_long, output_path, title)

###########################################
###########################################
###########################################

######################################### 6/26/2025 noise resonance found sth!! 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-6-26-pcb-bridge/'
# folder_feedline = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-24-pcb-testmore/'
# folder = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-23-pcb-cold/'
# folder_more_feedline = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-26-pcb-feedline/'
# folder_bridge = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-26-pcb-bridge/'

# b1_warm = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_room-test2'+'.csv'
# b1_cold = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_cold-test2'+'.csv'
# b1_warm_short = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_room_shortcable-nores'+'.csv'
# b1_warm_mid = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_room_midcable-nores'+'.csv'
# b1_warm_long_unstable = folder+'2025_6_23_xjw_noise_trace_ifbw_1000_hz_3dbm_B1_room_longcable-unstableres'+'.csv'
# b1_warm_long_cableonly = folder_more_feedline+'2025_6_26_cable_room_ifbw300hz_3dbm_ave10_long_mag'+'.csv'

# b1_feed = folder_feedline+'2025_6_24_xjw_noise_trace_ifbw_1000_hz_3dbm_feedline_res'+'.csv'

# file_list = [b1_warm, b1_warm_long_unstable, b1_feed, b1_warm_short, b1_warm_mid, b1_warm_long_cableonly]
# file_list_leg = ['b1_warm_long', 'b1_warm_long_unstable', 'b1_warm_long_feedonly', 'b1_warm_short', 'b1_warm_mid', 'b1_warm_long_cableonly'] 
# title = 'feedline_resonance'
# scan_range = [50, 80] # mhz

# freq_list, s21_list = overlay_traces(file_list, file_list_leg, title, scan_range, False, output_path+title)
# print("plot done! Feedline resonance cable dependent! ")

# file_list_nofeed = [b1_warm, b1_cold, b1_warm_short, b1_warm_mid, b1_warm_long_unstable]
# file_list_nofeed_leg = ['b1_warm', 'b1_cold', 'b1_warm_short', 'b1_warm_mid', 'b1_warm_long_unstable'] 
# title_nofeed = 'cable_depedent_resonance'
# scan_range = [50, 80] # mhz

# freq_list, s21_list = overlay_traces(file_list_nofeed, file_list_nofeed_leg, title_nofeed, scan_range, False, output_path+title_nofeed)
# print("plot done! Feedline resonance cable dependent! ")

# file_list = [b1_warm, b1_cold, b1_warm_short, b1_warm_mid, b1_warm_long_unstable, b1_feed]
# file_list_leg = ['b1_warm_long', 'b1_cold_long', 'b1_warm_short', 'b1_warm_mid', 'b1_warm_long_unstable', 'b1_warm_long_feedonly'] 
# title = 'cable_depedent_resonance_feedlineonly'
# scan_range = [50, 80] # mhz

# freq_list, s21_list = overlay_traces(file_list, file_list_leg, title, scan_range, False, output_path+title)
# print("plot done! Feedline resonance! ")

# # b1_room_noise = folder_feedline+'2025_6_24_xjw_noise_trace_ifbw_100_hz_3dbm_NOISE_0.001db_ave_10_room'+'.csv'

# feedline_imag = folder_more_feedline+'2025_6_26_feedline_room_ifbw300hz_3dbm_ave5_long_imag'+'.csv'
# feedline_real = folder_more_feedline+'2025_6_26_feedline_room_ifbw300hz_3dbm_ave5_long_real'+'.csv'

# b3_imag = folder_more_feedline+'2025_6_26_b3_room_ifbw300hz_3dbm_ave5_long_imag'+'.csv'
# b3_real = folder_more_feedline+'2025_6_26_b3_room_ifbw300hz_3dbm_ave5_long_real'+'.csv'
# title = 'junk'

# file_list_real = [feedline_real, b3_real]
# file_list_imag = [feedline_imag, b3_imag]
# file_list_leg = ['junk', 'junk'] 

# freq_list, s21_imag_feedline = overlay_traces(file_list_imag, file_list_leg, title, scan_range, False, output_path+title)
# freq_list, s21_real_feedline = overlay_traces(file_list_real, file_list_leg, title, scan_range, False, output_path+title)
# s21_feedline = s21_real_feedline[0] + 1j*s21_imag_feedline[0]
# s21_b3 = s21_real_feedline[1] + 1j*s21_imag_feedline[1]
# s21_list = [s21_feedline, s21_b3]
# file_leg = ['feedline_only_long_warm', 'b3_long_warm']
# title = 'S21_circles_feedline_b3'
# overlay_smith(s21_list, file_leg, title, output_path+title)
# print("plot done! Feedline resonance! ")

# b1_cold_long = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_4_cold_b1'+'.csv'
# b1_cold_real = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_10_cold_b1_real'+'.csv'
# b1_cold_imag = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_10_cold_b1_imag'+'.csv'

# file_list = [b1_cold_long]
# file_list_leg = ['b1_cold_long'] 
# title = 'b1_cold_longscan'
# scan_range = [10, 290] # mhz
# scan_range_y = [-3, -1.5] # mhz

# freq_list, s21_list = overlay_traces(file_list, file_list_leg, title, scan_range, scan_range_y, output_path+title)
# print("plot done! Feedline resonance! ")

# title = 'junk'
# file_list_real = [b1_cold_imag]
# file_list_imag = [b1_cold_real]
# file_list_leg = ['junk'] 

# freq_list, s21_imag_b1_cold = overlay_traces(file_list_imag, file_list_leg, title, scan_range, False, output_path+title)
# freq_list, s21_real_b1_cold = overlay_traces(file_list_real, file_list_leg, title, scan_range, False, output_path+title)
# s21_b1_cold = s21_real_b1_cold[0] + 1j*s21_imag_b1_cold[0]
# s21_list = [s21_b1_cold]
# file_leg = ['b1_cold_long']
# title = 'S21_circles_b1_cold'
# overlay_smith(s21_list, file_leg, title, output_path+title)
# print("plot done! Feedline resonance! ")

# # title = 'osmond_sweepfit'
# # freq_ghz = freq_list[0]*1e-9
# # sweep_fit(freq_ghz, s21_b1_cold, nsig=3, fwindow=5e-4, pdf_rewrite=True, additions=[], filename=output_path+title, start_f=1*1e6*1e-9, stop_f=300*1e6*1e-9)

# b2_warm_imag = folder_bridge+'30mhz_200mhz_shorted_notape_imag'+'.csv'
# b2_warm_real = folder_bridge+'30mhz_200mhz_shorted_notape_real'+'.csv'
# b2_warm_imag_tape1 = folder_bridge+'30mhz_200mhz_shorted_tape1_imag'+'.csv'
# b2_warm_real_tape1 = folder_bridge+'30mhz_200mhz_shorted_tape1_real'+'.csv'
# b2_warm_imag_tape2 = folder_bridge+'30mhz_200mhz_shorted_tape2_imag'+'.csv'
# b2_warm_real_tape2 = folder_bridge+'30mhz_200mhz_shorted_tape2_real'+'.csv'
# b2_warm_imag_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_imag'+'.csv'
# b2_warm_real_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_real'+'.csv'
# b2_warm_imag_tape4 = folder_bridge+'30mhz_200mhz_shorted_tape4_imag'+'.csv'
# b2_warm_real_tape4 = folder_bridge+'30mhz_200mhz_shorted_tape4_real'+'.csv'

# file_list_imag = [b2_warm_imag, b2_warm_imag_tape1, b2_warm_imag_tape2, b2_warm_imag_tape3, b2_warm_imag_tape4]
# file_list_real = [b2_warm_real, b2_warm_real_tape1, b2_warm_real_tape2, b2_warm_real_tape3, b2_warm_real_tape4]
# file_list_leg = ['junk', 'junk', 'junk', 'junk', 'junk'] 
# title = 'kid_bridge_tape_fr_change'
# scan_range = [30, 200] # mhz
# title = 'junk'

# freq_list, s21_imag_tape = overlay_traces(file_list_imag, file_list_leg, title, scan_range, False, output_path+title)
# freq_list, s21_real_tape = overlay_traces(file_list_real, file_list_leg, title, scan_range, False, output_path+title)
# s21_tapeno = 20*np.log10(np.sqrt(s21_real_tape[0]**2 + s21_imag_tape[0]**2))
# s21_tape1 = 20*np.log10(np.sqrt(s21_real_tape[1]**2 + s21_imag_tape[1]**2))
# s21_tape2 = 20*np.log10(np.sqrt(s21_real_tape[2]**2 + s21_imag_tape[2]**2))
# s21_tape3 = 20*np.log10(np.sqrt(s21_real_tape[3]**2 + s21_imag_tape[3]**2))
# s21_tape4 = 20*np.log10(np.sqrt(s21_real_tape[4]**2 + s21_imag_tape[4]**2))
# s21_list = [s21_tapeno, s21_tape1, s21_tape2, s21_tape3, s21_tape4]
# file_list_leg = ['b2_notape_shortcable', 'b2_tape1_shortcable', 'b2_tape2_shortcable', 'b2_tape3_shortcable', 'b2_tape4_shortcable'] 
# title = 's21_tape_fr_shift_shorted'
# overlay_s21mag(freq_list[0], s21_list, file_list_leg, title, scan_range, False, output_path+title)
# print("plot done! Feedline resonance cable dependent! ")

# s21_tapeno = s21_real_tape[0] + 1j*s21_imag_tape[0]
# s21_tape1 = s21_real_tape[1] + 1j*s21_imag_tape[1]
# s21_tape2 = s21_real_tape[2] + 1j*s21_imag_tape[2]
# s21_tape3 = s21_real_tape[3] + 1j*s21_imag_tape[3]
# s21_tape4 = s21_real_tape[4] + 1j*s21_imag_tape[4]
# s21_list = [s21_tapeno, s21_tape1, s21_tape2, s21_tape3, s21_tape4]
# title = 's21_tape_fr_shift_shorted_smith'
# overlay_smith(s21_list, file_list_leg, title, output_path+title)
# print("plot done! Feedline resonance cable dependent! ")

# b2_warm_imag_long = folder_bridge+'30mhz_200mhz_shorted_notape_imag_long'+'.csv'
# b2_warm_real_long = folder_bridge+'30mhz_200mhz_shorted_notape_real_long'+'.csv'
# b2_warm_imag_tape3_long = folder_bridge+'30mhz_200mhz_shorted_tape3_imag_long'+'.csv'
# b2_warm_real_tape3_long = folder_bridge+'30mhz_200mhz_shorted_tape3_real_long'+'.csv'
# b2_warm_imag_tape4_long = folder_bridge+'30mhz_200mhz_shorted_tape4_imag_long'+'.csv'
# b2_warm_real_tape4_long = folder_bridge+'30mhz_200mhz_shorted_tape4_real_long'+'.csv'

# file_list_imag = [b2_warm_imag_long, b2_warm_imag_tape3_long, b2_warm_imag_tape4_long]
# file_list_real = [b2_warm_real_long, b2_warm_real_tape3_long, b2_warm_real_tape4_long]
# file_list_leg = ['junk', 'junk', 'junk'] 
# title = 'kid_bridge_tape_fr_change_longcable'
# scan_range = [30, 200] # mhz
# title = 'junk'

# freq_list, s21_imag_tape = overlay_traces(file_list_imag, file_list_leg, title, scan_range, False, output_path+title)
# freq_list, s21_real_tape = overlay_traces(file_list_real, file_list_leg, title, scan_range, False, output_path+title)
# s21_tapeno_long = 20*np.log10(np.sqrt(s21_real_tape[0]**2 + s21_imag_tape[0]**2))
# s21_tape3_long = 20*np.log10(np.sqrt(s21_real_tape[1]**2 + s21_imag_tape[1]**2))
# s21_tape4_long = 20*np.log10(np.sqrt(s21_real_tape[2]**2 + s21_imag_tape[2]**2))
# s21_list = [s21_tapeno_long, s21_tape3_long, s21_tape4_long]
# file_list_leg = ['b2_notape_longcable', 'b2_tape3_longcable', 'b2_tape4_longcable'] 
# title = 's21_tape_fr_shift_shorted_longcable'
# overlay_s21mag(freq_list[0], s21_list, file_list_leg, title, scan_range, False, output_path+title)
# print("plot done! Feedline resonance cable dependent! ")

# b2_warm_imag_long_ln_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_imag_long_LN'+'.csv'
# b2_warm_real_long_ln_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_real_long_LN'+'.csv'

# file_list_imag = [b2_warm_imag_tape3, b2_warm_imag_tape3_long, b2_warm_imag_long_ln_tape3]
# file_list_real = [b2_warm_real_tape3, b2_warm_real_tape3_long, b2_warm_real_long_ln_tape3]
# file_list_leg = ['junk', 'junk', 'junk'] 
# title = 'kid_bridge_tape_fr_change_longcable_ln'
# scan_range = [30, 200] # mhz
# title = 'junk'

# freq_list, s21_imag_tape = overlay_traces(file_list_imag, file_list_leg, title, scan_range, False, output_path+title)
# freq_list, s21_real_tape = overlay_traces(file_list_real, file_list_leg, title, scan_range, False, output_path+title)
# s21_tape3_short = 20*np.log10(np.sqrt(s21_real_tape[0]**2 + s21_imag_tape[0]**2))
# s21_tape3_long = 20*np.log10(np.sqrt(s21_real_tape[1]**2 + s21_imag_tape[1]**2))
# s21_tape3_long_ln = 20*np.log10(np.sqrt(s21_real_tape[2]**2 + s21_imag_tape[2]**2))
# s21_list = [s21_tape3_short, s21_tape3_long, s21_tape3_long_ln]
# file_list_leg = ['b2_tape3_shortcable', 'b2_tape3_longcable', 'b2_tape3_longcable_ln'] 
# title = 's21_tape_fr_shift_shorted_longcable_ln'
# overlay_s21mag(freq_list[0], s21_list, file_list_leg, title, scan_range, False, output_path+title)
# print("plot done! Feedline resonance cable dependent! ")

###########################################
###########################################
###########################################

# ######################################### 6/27/2025 a few more plots:)) 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-6-27-pcb-moreschecks/'
# folder_feedline = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-24-pcb-testmore/'
# folder_more_feedline = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-26-pcb-feedline/'

# b1_cold_long = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_4_cold_b1'+'.csv'
# feedline_warm_long = folder_more_feedline+'2025_6_27_feedline_room_ifbw300hz_3dbm_ave10_long_mag_10mhz_300mhz_16000pts'+'.csv'
# feedline_warm_short = folder_more_feedline+'2025_6_27_feedline_room_ifbw300hz_3dbm_ave10_long_mag_10mhz_300mhz_16000pts_shortcable'+'.csv'

# file_list = [b1_cold_long, feedline_warm_long, feedline_warm_short]
# file_list_leg = ['b1_cold_long', 'feedline_warm_long', 'feedline_warm_short']
# title = 'scan_b1_feedline'
# scan_range = [10, 290] # mhz
# scan_range_y = [-2.6, 0] # mhz

# freq_list, s21_list = overlay_traces(file_list, file_list_leg, title, scan_range, scan_range_y, output_path+title)
# print("plot done! Feedline resonance! ")

# b1_cold_real = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_10_cold_b1_real'+'.csv'
# b1_cold_imag = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_10_cold_b1_imag'+'.csv'

# title = 'junk'
# file_list_real = [b1_cold_imag]
# file_list_imag = [b1_cold_real]
# file_list_leg = ['junk'] 

# freq_list, s21_imag_b1_cold = overlay_traces(file_list_imag, file_list_leg, title, scan_range, False, output_path+title)
# freq_list, s21_real_b1_cold = overlay_traces(file_list_real, file_list_leg, title, scan_range, False, output_path+title)
# s21_b1_cold = s21_real_b1_cold[0] + 1j*s21_imag_b1_cold[0]

# peak_search(freq_list[0], s21_b1_cold, save_path, fwindow=5e-4, start_f=None, stop_f=None, nsig=3, max_N_peaks=10)

###########################################
###########################################
###########################################

# # ######################################### 6/30/2025 resistor test:)) 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-6-30-pcb-resistor/'
# folder_resistor = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-30-pcb-resistor/'
# folder_bridge = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-26-pcb-bridge/'

# b2_warm_imag_tape3_long = folder_bridge+'30mhz_200mhz_shorted_tape3_imag_long'+'.csv'
# b2_warm_real_tape3_long = folder_bridge+'30mhz_200mhz_shorted_tape3_real_long'+'.csv'
# b2_imag_long_ln_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_imag_long_LN'+'.csv'
# b2_real_long_ln_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_real_long_LN'+'.csv'
# b2_warm_imag_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_imag'+'.csv'
# b2_warm_real_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_real'+'.csv'
# b2_warm_freq, b2_warm_imag_tape3_long = load_noise_data(b2_warm_imag_tape3_long)
# b2_warm_freq, b2_warm_real_tape3_long = load_noise_data(b2_warm_real_tape3_long)
# b2_warm_tape3_long = 20*np.log10(np.sqrt(b2_warm_real_tape3_long**2 + b2_warm_imag_tape3_long**2))
# b2_warm_freq, b2_ln_imag_tape3_long = load_noise_data(b2_imag_long_ln_tape3)
# b2_warm_freq, b2_ln_real_tape3_long = load_noise_data(b2_real_long_ln_tape3)
# b2_ln_tape3_long = 20*np.log10(np.sqrt(b2_ln_real_tape3_long**2 + b2_ln_imag_tape3_long**2))
# b2_warm_freq, b2_warm_imag_tape3_short = load_noise_data(b2_warm_imag_tape3)
# b2_warm_freq, b2_warm_real_tape3_short = load_noise_data(b2_warm_real_tape3)
# b2_warm_tape3_short = 20*np.log10(np.sqrt(b2_warm_real_tape3_short**2 + b2_warm_imag_tape3_short**2))

# b2_solder_preres = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PreResolder'+'.csv'
# b2_solder_preres_freq, b2_solder_org_s21mag = load_noise_data(b2_solder_preres)

# freq_list = [b2_warm_freq, b2_warm_freq, b2_warm_freq, b2_solder_preres_freq]
# s21_list = [b2_warm_tape3_short, b2_warm_tape3_long, b2_ln_tape3_long, b2_solder_org_s21mag]
# file_list_leg = ['b2_warm_short_tape3', 'b2_warm_long_tape3', 'b2_ln_long_tape3', 'b2_warm_short_presolder_tape3'] 
# title = 'b2_bridge_change'
# scan_range = [30, 300] # mhz
# scan_range_y = [-7, -1]

# overlay_s21mag(freq_list, s21_list, file_list_leg, title, scan_range, False, output_path+title)
# print("plot done! Reproducing results??")

# b2_solder_bread = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder'+'.csv'
# b2_solder_bread_freq, b2_solder_bread_s21mag = load_noise_data(b2_solder_bread)

# b2_solder_notape = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_NoTape_0Ohm'+'.csv'
# b2_solder_notape_freq, b2_solder_notape_s21mag = load_noise_data(b2_solder_notape)

# b2_solder_org_real = folder_bridge+'30mhz_200mhz_shorted_notape_imag'+'.csv'
# b2_solder_org_freq, b2_solder_org_real = load_noise_data(b2_solder_org_real)
# b2_solder_org_imag = folder_bridge+'30mhz_200mhz_shorted_notape_real'+'.csv'
# b2_solder_org_freq, b2_solder_org_imag = load_noise_data(b2_solder_org_imag)
# b2_solder_org_mag = 20*np.log10(np.sqrt(b2_solder_org_real**2 + b2_solder_org_imag**2))

# freq_list = [b2_solder_preres_freq, b2_solder_bread_freq, b2_solder_bread_freq, b2_solder_org_freq]
# s21_list = [b2_solder_org_s21mag, b2_solder_bread_s21mag, b2_solder_notape_s21mag, b2_solder_org_mag]
# file_list_leg = ['b2_presolder_tape3', 'b2_breadboard_tape3', 'b2_breadboard_notape', 'b2_short_notape'] 
# title = 'b2_bridge_breadboard'
# scan_range = [30, 100] # mhz
# scan_range_y = [-10, 0]

# overlay_s21mag(freq_list, s21_list, file_list_leg, title, scan_range, scan_range_y, output_path+title)
# print("plot done!")

# b2_0ohm_imag = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Imag_NoTape_0Ohm'+'.csv'
# b2_0ohm_freq, b2_0ohm_imag = load_noise_data(b2_0ohm_imag)
# b2_0ohm_real = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Real_NoTape_0Ohm'+'.csv'
# b2_0ohm_freq, b2_0ohm_real = load_noise_data(b2_0ohm_real)
# b2_0ohm_mag = 20*np.log10(np.sqrt(b2_0ohm_real**2 + b2_0ohm_imag**2))
# b2_0ohm_z = b2_0ohm_real + 1j*b2_0ohm_imag

# b2_10ohm_imag = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Imag_NoTape_10ohm'+'.csv'
# b2_10ohm_freq, b2_10ohm_imag = load_noise_data(b2_10ohm_imag)
# b2_10ohm_real = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Real_NoTape_10ohm'+'.csv'
# b2_10ohm_freq, b2_10ohm_real = load_noise_data(b2_10ohm_real)
# b2_10ohm_mag = 20*np.log10(np.sqrt(b2_10ohm_real**2 + b2_10ohm_imag**2))
# b2_10ohm_z = b2_10ohm_real + 1j*b2_10ohm_imag

# b2_20ohm_imag = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Imag_NoTape_20ohm'+'.csv'
# b2_20ohm_freq, b2_20ohm_imag = load_noise_data(b2_20ohm_imag)
# b2_20ohm_real = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Real_NoTape_20ohm'+'.csv'
# b2_20ohm_freq, b2_20ohm_real = load_noise_data(b2_20ohm_real)
# b2_20ohm_mag = 20*np.log10(np.sqrt(b2_20ohm_real**2 + b2_20ohm_imag**2))
# b2_20ohm_z = b2_20ohm_real + 1j*b2_20ohm_imag

# b2_1kohm_imag = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Imag_NoTape_1kohm'+'.csv'
# b2_1kohm_freq, b2_1kohm_imag = load_noise_data(b2_1kohm_imag)
# b2_1kohm_real = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Real_NoTape_1kohm'+'.csv'
# b2_1kohm_freq, b2_1kohm_real = load_noise_data(b2_1kohm_real)
# b2_1kohm_mag = 20*np.log10(np.sqrt(b2_1kohm_real**2 + b2_1kohm_imag**2))
# b2_1kohm_z = b2_1kohm_real + 1j*b2_1kohm_imag

# b2_open_imag = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Imag_NoTape_open'+'.csv'
# b2_open_freq, b2_open_imag = load_noise_data(b2_open_imag)
# b2_open_real = folder_resistor+'2025_6_30_IFBW300Hz_30MHzto300MHz_PostResolder_Real_NoTape_open'+'.csv'
# b2_open_freq, b2_open_real = load_noise_data(b2_open_real)
# b2_open_mag = 20*np.log10(np.sqrt(b2_open_real**2 + b2_open_imag**2))
# b2_open_z = b2_open_real + 1j*b2_open_imag

# s21_list = [b2_solder_notape_s21mag, b2_10ohm_mag, b2_20ohm_mag, b2_1kohm_mag, b2_open_mag]
# freq_list = [b2_0ohm_freq, b2_10ohm_freq, b2_20ohm_freq, b2_1kohm_freq, b2_open_freq]
# file_list_leg = ['b2_0ohm_notape', 'b2_10ohm', 'b2_20ohm', 'b2_1kohm', 'b2_open'] 
# title = 'b2_resistance_change'
# scan_range = [30, 100] # mhz
# scan_range_y = [-10, 0]

# overlay_s21mag(freq_list, s21_list, file_list_leg, title, scan_range, scan_range_y, output_path+title)
# print("plot done!")

# title = 'b2_resistance_change_polar'
# s21_list = [b2_0ohm_z, b2_10ohm_z, b2_20ohm_z]
# freq_list = [b2_0ohm_freq, b2_10ohm_freq, b2_20ohm_freq]
# masked_freq_list, masked_s21_list = mask_frequency_range_list(freq_list, s21_list, 40*1e6, 100*1e6)
# file_list_leg = ['b2_0ohm_notape', 'b2_10ohm', 'b2_20ohm'] 
# overlay_smith(s21_list, file_list_leg, title, output_path+title)
# title = 'b2_resistance_change_polar_masked'
# overlay_smith(masked_s21_list, file_list_leg, title, output_path+title)
# print("plot done!")

# # title = 'b2_fit_vna'
# # overlay_fit(masked_s21_list, masked_freq_list, file_list_leg, title, output_path+title)
# # print("plot done!")

# folder_son = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-30-son/'
# son_bridge = folder_son+'2025-6-26-paa-fr4-bridge_s21'+'.csv'
# freq_mhz_son, s21_real_son, s21_imag_son = read_s21_son(son_bridge)
# s21_mag_son = 20*np.log10(np.sqrt(s21_real_son**2 + s21_imag_son**2))
# folder_feedline = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-24-pcb-testmore/'
# b1_cold_long = folder_feedline+'2025_6_24_xjw_vna_ifbw_100_hz_3dbm_ave_4_cold_b1'+'.csv'
# b1_cold_long_freq, b1_cold_long_s21 = load_noise_data(b1_cold_long)

# s21_list = [s21_mag_son, b1_cold_long_s21]
# freq_list = [freq_mhz_son*1e6, b1_cold_long_freq]
# file_list_leg = ['sonnet sim', 'b1_cold_long'] 
# title = 'no_resonance_expected'
# scan_range = [30, 200] # mhz
# # scan_range_y = [-10, 0]

# overlay_s21mag(freq_list, s21_list, file_list_leg, title, scan_range, False, output_path+title)
# print("plot done!")

# son_bridge_short = folder_son+'2025-6-26-paa-fr4-bridge_layers_s21'+'.csv'
# freq_mhz_son_short, s21_real_son_short, s21_imag_son_short = read_s21_son(son_bridge_short)
# s21_mag_son_short = 20*np.log10(np.sqrt(s21_real_son_short**2 + s21_imag_son_short**2))
# s21_list = [s21_mag_son_short, b2_solder_org_mag]
# freq_list = [freq_mhz_son_short*1e6, b2_solder_org_freq]
# file_list_leg = ['sonnet_sim', 'b2_short_notape'] 
# title = 'resonance_found_with_direct_coupling'
# scan_range = [30, 150] # mhz
# # scan_range_y = [-10, 0]

# overlay_s21mag(freq_list, s21_list, file_list_leg, title, scan_range, False, output_path+title)

# s21_list = [s21_real_son_short+1j*s21_imag_son_short, b2_solder_org_real+1j*b2_solder_org_imag]
# overlay_smith(s21_list, file_list_leg, title, output_path+title+'polar')
# print("plot done!")

# masked_freq_list, masked_s21_list = mask_frequency_range_list(freq_list, s21_list, 40*1e6, 140*1e6)
# title = 'b2_son_polar_masked'
# overlay_smith(masked_s21_list, file_list_leg, title, output_path+title)
# print("plot done!")

# title = 'b2_son_fit_vna'
# overlay_fit(masked_s21_list[:1], masked_freq_list[:1], file_list_leg, title, output_path+title)
# print("plot done!")

# b2_warm_imag = folder_bridge+'30mhz_200mhz_shorted_notape_imag'+'.csv'
# b2_warm_freq, b2_warm_imag_tape0 = load_noise_data(b2_warm_imag)
# b2_warm_real = folder_bridge+'30mhz_200mhz_shorted_notape_real'+'.csv'
# b2_warm_freq, b2_warm_real_tape0 = load_noise_data(b2_warm_real)
# b2_warm_imag_tape1 = folder_bridge+'30mhz_200mhz_shorted_tape1_imag'+'.csv'
# b2_warm_freq, b2_warm_imag_tape1 = load_noise_data(b2_warm_imag_tape1)
# b2_warm_real_tape1 = folder_bridge+'30mhz_200mhz_shorted_tape1_real'+'.csv'
# b2_warm_freq, b2_warm_real_tape1 = load_noise_data(b2_warm_real_tape1)
# b2_warm_imag_tape2 = folder_bridge+'30mhz_200mhz_shorted_tape2_imag'+'.csv'
# b2_warm_freq, b2_warm_imag_tape2 = load_noise_data(b2_warm_imag_tape2)
# b2_warm_real_tape2 = folder_bridge+'30mhz_200mhz_shorted_tape2_real'+'.csv'
# b2_warm_freq, b2_warm_real_tape2 = load_noise_data(b2_warm_real_tape2)
# b2_warm_imag_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_imag'+'.csv'
# b2_warm_freq, b2_warm_imag_tape3 = load_noise_data(b2_warm_imag_tape3)
# b2_warm_real_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_real'+'.csv'
# b2_warm_freq, b2_warm_real_tape3 = load_noise_data(b2_warm_real_tape3)
# b2_warm_imag_tape4 = folder_bridge+'30mhz_200mhz_shorted_tape4_imag'+'.csv'
# b2_warm_freq, b2_warm_imag_tape4 = load_noise_data(b2_warm_imag_tape4)
# b2_warm_real_tape4 = folder_bridge+'30mhz_200mhz_shorted_tape4_real'+'.csv'
# b2_warm_freq, b2_warm_real_tape4 = load_noise_data(b2_warm_real_tape4)

# s21_list = [b2_warm_real_tape0+1j*b2_warm_imag_tape0,
# b2_warm_real_tape1+1j*b2_warm_imag_tape1,
# b2_warm_real_tape2+1j*b2_warm_imag_tape2,
# b2_warm_real_tape3+1j*b2_warm_imag_tape3,
# b2_warm_real_tape4+1j*b2_warm_imag_tape4]
# freq_list = [b2_warm_freq, b2_warm_freq, b2_warm_freq, b2_warm_freq, b2_warm_freq]
# file_list_leg = ['b2_notape_shortcable', 'b2_tape1_shortcable', 'b2_tape2_shortcable', 'b2_tape3_shortcable', 'b2_tape4_shortcable'] 
# title = 'check_q_when_capacitance_changes'
# overlay_smith(s21_list, file_list_leg, title, output_path+title)
# print("plot done!")

# masked_freq_list, masked_s21_list = mask_frequency_range_list(freq_list, s21_list, 50*1e6, 80*1e6)
# title = 'check_q_when_capacitance_changes'
# overlay_smith(masked_s21_list, file_list_leg, title, output_path+title)
# print("plot done!")

# # title = 'b2_capacitance_fit'
# # overlay_fit(masked_s21_list, masked_freq_list, file_list_leg, title, output_path+title)
# # print("plot done!")

# title = 'toy_s21_checkqc_qi'
# plot_s21_magnitude(output_path+title)

###########################################
###########################################
###########################################

# # ######################################### 7/1/2025 more qc qi test:)) 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-7-1-pcb-lnqi/'
# folder_bridge = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-26-pcb-bridge/'
# b2_warm_imag_tape3_long = folder_bridge+'30mhz_200mhz_shorted_tape3_imag_long'+'.csv'
# b2_warm_real_tape3_long = folder_bridge+'30mhz_200mhz_shorted_tape3_real_long'+'.csv'
# b2_warm_freq, b2_warm_imag = load_noise_data(b2_warm_imag_tape3_long)
# b2_warm_freq, b2_warm_real = load_noise_data(b2_warm_real_tape3_long)

# b2_warm_imag_long_ln_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_imag_long_LN'+'.csv'
# b2_warm_real_long_ln_tape3 = folder_bridge+'30mhz_200mhz_shorted_tape3_real_long_LN'+'.csv'
# b2_ln_freq, b2_ln_imag = load_noise_data(b2_warm_imag_long_ln_tape3)
# b2_ln_freq, b2_ln_real = load_noise_data(b2_warm_real_long_ln_tape3)

# s21_list = [20*np.log10(np.sqrt(b2_warm_real**2 + b2_warm_imag**2)), 
# 20*np.log10(np.sqrt(b2_ln_real**2 + b2_ln_imag**2))]

# freq_list = [b2_warm_freq, b2_ln_freq]
# file_list_leg = ['b2_warm_long_tape3', 'b2_ln_long_tape3'] 
# title = 'ln_not_changing_qi'
# scan_range = [60, 100] # mhz
# scan_range_y = [-8, 0]

# overlay_s21mag(freq_list, s21_list, file_list_leg, title, scan_range, scan_range_y, output_path+title)
# print("plot done!")

# s21_list = [b2_warm_real+1j*b2_warm_imag, b2_ln_real+1j*b2_ln_imag]
# overlay_smith(s21_list, file_list_leg, title, output_path+title+'_polar')
# print("plot done!")

# masked_freq_list, masked_s21_list = mask_frequency_range_list(freq_list, s21_list, 60*1e6, 120*1e6)
# overlay_smith(masked_s21_list, file_list_leg, title, output_path+title+'_masked')
# print("plot done!")

# overlay_fit(masked_s21_list, masked_freq_list, file_list_leg, title, output_path+title+'_fit')
# print("plot done!")

###########################################
###########################################
###########################################

# # # ######################################### 7/2/2025 lower conductivity not quite reaching low q_i seen? 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-7-2-pcb-qilow/'
# folder_son = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-30-son/'
# folder_bridge = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-6-26-pcb-bridge/'

# son_bridge_short = folder_son+'2025-6-26-paa-fr4-bridge_layers_s21'+'.csv'
# son_bridge_short_low_cond = folder_son+'2025-6-26-paa-fr4-bridge_layers_lower_cond_s21'+'.csv'
# freq_mhz_son_short, s21_real_son_short, s21_imag_son_short = read_s21_son(son_bridge_short)
# freq_mhz_son_short_low_cond, s21_real_son_short_low_cond, s21_imag_son_short_low_cond = read_s21_son(son_bridge_short_low_cond)

# b2_solder_org_real = folder_bridge+'30mhz_200mhz_shorted_notape_imag'+'.csv'
# b2_solder_org_freq, b2_solder_org_real = load_noise_data(b2_solder_org_real)
# b2_solder_org_imag = folder_bridge+'30mhz_200mhz_shorted_notape_real'+'.csv'
# b2_solder_org_freq, b2_solder_org_imag = load_noise_data(b2_solder_org_imag)
# b2_solder_org_mag = 20*np.log10(np.sqrt(b2_solder_org_real**2 + b2_solder_org_imag**2))

# s21_mag_son_short = 20*np.log10(np.sqrt(s21_real_son_short**2 + s21_imag_son_short**2))
# s21_mag_son_short_low_cond = 20*np.log10(np.sqrt(s21_real_son_short_low_cond**2 + s21_imag_son_short_low_cond**2))
# s21_list = [s21_mag_son_short, s21_mag_son_short_low_cond, b2_solder_org_mag]
# freq_list = [freq_mhz_son_short*1e6, freq_mhz_son_short_low_cond*1e6, b2_solder_org_freq]
# file_list_leg = ['sonnet_sim', 'sonnet_sim_worse_cond', 'b2_short_notape'] 
# title = 'son_sim_lower_cond_explain_some_low_qi'
# scan_range = [30, 150] # mhz
# # scan_range_y = [-10, 0]

# overlay_s21mag(freq_list, s21_list, file_list_leg, title, scan_range, False, output_path+title)
# print("plot done!")

# s21_list = [s21_real_son_short+1j*s21_imag_son_short, s21_real_son_short_low_cond+1j*s21_imag_son_short_low_cond, b2_solder_org_real+1j*b2_solder_org_imag]
# overlay_smith(s21_list, file_list_leg, title, output_path+title+'_polar')
# print("plot done!")

# masked_freq_list, masked_s21_list = mask_frequency_range_list(freq_list, s21_list, 30*1e6, 150*1e6)
# title = 'b2_son_polar_masked_low_cond'
# overlay_smith(masked_s21_list, file_list_leg, title, output_path+title)
# print("plot done!")

# title = 'b2_son_fit_vna_low_cond'
# overlay_fit(masked_s21_list[:2], masked_freq_list[:2], file_list_leg, title, output_path+title)
# print("plot done!")


###########################################
###########################################
###########################################

# # ######################################### 7/6/2025 is there really a resonance?? 
output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-7-6-son-reswhere/'
folder_son = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/data/2025-7-6-son-data-repstrange/'

son_repstrange = folder_son+'2025-7-2-paa-lossless_idealcap_rep_s21'+'.csv'
freq_mhz, s21_real, s21_imag, s21_mag, s21_z = read_s21_son(son_repstrange)

s21_list = [s21_mag]
freq_list = [freq_mhz*1e6]
file_list_leg = ['sonnet_sim'] 
title = 'son_sim_cannot_rep'
scan_range = [0, 1000] # mhz
# scan_range_y = [-10, 0]

overlay_s21mag(freq_list, s21_list, file_list_leg, title, scan_range, False, output_path+title)
print("plot done!")

s21_list = [s21_z]
overlay_smith(s21_list, file_list_leg, title, output_path+title+'_polar')
print("plot done!")

masked_freq_list, masked_s21_list = mask_frequency_range_list(freq_list, s21_list, 60*1e6, 360*1e6)
overlay_smith(masked_s21_list, file_list_leg, title, output_path+title+'_masked')
print("plot done!")

overlay_fit(s21_list, freq_list, file_list_leg, title, output_path+title+'_fit', 220*1e-3)
print("plot done!")
