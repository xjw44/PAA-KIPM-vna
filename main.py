import sys
sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna')
from utils.plotVNA import *
from utils.findPeak import *
from utils.fitres import *
from utils.toy_S21 import *
import configData
from MB.plotMB import *

# # ######################################### 7/10/2025 find the qc, fr, mb
# overlay_s21mag(configData.son_files["freq_list"], 
# 	configData.son_files["s21_list"], 
# 	configData.son_files["file_list_leg"], 
# 	configData.son_files["title"], 
# 	configData.son_files["scan_range"], False, 
# 	configData.son_files["output_path"]+configData.son_files["title"])
# print("plot done!")

# overlay_smith(configData.son_files["masked_z_list"], 
# 	configData.son_files["file_list_leg"], 
# 	configData.son_files["title"], 
# 	configData.son_files["output_path"]+configData.son_files["title"]+'_polar_mask')
# print("plot done!")

# overlay_fit(configData.son_files["masked_z_list"], 
# 	configData.son_files["masked_freq_list"],
# 	configData.son_files["file_list_leg"], 
# 	configData.son_files["title"], 
# 	configData.son_files["output_path"]+configData.son_files["title"]+'_fit',
# 	configData.son_files["fr_est"])
# print("plot done!")

# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna/output/2025-7-10-hf-m-b/'
# title = 'nqp-exp'
# plotN_qp(output_path+title)
# print("plot done!")

# overlay_fit(s21_list, freq_list, file_list_leg, title, output_path+title+'_fit', 220*1e-3)
# print("plot done!")


