import sys
# sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb')
from mb.plotMB import *
from mb.mbEquations import *
from eff.plot_eff import *
from res.plotRes import *
from res.plotScaling import *
from config.eff_config import *

from pathlib import Path
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

# # ######################################### 9/11/2025 cpad!! sth about amp noise is very wrong... 
# title = 'res-all-fin'
# # compare_resolution(output_path+title)
# compare_resolution_sub(output_path+title)
# print("plot done!")

# output_path = 'output/2025-9-10-cpadmore/'
# title = 'thre-noise-eff'
# plot_threshold(output_path+title)
# print("plot done!")

# title = 'nqp-range'
# scale_nqp_plot(output_path+title)
# print("plot done!")

# title = 'taur-range'
# scale_taur_plot(output_path+title)
# print("plot done!")

# title = 'qi-range'
# scale_qi_plot(output_path+title)
# print("plot done!")

# # ######################################### 9/11/2025 cpad!! sth about amp noise is very wrong... 
# output_path = 'output/2025-9-11-cpadval/'
# title = 'res-all-fin'
# compare_resolution_sub(output_path+title)
# print("plot done!")

# title = 'nqp-range'
# scale_nqp_plot(output_path+title)
# print("plot done!")

# title = 'nkids'
# scale_nkids_plot(output_path+title)
# print("plot done!")

# title = 'vol-scale'
# scale_vol_plot(output_path+title)
# print("plot done!")

# # ######################################### 9/22/2025 cpad!! sth about amp noise is very wrong... 
# output_path = 'output/2025-9-22-cpadfin/'
# # title = 'vol-scale-ticks'
# # title = 'vol-scale-checkamp'
# title = 'vol-scale-checktls'
# # title = 'vol-scale-checkamp-vol'
# scale_vol_plot(output_path+title)
# print("plot done!")

# title = 'qc-scale'
# scale_qc_plot(output_path+title)
# print("plot done!")

# title = 'pbif-scale'
# plot_eabs_pbif(output_path+title)
# print("plot done!")

# title = 'res-all-fin'
# compare_resolution_sub(output_path+title)
# print("plot done!")

# title = 'psd-sum'
# df_psd_kid, df_psd_paa, df_psd_music = df_psd_all()
# plot_psd_sum(output_path+title, df_psd_paa)
# print("plot done!")

# # ######################################### 9/23/2025 cpad!! could try to wrap up plotting today  
# output_path = 'output/2025-9-23-cpadfin/'

# df_psd_kid, df_psd_paa, df_psd_music = df_psd_all()
# title = 'psd-sum-paa'
# plot_psd_sum(output_path+title, df_psd_paa)
# print("plot done!")

# title = 'psd-sum-kid-testmatching'
# plot_psd_sum(output_path+title, df_psd_kid)
# print("plot done!")

# title = 'res-matching-kipm'
# compare_resolution_sub(output_path+title)
# print("plot done!")

# title = 'psd-sum-music'
# plot_psd_sum(output_path+title, df_psd_music)
# print("plot done!")

# title = 'dm-threshold'
# plot_dm_thresholds(output_path+title)
# print("plot done!")

# # ######################################### 9/24/2025 cpad!! presentation starts
# output_path = 'output/2025-9-24-cpadslides/'

# title = 'qi-range'
# scale_qi_plot(output_path+title)
# print("plot done!")

# title = 'qc-scale'
# scale_qc_plot(output_path+title)
# print("plot done!")

# title = 'vol-scale'
# scale_vol_plot(output_path+title)
# print("plot done!")

# title = 'res-qi-fr'
# plot_energy_response(output_path+title)
# print("plot done!")

# # ######################################### 9/30/2025 cpad!! more plots 
# output_path = 'output/2025-9-30-cpadfin/'
# title = 'res-all-fin-bar'
# # compare_resolution_sub_bar(output_path+title)
# print("plot done!")

# title = 'res-all-fin-bar-over'
# compare_resolution_overlay(output_path+title)
# print("plot done!")

# title = 'ph-eff-rep'
# plot_ph_eff_area_clean(output_path+title)
# print("plot done!")

# # ######################################### 10/1/2025 cpad!! finishing truly:))
# output_path = 'output/2025-10-1-cpadfin/'
# title = 'amp-comp-bars'
# compare_resolution_sub_bar(output_path+title)
# print("plot done!")

# title = 'gr-comp-bars'
# compare_resolution_sub_bar(output_path+title)
# print("plot done!")

# df_psd_kid, df_psd_paa, df_psd_music = df_psd_all()
# title = 'psd-sum-kid-testmatching'
# plot_psd_sum(output_path+title, df_psd_kid)
# print("plot done!")

# # ######################################### 10/2/2025 cpad!! finishing truly:))
# output_path = 'output/2025-10-1-cpadfin/'
# title = 'amp-comp-bars'
# compare_resolution_sub_bar(output_path+title)
# print("plot done!")





# # ######################################### 9/4/2025 finally there:)) just need to double check a few things:)) 
# output_path = 'output/2025-9-4-confpower/'
# title = 'res-all-volpfeed'
# compare_resolution(output_path+title)
# print("plot done!")

# title = 'amp-res-eabs_worsefeed'
# amp_res_vs_eabs(output_path+title, debug=False)
# print("plot done!")

# # ######################################### 9/4/2025 finally there:)) just need to double check a few things:)) 
# output_path = 'output/2025-9-8-conffin/'
# # title = 'res-hf-plot'
# # title = 'res-hf-plot-lines'
# title = 'res-hf-plot-lines-comp'
# plot_energy_response(output_path+title)
# print("plot done!")

# output_path = 'output/2025-9-8-psdall/'
# title = 'psd-all-tls-vol'
# plot_psd_all(output_path+title)
# print("plot done!")

# output_path = 'output/2025-9-8-resall/'
# # title = 'amp-res-fin'
# # amp_res_vs_eabs(output_path+title, debug=False)
# # print("plot done!")