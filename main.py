import sys
# sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb')
from mb.plotMB import *
from mb.mbEquations import *
from eff.plot_eff import *
from res.plotRes import *

from pathlib import Path
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

# # ######################################### 8/12/2025 really want to finish the conf note now!! 
# output_path = 'output/2025-8-12-confpsd/'
# # title = 'psd-all-gr-sohigh'
# # title = 'psd-all-amp-solow'
# # title = 'psd-all-tls-changes'
# title = 'psd-all-tls-vol'
# plot_psd_all(output_path+title)
# print("plot done!")

# title = 'kappa-exp-newmusic-strange'
# plot_kappas(output_path+title)
# print("plot done!")

# # ######################################### 8/13/2025 really want to finish the conf note now!! 
# output_path = 'output/2025-8-14-eabsandeff/'
# # title = 'amp-res-eabs'
# title = 'amp-res-eabs_largeqc'
# amp_res_vs_eabs(output_path+title, debug=True)
# print("plot done!")

# title = 'kappa-exp-newmusic-strange'
# plot_kappas(output_path+title)
# print("plot done!")

# # ######################################### 8/19/2025 last stretch! energy response should be done soon! 
output_path = 'output/2025-8-19-eabs/'
# title = 'eabs-exp-fin'
# plot_energy_response(output_path+title)
# print("plot done!")




# # ######################################### 7/14/2025 check all possible kipm response:)) 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/output/2025-7-14-mb/'
# title = 'nqp0-exp'
# plotN_qp(output_path+title)
# print("plot done!")

# title = 'r-exp'
# plot_R_const(output_path+title)
# print("plot done!")

# title = 'taur-exp'
# plot_tau_r(output_path+title)
# print("plot done!")

# title = 'kappa-exp'
# plot_kappas(output_path+title)
# print("plot done!")

# title = 'qi-exp'
# plot_dissipation(output_path+title)
# print("plot done!")

# title = 'fr-exp'
# plot_frequency(output_path+title)
# print("plot done!")

# # ######################################### 7/14/2025 check all possible kipm response:)) 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/output/2025-7-16-mb/'

# title = 's21-exp'
# plot_s21(output_path+title)
# print("plot done!")

# title = 'eabs-exp-log'
# plot_energy_response(output_path+title, plot_log=True)
# print("plot done!")

# title = 'ph-eff-vol'
# plot_ph_eff_vol(output_path+title)
# print("plot done!")

# # ######################################### 7/22/2025 check eff and res and be done!! 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/output/2025-7-22-effres/'
# title = 'ph-eff-vol'
# title = 'ph-eff-vol-area'
# plot_ph_eff_vol(output_path+title)
# print("plot done!")

# title = 'ph-eff-lifetime_log'
# plot_ph_eff_area(output_path+title, plot_log=True)
# print("plot done!")

# title = 'ph-eff-lifetime_log_smallchip'
# plot_ph_eff_area(output_path+title, plot_log=True)
# print("plot done!")

# title = 'ph-eff-vol-area_thickchip'
# plot_ph_eff_vol(output_path+title)
# print("plot done!")

# title = 'psd-all'
# plot_psd_all(output_path+title)
# print("plot done!")

# # ######################################### 7/22/2025 check eff and res and be done!! 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/output/2025-7-25-tls/'

# title = 'psd-all'
# plot_psd_all(output_path+title)
# print("plot done!")

# title = 'signal-simple'
# plot_signal_time(output_path+title)
# print("plot done!")

# title = 'res-all'
# compare_resolution(output_path+title)
# print("plot done!")

# title = 'eabs-exp'
# plot_energy_response(output_path+title)
# print("plot done!")

# # ######################################### 7/22/2025 check eff and res and be done!! 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/output/2025-7-29-resfin/'

# title = 'psd-all'
# plot_psd_all(output_path+title)
# print("plot done!")

# title = 'res-all'
# compare_resolution(output_path+title)
# print("plot done!")

# # ######################################### 7/22/2025 check eff and res and be done!! 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/output/2025-7-30-reseabs/'

# title = 'res-all'
# compare_resolution(output_path+title)
# print("plot done!")

# title = 'gr-res-eabs'
# gr_res_vs_eabs(output_path+title)
# print("plot done!")

# title = 'tls-res-eabs'
# tls_res_vs_eabs(output_path+title)
# print("plot done!")

# # ######################################### 7/22/2025 check eff and res and be done!! 
# output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/output/2025-7-31-ltdfinfin/'

# title = 'tls-res-eabs-save-theta-copy'
# tls_res_vs_eabs(output_path+title)
# print("plot done!")

# title = 'amp-res-eabs'
# amp_res_vs_eabs(output_path+title)
# print("plot done!")



