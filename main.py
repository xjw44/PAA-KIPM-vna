import sys
# sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb')
from mb.plotMB import *
from mb.mbEquations import *
from eff.plot_eff import *
from res.plotRes import *

from pathlib import Path
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

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
# title = 'amp-res-fin'
# amp_res_vs_eabs(output_path+title, debug=False)
# print("plot done!")

# title = 'res-all-fin'
# # compare_resolution(output_path+title)
# compare_resolution_sub(output_path+title)
# print("plot done!")

output_path = 'output/2025-9-10-cpadmore/'
title = 'thre-noise-eff'
plot_psd_all(output_path+title)
print("plot done!")



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
# output_path = 'output/2025-8-19-eabs/'
# title = 'eabs-exp-fin'
# plot_energy_response(output_path+title)
# print("plot done!")
