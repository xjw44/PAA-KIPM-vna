import sys
sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb')
from mb.plotMB import *
from mb.mbEquations import *

# # ######################################### 7/14/2025 check all possible kipm response:)) 
output_path = '/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/output/2025-7-14-mb/'
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

title = 's21-exp'
plot_s21(output_path+title)
print("plot done!")


