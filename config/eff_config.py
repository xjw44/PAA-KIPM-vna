vol_kid = {
# "dmle2_single": {
# 	"v_al": 8*1e5*1e-18, # m**3
# 	"v_al_act": 2.53*1e4*1e-18,
# 	"t_al": 30*1e-9, # nm
# 	"eff_measured": 1.1*1e-2,
# 	"a_tot": (22*1e-3)**2, # mm
# },
# "OW200127": {
# 	"v_al": 21*1e3*1e-18,
# 	"v_al_act": 21*1e3*1e-18,
# 	"t_al": 30*1e-9, # nm
# 	# "v_nb": 1.23*1e6*1e-18,
# 	"v_nb": 8.6*1e6*1e-18,
# 	"eff_measured": 0.78*1e-2,
# 	"a_tot": (22*1e-3)**2, # mm
# },
"B240103": {
	"v_al": 2.19*1e4*1e-18,
	"v_al_act": 1.19*1e4*1e-18,
	"t_al": 30*1e-9, # nm
	# "v_nb": 1.8*1e5*1e-18,
	"v_nb": 0.25*1e6*1e-18,
	"eff_measured": 1.1*1e-2,
	"a_tot": (22*1e-3)**2, # mm
},
"B240103_smallchip": {
	"v_al": 2.19*1e4*1e-18,
	"v_al_act": 1.19*1e4*1e-18,
	"t_al": 30*1e-9, # nm
	# "v_nb": 1.8*1e5*1e-18,
	"v_nb": 0.25*1e6*1e-18*0.5,
	"a_tot": (22*1e-3)**2, # mm
},
"B240103_thickerchip": {
	"v_al": (1.19*20+1)*1e4*1e-18,
	"v_al_act": 1.19*1e4*1e-18*20,
	"t_al": 30*1e-9, # nm
	# "v_nb": 1.8*1e5*1e-18,
	"v_nb": 0.25*1e6*1e-18,
	"a_tot": (22*1e-3)**2, # mm
},
# "paa": {
# 	"v_al": 101*100*100*1e-12*600*1e-9,
# 	"v_al_act": 101*100*100*1e-12*600*1e-9,
# 	"v_nb": 22*1e-3*220*1e-6*70*1e-9,
# 	"t_al": 600*1e-9, # nm
# 	"a_tot": (22*1e-3)**2, # mm
# },
}

area_ph = {
"HVeV-NFC": {
	"a_frac": 0.5,
	"eff_measured": 0.95,
},
}

area_kid = {
# "dmle2_single": {
# 	"v_al": 8*1e5*1e-18, # m**3
# 	"v_al_act": 2.53*1e4*1e-18,
# 	"t_al": 30*1e-9, # nm
# 	"eff_measured": 1.1*1e-2,
# 	"a_tot": (22*1e-3)**2, # mm
# },
# "OW200127": {
# 	"v_al": 21*1e3*1e-18,
# 	"v_al_act": 21*1e3*1e-18,
# 	"t_al": 30*1e-9, # nm
# 	# "v_nb": 1.23*1e6*1e-18,
# 	"v_nb": 8.6*1e6*1e-18,
# 	"eff_measured": 0.78*1e-2,
# 	"a_tot": (22*1e-3)**2, # mm
# },
"B240103": {
	"v_al": 2.19*1e4*1e-18,
	"v_al_act": 1.19*1e4*1e-18,
	"t_al": 30*1e-9, # nm
	# "v_nb": 1.8*1e5*1e-18,
	"v_nb": 0.25*1e6*1e-18,
	"eff_measured": 1.1*1e-2,
	"a_tot": (22*1e-3)**2*2, # mm
},
"paa": {
	"v_al": 101*100*100*1e-12*600*1e-9*10,
	"v_al_act": 101*100*100*1e-12*600*1e-9*10,
	"v_nb": 22*1e-3*220*1e-6*70*1e-9,
	"t_al": 600*1e-9, # nm
	"a_tot": (10*1e-3)**2*2, # mm
},
}

# xi_nb_al_ratio=0.1292
# r_e_loss_xi_al=118000*1e-18
xi_nb_al_ratio=0.09531
r_e_loss_xi_al=413000*1e-18
eta_pb = 0.46

n_bar_abs = 0.91 # si to al 
c_s_ph = 5880 # m/s 

eta_sub_nom = 1*1e-3 # mm
eta_sub_scdms = 4*1e-3 # mm

tau_life_ph_nom = 400*1e-6 # ms

debug = True
if debug:
	for name, vals in area_kid.items():
	    v_al_act = vals.get("v_al_act")
	    t_al = vals.get("t_al")
	    if v_al_act and t_al:
	        print(f"{name}: v_al_act/t_al = {v_al_act/t_al:.3e} m²")

