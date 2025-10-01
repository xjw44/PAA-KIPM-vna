import sys
# sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/res')
# sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/mb')
from pathlib import Path

# Path to the directory where this script lives
here = Path(__file__).resolve().parent

# Append ../res and ../mb relative to config/
sys.path.append(str(here.parent / "res"))
sys.path.append(str(here.parent / "mb"))

from mbEquations import *
from resEquations import *

debug = False

total_eff_exp = 0.6
total_eff_obs_kid = 0.01

nqp_target = 20*1e18 # m**-3
nqp_music = 1000*1e18 # m**-3
tau_r_target = 2*1e-3 # ms
tau_r_music = 50*1e-6 # ms

# rho_nom_al = 1/27*1e-6 # ohm*mum
# rho_nom_hf = 33*1e-8 # ohm*mum
rho_nom_al = 0.2*1e-8 # ohm*mum
rho_nom_hf = 70*1e-8 # ohm*mum

delta_0_hf = 38*1e-6 # ev
delta_0_al = 180*1e-6 # ev

alpha_kid = 0.038
gamma_nom = -1 # 
alpha_paa = 0.46
alpha_music = 0.26
alpha_gamma_kid = alpha_kid*abs(gamma_nom)
alpha_gamma_paa = alpha_paa*abs(gamma_nom)

N_0_al = 1.72E28 # 1/(m^3*eV), Single-spin density of states (aluminum, from Jiansong's Thesis)
N_0_hf = 3.6*1e28 # 1/(m^3*eV), Single-spin density of states (aluminum, from Jiansong's Thesis)

t_eff_al = 170*1e-3 # k 
t_eff_hf = 38*1e-3 # k 
t_eff_music = 250*1e-3 # k 

f_0_nom = 4.0*1e9 # hz 
f_0_music = 300*1e6 # hz 
f_0_kid = 4e9   # Hz

k2_paa = kappa_2(t_eff_hf, f_0_nom, delta_0_hf, N_0_hf)   # [m³]
k2_kid = kappa_2(t_eff_al, f_0_nom, delta_0_al, N_0_al)   # [m³]
k2_music = kappa_2(t_eff_music, f_0_music, delta_0_al, N_0_al)   # [m³]
k1_paa = kappa_1(t_eff_hf, f_0_nom, delta_0_hf, N_0_hf)
k1_kid = kappa_1(t_eff_al, f_0_nom, delta_0_al, N_0_al)
k1_music = kappa_1(t_eff_music, f_0_music, delta_0_al, N_0_al)

qi0_nom = 1*1e6
qc0_nom = 100*1e3
qr0_nom = calculate_Qr(qi0_nom, qc0_nom)

qi0_music = 1*1e5
qc0_music = 1e4
qr0_music = calculate_Qr(qi0_music, qc0_music)

c_kid = 10.5e-12  # Farads
c_music = 21*1e-12  # Farads
c_paa = 0.1*1e-12  # Farads

indhenry_kid = calculate_inductance(c_kid, f_0_kid, debug=False) # henry
# l_paa = calculate_inductance(c_paa, f_0_kid, debug=True) # henry
indhenry_paa = 16*1e-9 # h
indhenry_music = 11*1e-9 # nh
fr_check_music = calculate_resonant_frequency(c_music, indhenry_music)
fr_check_paa = calculate_resonant_frequency(c_paa, indhenry_paa)

#######################################################
#######################################################
#######################################################

vol_kid = 11900*1e-18 # m**3
vol_paa = 21*1e-18 # m**3
vol_music = 750*1e-18 # m**3

t_ind_paa = 100*1e-9 # nm
t_ind_kid = 30*1e-9 # nm
t_ind_music = 100*1e-9 # nm

l_abs = 100*1e-6 # mum
l_sub = 1*1e-2 # cm
l_abs_num = 103

w_ind_paa = 1*1e-6 # nm
w_ind_kid = 80*1e-6 # nm
w_ind_music = 1*1e-6 # nm

l_ind_paa = vol_paa/w_ind_paa/t_ind_paa
l_ind_kid = vol_kid/w_ind_kid/t_ind_kid
l_ind_music = vol_music/w_ind_music/t_ind_music

indhenry_kid_pre = total_inductance_L(rho_nom_al, delta_0_al, alpha_kid, l_ind_kid, w_ind_kid, t_ind_kid, Delta_in_eV=True, debug=False)
indhenry_paa_pre = total_inductance_L(rho_nom_hf, delta_0_hf, alpha_paa, l_ind_paa, w_ind_paa, t_ind_paa, Delta_in_eV=True, debug=False)
indhenry_music_pre = total_inductance_L(rho_nom_al, delta_0_al, alpha_music, l_ind_music, w_ind_music, t_ind_music, Delta_in_eV=True, debug=False)

indhenry_sq = 16*1e-12 #ph/sq
indhenry_geo_paa = indhenry_paa*(1-alpha_paa)

# pfeed_paa = compute_p_feed(f_0_nom, l_paa, qi0_nom, N_0_hf, delta_0_hf, rho_nom_al, w_ind_paa, t_ind_paa)
# pfeed_kid = compute_p_feed(f_0_nom, l_kid, qi0_nom, N_0_al, delta_0_al, rho_nom_al, w_ind_kid, t_ind_kid)
# pfeed_music = compute_p_feed(f_0_music, l_music, qi0_music, N_0_al, delta_0_al, rho_nom_al, w_ind_music, t_ind_music)
# pfeed_paa_dBm = power_to_dbm(pfeed_paa)
# pfeed_kid_dBm = power_to_dbm(pfeed_kid)
# pfeed_music_dBm = power_to_dbm(pfeed_music)

pfeed_paa_vol = P_bif(N_0_hf, delta_0_hf, vol_paa, f_0_nom, qc0_nom, alpha_paa, qr0_nom, debug=False)
pfeed_paa_vol_dBm = power_to_dbm(pfeed_paa_vol, debug=debug)
pfeed_kid_vol = P_bif(N_0_al, delta_0_al, vol_kid, f_0_nom, qc0_nom, alpha_kid, qr0_nom, debug=debug)
pfeed_kid_vol_dBm = power_to_dbm(pfeed_kid_vol, debug=debug)
pfeed_kid_obs_dbm = -80 
pfeed_kid_obs = dbm_to_power(pfeed_kid_obs_dbm)
pfeed_music_vol = P_bif(N_0_al, delta_0_al, vol_music, f_0_music, qc0_music, alpha_music, qr0_music, debug=debug)
pfeed_music_vol_dBm = power_to_dbm(pfeed_music_vol, debug=debug)

nqp_r = nqp_ratio(pfeed_paa_vol, vol_paa, pfeed_kid_vol, vol_kid, debug=debug)
fano_err = fano_sigma(0.5, delta_0_hf, F=0.2, debug=debug)

tn_nom_hemt = 5 # k 
tn_nom = 0.3 # k 

tls_beta = 2 
tls_n = 0.5
a_c_paa = 0.004*1e-6 # m**2 
a_c_music = 0.6*1e-6 # m**2 
t_a_si_paa = 1000*1e-9 # m
t_a_si_music = 800*1e-9 # m
eps_r = 11.7 # silicon 
c_check_music = calculate_capacitance(eps_r, a_c_music, t_a_si_music)
c_check_paa = calculate_capacitance(eps_r, a_c_paa, t_a_si_paa)
v_c_paa = a_c_paa*t_a_si_paa
v_c_music = a_c_music*t_a_si_music

# wres_paa = compute_W_res(l_paa, N_0_hf, delta_0_hf, rho_nom_al, w_ind_paa, t_ind_paa)
wres_paa_vol = W_er(qr0_nom, qc0_nom, f_0_nom, pfeed_paa_vol)

# wres_music = compute_W_res(l_music, N_0_al, delta_0_al, rho_nom_al, w_ind_music, t_ind_music)
wres_music_vol = W_er(qr0_music, qc0_music, f_0_music, pfeed_music_vol)

# eres_paa_err = compute_E_field(wres_paa, c_paa, t_a_si_paa)
# eres_music_err = compute_E_field(wres_music, c_music, t_a_si_music)
eres_paa = compute_E_field(wres_paa_vol, c_paa, t_a_si_paa)
eres_music = compute_E_field(wres_music_vol, c_music, t_a_si_music)

froll_paa = compute_f_rolloff(f_0_nom, qr0_nom)
froll_music = compute_f_rolloff(f_0_music, qr0_music)

j_dff_tls_music_1khz = 10**-21 # 1/hz 
j_dff_tls_kipm_1khz = 10**-19 # 1/hz 
j_dff_tls_paa_1khz = update_tls_psd(t_eff_hf, t_eff_music, v_c_paa, v_c_music, 
        eres_paa, eres_music, j_dff_tls_music_1khz, tls_beta)
deltaf_paa = compute_delta_f(tau_r_target)

#######################################################
#######################################################
#######################################################

gr_res_paa = 2*1e-3 # mev 

f_kid = np.linspace(4*1e9*(1-0.001), 4*1e9*(1+0.001), 10000)  # hz
f_paa = np.linspace(4*1e9-0.005*1e9, 4*1e9+0.005*1e9, 10000)  # hz

gr_paa = compute_gr_resolution(nqp_target, vol_paa, delta_0_hf)
gr_kid = compute_gr_resolution(nqp_target, vol_kid, delta_0_al)
gr_paa_diss, gr_paa_freq = convert_eabs_res_to_s21_res(gr_paa, gr_paa, vol_paa, delta_0_hf, 
        alpha_paa, gamma_nom, k1_paa, k2_paa, qc0_nom, qr0_nom, debug=False)

# amp_ds21_res_paa = compute_amp_resolution(tau_r_target, tn_nom, pfeed_paa, debug=False)
# amp_ds21_res_kid = compute_amp_resolution(tau_r_target, tn_nom, pfeed_kid, debug=False)
amp_ds21_res_paa_vol = compute_amp_resolution(tau_r_target, tn_nom, pfeed_paa_vol, debug=debug)

# tn_ow = 4 # k 
# tau_r_ow = 100*1e-6 # mus
# pfeed_ow = dbm_to_power(-80)
# amp_ds21_res_val_ow = compute_amp_resolution(tau_r_ow, tn_ow, pfeed_ow, debug=True)
amp_ds21_res_kid_vol = compute_amp_resolution(tau_r_target, tn_nom, pfeed_kid_vol, debug=debug)
amp_ds21_res_kid_obs = compute_amp_resolution(tau_r_target, tn_nom_hemt, pfeed_kid_obs, debug=debug)

# amp_eabs_res_paa_diss, amp_eabs_res_paa_freq = convert_amp_res_to_eabs_res(amp_ds21_res_paa,
#                 vol_paa, delta_0_hf, alpha_paa, gamma_nom, k1_paa, k2_paa, qc0_nom, qr0_nom, debug=False)
# amp_eabs_res_kid_diss, amp_eabs_res_kid_freq = convert_amp_res_to_eabs_res(amp_ds21_res_kid,
#                 vol_kid, delta_0_al, alpha_kid, gamma_nom, k1_kid, k2_kid, qc0_nom, qr0_nom, debug=False)

amp_eabs_res_paa_diss_vol, amp_eabs_res_paa_freq_vol = convert_amp_res_to_eabs_res(amp_ds21_res_paa_vol,
                vol_paa, delta_0_hf, alpha_paa, gamma_nom, k1_paa, k2_paa, qc0_nom, qr0_nom, debug=debug)
amp_eabs_res_kid_diss_vol, amp_eabs_res_kid_freq_vol = convert_amp_res_to_eabs_res(amp_ds21_res_kid_vol,
                vol_kid, delta_0_al, alpha_kid, gamma_nom, k1_kid, k2_kid, qc0_nom, qr0_nom, debug=debug)
amp_eabs_res_kid_diss_obs, amp_eabs_res_kid_freq_obs = convert_amp_res_to_eabs_res(amp_ds21_res_kid_obs,
                vol_kid, delta_0_al, alpha_kid, gamma_nom, k1_kid, k2_kid, qc0_nom, qr0_nom, debug=debug)
# amp_eabs_res_ow_diss_vol, amp_eabs_res_ow_freq_vol = convert_amp_res_to_eabs_res(amp_ds21_res_val_ow,
#                 vol_kid, delta_0_al, alpha_kid, gamma_nom, k1_kid, k2_kid, qc0_nom, qr0_nom, debug=True)

#######################################################
#######################################################
#######################################################

tls_dff_paa = tls_variance(tau_r_target, j_dff_tls_paa_1khz, froll_paa, deltaf_paa)
tls_eabs_paa = convert_tls_res_to_eabs_res(tls_dff_paa, vol_paa, delta_0_hf, 
        alpha_paa, gamma_nom, k2_paa)
tls_s21_diss_paa, tls_s21_freq_paa = convert_eabs_res_to_s21_res(tls_eabs_paa, tls_eabs_paa, vol_paa, delta_0_hf, 
        alpha_paa, gamma_nom, k1_paa, k2_paa, qc0_nom, qr0_nom, debug=False)

tls_dff_kid = tls_variance(tau_r_target, j_dff_tls_music_1khz, froll_paa, deltaf_paa)
tls_eabs_kid = convert_tls_res_to_eabs_res(tls_dff_kid, vol_kid, delta_0_al, 
        alpha_kid, gamma_nom, k2_kid)

tot_freq_paa, tot_diss_paa = compute_total_resolution(gr_paa, 
        amp_eabs_res_paa_freq_vol, amp_eabs_res_paa_diss_vol, tls_eabs_paa)
tot_freq_paa_s21 = compute_total_resolution_list([gr_paa_freq, 
        amp_ds21_res_paa_vol, tls_s21_freq_paa])
tot_diss_paa_s21 = compute_total_resolution_list([gr_paa_diss, 
        amp_ds21_res_paa_vol])

tot_freq_kid, tot_diss_kid = compute_total_resolution(gr_kid, 
        amp_eabs_res_kid_freq_vol, amp_eabs_res_kid_diss_vol, tls_eabs_kid)
tot_freq_kid_obs, tot_diss_kid_obs = compute_total_resolution(gr_kid, 
        amp_eabs_res_kid_freq_obs, amp_eabs_res_kid_diss_obs, tls_eabs_kid)


#######################################################
#######################################################
#######################################################

tc_hf = delta_to_tc(delta_0_hf)
tc_al = delta_to_tc(delta_0_al)

z1_al = 1.43 # renormalization factor 
# al_b = 0.317*1e-3 # mev**-2 
al_b = 317 # ev**-2 

target_real = 1 - qr0_nom / qc0_nom

if debug: 
    print(f"[DEBUG] qi0_nom = {qi0_nom}")
    print(f"[DEBUG] qr0_nom = {qr0_nom}")
    print(f"[DEBUG] l_kid = {l_kid}")
    print(f"[DEBUG] k2_paa   = {k2_paa}")
    print(f"[DEBUG] k2_kid   = {k2_kid}")
    print(f"[DEBUG] k2_music = {k2_music}")
    print(f"[DEBUG] k1_paa   = {k1_paa}")
    print(f"[DEBUG] k1_kid   = {k1_kid}")
    print(f"[DEBUG] k1_music = {k1_music}")
    print(f"[DEBUG] j_dff_tls_paa_1khz = {j_dff_tls_paa_1khz}")
    print(f"[DEBUG] v_c_paa   = {v_c_paa}")
    print(f"[DEBUG] v_c_music = {v_c_music}")
    print(f"[DEBUG] fr_check_paa   = {fr_check_paa}")
    print(f"[DEBUG] fr_check_music = {fr_check_music}")
    print(f"[DEBUG] c_check_paa   = {c_check_paa}")
    print(f"[DEBUG] c_check_music = {c_check_music}")
    print(f"[DEBUG] wres_paa   = {wres_paa}")
    print(f"[DEBUG] wres_music = {wres_music}")
    print(f"[DEBUG] wres_paa_vol   = {wres_paa_vol}")
    print(f"[DEBUG] wres_music_vol = {wres_music_vol}")
    print(f"[DEBUG] eres_paa   = {eres_paa}")
    print(f"[DEBUG] eres_music = {eres_music}")
    print(f"[DEBUG] eres_paa_err   = {eres_paa_err}")
    print(f"[DEBUG] eres_music_err = {eres_music_err}")
    print(f"[DEBUG] froll_paa   = {froll_paa}")
    print(f"[DEBUG] froll_music = {froll_music}")
    for name, val, val_dbm in [
        ("pfeed_paa",   pfeed_paa,   pfeed_paa_dBm),
        ("pfeed_kid",   pfeed_kid,   pfeed_kid_dBm),
        ("pfeed_music", pfeed_music, pfeed_music_dBm),
    ]:
        print(f"[DEBUG] {name:<12} = {val:.6e} W  ({val_dbm:.2f} dBm)")

