import sys
sys.path.append('/central/home/xjw/workdir/qkid/PAA-KIPM-vna-mb/mb')
from mbEquations import *
from resEquations import *

nqp_target = 20*1e-18 # m**-3
tau_r_target = 2*1e-3 # ms

vol_kid = 11900*1e-18 # m**3
vol_paa = 21*1e-18 # m**3
vol_music = 750*1e-18 # m**3

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

k2_paa = kappa_2(t_eff_hf, f_0_nom, delta_0_hf, N_0_hf)   # [m³]
k2_kid = kappa_2(t_eff_al, f_0_nom, delta_0_al, N_0_al)   # [m³]
k2_music = kappa_2(t_eff_music, f_0_music, delta_0_al, N_0_al)   # [m³]
k1_paa = kappa_1(t_eff_hf, f_0_nom, delta_0_hf, N_0_hf)
k1_kid = kappa_1(t_eff_al, f_0_nom, delta_0_al, N_0_al)
k1_music = kappa_1(t_eff_music, f_0_music, delta_0_al, N_0_al)

qi0_nom = 1*1e6
qc0_nom = 300*1e3
qr0_nom = calculate_Qr(qi0_nom, qc0_nom)

qi0_music = 1*1e5
qc0_music = 1e4
qr0_music = calculate_Qr(qi0_music, qc0_music)

c_kid = 10.5e-12  # Farads
c_music = 26*1e-12  # Farads
c_paa = 0.1*1e-12  # Farads

f_r_kid = 1.1e9   # Hz
l_kid = calculate_inductance(c_kid, f_r_kid) # henry
l_paa = 16*1e-9 # h
l_music = 11*1e-9 # nh

rho_nom_al = 1/27*1e-6 # ohm*mum

t_ind_paa = 100*1e-9 # nm
t_ind_kid = 30*1e-9 # nm
t_ind_music = 100*1e-9 # nm

w_ind_paa = 1*1e-6 # nm
w_ind_kid = 80*1e-6 # nm
w_ind_music = 1*1e-6 # nm

pfeed_paa = compute_p_feed(f_0_nom, l_paa, qi0_nom, N_0_hf, delta_0_hf, rho_nom_al, w_ind_paa, t_ind_paa, debug=True)
pfeed_kid = compute_p_feed(f_0_nom, l_kid, qi0_nom, N_0_al, delta_0_al, rho_nom_al, w_ind_kid, t_ind_kid, debug=True)
pfeed_music = compute_p_feed(f_0_music, l_music, qi0_music, N_0_al, delta_0_al, rho_nom_al, w_ind_music, t_ind_music, debug=True)
pfeed_paa_dBm = power_to_dbm(pfeed_paa)
pfeed_kid_dBm = power_to_dbm(pfeed_kid)
pfeed_music_dBm = power_to_dbm(pfeed_music)

tn_nom = 5 # k 

tls_beta = 2 
tls_n = 1
a_c_paa = 0.004*1e-6 # m**2 
a_c_music = 0.2*1e-6 # m**2 
t_a_si_paa = 800*1e-9 # m
t_a_si_music = 1000*1e-9 # m
v_c_paa = a_c_paa*t_a_si_paa
v_c_music = a_c_music*t_a_si_music

wres_paa = compute_W_res(l_paa, N_0_hf, delta_0_hf, rho_nom_al, w_ind_paa, t_ind_paa)
wres_music = compute_W_res(l_music, N_0_al, delta_0_al, rho_nom_al, w_ind_music, t_ind_music)
eres_paa = compute_E_field(wres_paa, c_paa, t_a_si_paa)
eres_music = compute_E_field(wres_music, c_music, t_a_si_music)

froll_paa = compute_f_rolloff(f_0_nom, qr0_nom)
froll_music = compute_f_rolloff(f_0_music, qr0_music)

j_dff_tls_music_1khz = 10**-21 # 1/hz 

label_nqp_target = rf'$n_{{qp,0}}$ = {nqp_target*1e18:.0f} $\mathrm{{\mu m^{{-3}}}}$'
label_tau_r_target = rf'$\tau_r = {tau_r_target*1e3:0.2f}\,\mathrm{{ms}}$'
label_delta_0_al = rf"$\Delta_0 = {delta_0_al*1e6:.0f}\,\mathrm{{\mu eV}}$"
label_delta_0_hf = rf"$\Delta_0 = {delta_0_hf*1e6:.0f}\,\mathrm{{\mu eV}}$"
label_vol_kid = rf"$V_{{ind}} = {vol_kid*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"
label_vol_paa = rf"$V_{{ind}} = {vol_paa*1e+18:.2e}\,\mathrm{{\mu m^{{3}}}}$"
label_alpha_kid = rf"$\alpha = {alpha_kid}$"
label_alpha_paa = rf"$\alpha = {alpha_paa}$"
label_gamma_nom = rf"$\gamma = {gamma_nom}$"
label_qi0_nom = rf"$Q_{{i,0}} = {qi0_nom}$"
label_qc0_nom = rf"$Q_{{c}} = {qc0_nom}$"
label_qr0_nom = rf"$Q_{{r,0}} = {qr0_nom:.2g}$"
label_f_0_nom = rf"$f_{{r,0}} = {f_0_nom*1e-9:.0f}\,\mathrm{{GHz}}$"
label_t_eff_al = rf'$T_{{eff}}^{{Al}}$ = {t_eff_al*1e3} mK'
label_t_eff_hf = rf'$T_{{eff}}^{{Al}}$ = {t_eff_hf*1e3} mK'
label_N_0_al = rf"$N_0 = {N_0_al*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"
label_N_0_hf = rf"$N_0 = {N_0_hf*1e-18:.2e}\,\mathrm{{eV^{{-1}}\mu m^{{-3}}}}$"
label_k1_paa = rf'$\kappa_{{1}}$ = {k1_paa*1e18:.1e} $\mathrm{{\mu m^3}}$'
label_k1_kid = rf'$\kappa_{{1}}$ = {k1_kid*1e18:.1e} $\mathrm{{\mu m^3}}$'
label_k2_paa = rf'$\kappa_{{2}}$ = {k2_paa*1e18:.1e} $\mathrm{{\mu m^3}}$'
label_k2_kid = rf'$\kappa_{{2}}$ = {k2_kid*1e18:.1e} $\mathrm{{\mu m^3}}$'
label_l_paa = rf'$L_{{tot}}$ = {l_paa*1e9:.2e} nH'
label_l_kid = rf'$L_{{tot}}$ = {l_kid*1e9:.2e} nH'
label_t_ind_kid = rf"$t_{{ind}} = {t_ind_kid*1e9:.2e}\,\mathrm{{nm}}$"
label_t_ind_paa = rf"$t_{{ind}} = {t_ind_paa*1e9:.2e}\,\mathrm{{nm}}$"
label_w_ind_kid = rf"$w_{{ind}} = {w_ind_kid*1e6:.2e}\,\mathrm{{\mu m}}$"
label_w_ind_paa = rf"$w_{{ind}} = {w_ind_paa*1e6:.2e}\,\mathrm{{\mu m}}$"
label_pfeed_kid = rf"$P_{{feed}} = {pfeed_kid_dBm:.2e}$ dBm"
label_pfeed_paa = rf"$P_{{feed}} = {pfeed_paa_dBm:.2e}$ dBm"
label_pfeed_music = rf"$P_{{feed}} = {pfeed_music_dBm:.2e}$ dBm"
label_rho_nom_al = rf"$\rho_{{n,Al}} = {rho_nom_al*1e6:.2e}\,\mathrm{{\Omega\mu m}}$ "
label_tn_nom = rf"$T_{{N}} = {tn_nom:.0f}\,\mathrm{{K}}$ "

label_tls_beta = rf"$\beta_{{\mathrm{{TLS}}}} = {tls_beta}$"
label_tls_n = rf"$n_{{\mathrm{{TLS}}}} = {tls_n}$"
label_a_c_paa = rf"$A_C^{{\mathrm{{PAA}}}} = {a_c_paa*1e12:.2f}\,\mathrm{{\mu m^2}}$"
label_a_c_music = rf"$A_C^{{\mathrm{{MUSIC}}}} = {a_c_music*1e12:.2f}\,\mathrm{{\mu m^2}}$"
label_t_a_si_paa = rf"$t_{{\mathrm{{a-Si}}}}^{{\mathrm{{PAA}}}} = {t_a_si_paa*1e9:.0f}\,\mathrm{{nm}}$"
label_t_a_si_music = rf"$t_{{\mathrm{{a-Si}}}}^{{\mathrm{{MUSIC}}}} = {t_a_si_music*1e9:.0f}\,\mathrm{{nm}}$"
label_wres_paa = rf"$W_{{\mathrm{{res}}}}^{{\mathrm{{PAA}}}} = {wres_paa:.2e}\,\mathrm{{J}}$"
label_wres_music = rf"$W_{{\mathrm{{res}}}}^{{\mathrm{{MUSIC}}}} = {wres_music:.2e}\,\mathrm{{J}}$"
label_eres_paa = rf"$E_{{\mathrm{{res}}}}^{{\mathrm{{PAA}}}} = {eres_paa:.2e}\,\mathrm{{V/m}}$"
label_eres_music = rf"$E_{{\mathrm{{res}}}}^{{\mathrm{{MUSIC}}}} = {eres_music:.2e}\,\mathrm{{V/m}}$"
label_froll_paa = rf"$f_{{\mathrm{{roll}}}}^{{\mathrm{{PAA}}}} = {froll_paa:.2e}\,\mathrm{{Hz}}$"
label_froll_music = rf"$f_{{\mathrm{{roll}}}}^{{\mathrm{{MUSIC}}}} = {froll_music:.2e}\,\mathrm{{Hz}}$"

label_l_music = rf"$L^{{\mathrm{{MUSIC}}}} = {l_music*1e9:.2f}\,\mathrm{{nH}}$"
label_t_ind_music = rf"$t_{{\mathrm{{ind}}}}^{{\mathrm{{MUSIC}}}} = {t_ind_music*1e9:.0f}\,\mathrm{{nm}}$"
label_w_ind_music = rf"$w_{{\mathrm{{ind}}}}^{{\mathrm{{MUSIC}}}} = {w_ind_music*1e6:.1f}\,\mathrm{{\mu m}}$"
label_c_music = rf"$C_{{\mathrm{{MUSIC}}}} = {c_music*1e12:.0f}\,\mathrm{{pF}}$"
label_c_paa = rf"$C_{{\mathrm{{PAA}}}} = {c_paa*1e12:.1f}\,\mathrm{{pF}}$"
label_t_eff_music = rf"$T_{{\mathrm{{eff,MUSIC}}}} = {t_eff_music*1e3:.0f}\,\mathrm{{mK}}$"
label_j_dff_tls_music_1khz = rf"$J^{{\delta f_r/f_{{r,0}}}}_{{\mathrm{{TLS,MUSIC}}}}(1\,\mathrm{{kHz}}) = {j_dff_tls_music_1khz:.0e}\,\mathrm{{Hz^{{-1}}}}$"

tc_hf = delta_to_tc(delta_0_hf)
tc_al = delta_to_tc(delta_0_al)

z1_al = 1.43 # renormalization factor 
# al_b = 0.317*1e-3 # mev**-2 
al_b = 317 # ev**-2 

target_real = 1 - qr0_nom / qc0_nom

