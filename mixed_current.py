import numpy as np

MFP = lambda mach, gamma, r: (
        mach * np.sqrt(gamma_t / r) / (1 + (gamma - 1) / 2 * mach ** 2) ** ((gamma + 1) / (2 * (gamma - 1))))

gamma_c = 1.4

c_pc = 1004.8

P_0 = 19391

T_0 = 216.65

H_0 = 12000

gamma_t = 1.3

gamma_AB = 1.3

c_pt = 1235.1

c_pAB = 1235.1

pi_dmax = 0.98

M_5 = 0.4

eta_mH = 0.99

eta_mL = 0.99

P_0_P_9_ratio = 0.9

e_c = 0.9

e_f = 0.9

e_tH = 0.91

e_tL = 0.91

pi_n = 0.98

h_pr = 42798.4 * 1000

T_t4 = 1666.7

T_t7 = 2000

eta_b = 0.99

pi_b = 0.96

eta_AB = 0.95

pi_AB = 0.94

pi_c = 20  # variando até 40

pi_f = 2  # variando até 5

M_0 = 0.9  # 2.0

'''ENTRADA DE AR (0-1)'''

tau_r = 1 + (gamma_c - 1) / 2 * M_0 ** 2

pi_r = tau_r ** (gamma_c / (gamma_c - 1))

'''ENTRADA DE AR (1-2)'''

if 1 < M_0 < 3.5:
    eta_r = 1 - 0.075 * (M_0 - 1) ** 1.35
else:
    eta_r = 1

pi_d = eta_r * pi_dmax

tau_d = pi_d ** ((gamma_c - 1) / gamma_c)

'''FAN 2-2.5 E 2-2.13'''

tau_f = pi_f ** ((gamma_c - 1) / (gamma_c * e_f))

'''COMPRESSOR DE ALTA (2.5-3)'''

pi_cH = pi_c / pi_f

tau_cH = pi_cH ** ((gamma_c - 1) / (gamma_c * e_c))

eta_cH = (pi_cH ** ((gamma_c - 1) / gamma_c) - 1) / (tau_cH - 1)

'''COMBUSTOR (3-4)'''

tau_lambda = c_pt * T_t4 / (c_pc * T_0)

# fuel - air ratio
f_b = (tau_lambda - tau_r * tau_d * tau_f * tau_cH) / (eta_b * h_pr / (c_pc * T_0) - tau_lambda)

'''TURBINA DE ALTA PRESSÃO (4-4.5)'''

tau_tH = 1 - ((tau_cH - 1) * tau_r * tau_d * tau_f) / ((1 + f_b) * eta_mH * tau_lambda)

pi_tH = tau_tH ** (gamma_t / ((gamma_t - 1) * e_tH))

eta_tH = (tau_tH - 1) / (pi_tH ** ((gamma_t - 1) / gamma_t) - 1)

'''TURBINA DE BAIXA PRESSAO (4.5 - 5)'''

tau_tL = (1 / (pi_cH * pi_b * pi_tH)) ** (((gamma_t - 1) * e_tL) / gamma_t)

pi_tL = tau_tL ** (gamma_t / ((gamma_t - 1) * e_tL))

'''BALANÇO DE ENERGIA TURBINA DE BAIXA'''

alpha = (eta_mL * (1 + f_b) * (tau_lambda / tau_r) * (
        1 - (pi_f / (pi_cH * pi_b)) ** ((gamma_t - 1) * e_tL / gamma_t)) - (tau_cH - 1)) / (tau_f - 1)

alpha_prime = alpha / (1 + f_b)

'''MISTURADOR'''

P_t16_P_t5_ratio = 1

M_16 = np.sqrt(
    2 / (gamma_c - 1) * ((P_t16_P_t5_ratio * (1 + (gamma_t - 1) / 2 * M_5 ** 2)) ** (gamma_t / (gamma_t - 1))) ** (
            (gamma_c - 1) / gamma_c) - 1)
