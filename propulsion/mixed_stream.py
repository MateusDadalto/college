from numpy import sqrt, NaN

MFP = lambda mach, gamma, r: (
        mach * sqrt(gamma_t / r) / (1 + (gamma - 1) / 2 * mach ** 2) ** ((gamma + 1) / (2 * (gamma - 1))))

phi = lambda mach, gamma: mach ** 2 * (1 + (gamma - 1) / 2 * mach ** 2) / (1 + gamma * mach ** 2) ** 2

gamma_c = 1.4

c_pc = 1004.8

P_0 = 19391

T_0 = 216.65

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

h_PR = 42798.4 * 1000

T_t4 = 1666.7

T_t7 = 2000

eta_b = 0.99

pi_b = 0.96

eta_AB = 0.95

pi_AB = 0.94

R_c = 287

R_t = 291


def calculate(pi_c, pi_f, M_0, has_after_burner=True):
    # pi_c = 20  # variando até 40

    # pi_f = 2  # variando até 5

    # M_0 = 0.9  # 2.0

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

    tau_cH = pi_cH ** ((gamma_c - 1) / (gamma_c * e_f))

    eta_cH = (pi_cH ** ((gamma_c - 1) / gamma_c) - 1) / (tau_cH - 1)

    '''COMBUSTOR (3-4)'''

    tau_lambda = c_pt * T_t4 / (c_pc * T_0)

    # fuel - air ratio
    f_b = (tau_lambda - tau_r * tau_d * tau_f * tau_cH) / (eta_b * h_PR / (c_pc * T_0) - tau_lambda)

    '''TURBINA DE ALTA PRESSÃO (4-4.5)'''

    tau_tH = 1 - ((tau_cH - 1) * tau_r * tau_d * tau_f) / ((1 + f_b) * eta_mH * tau_lambda)

    pi_tH = tau_tH ** (gamma_t / ((gamma_t - 1) * e_tH))

    eta_tH = (tau_tH - 1) / (pi_tH ** ((gamma_t - 1) / gamma_t) - 1)

    '''TURBINA DE BAIXA PRESSAO (4.5 - 5)'''

    tau_tL = (1 / (pi_cH * pi_b * pi_tH)) ** (((gamma_t - 1) * e_tL) / gamma_t)

    pi_tL = tau_tL ** (gamma_t / ((gamma_t - 1) * e_tL))

    '''BALANCO DE ENERGIA TURBINA DE BAIXA'''

    # alpha = (eta_mL * (1 + f_b) * (tau_lambda / tau_r) * (
    #         1 - tau_tL * tau_tH) -     (tau_cH * tau_f - 1)) / (tau_f - 1)

    alpha = ((1 + f_b) * (1 - tau_tL) / (tau_f - 1) * eta_mL * eta_mH * tau_lambda * tau_tH / (tau_r * tau_d)) - 1

    alpha_prime = alpha / (1 + f_b)

    '''MISTURADOR'''

    P_t16_P_t5_ratio = 1

    M_16 = sqrt(
        2 / (gamma_c - 1) * ((P_t16_P_t5_ratio * (1 + (gamma_t - 1) / 2 * M_5 ** 2)) ** (gamma_t / (gamma_t - 1))) ** (
                (gamma_c - 1) / gamma_c) - 1)

    c_p6 = (c_pt + alpha_prime * c_pc) / (1 + alpha_prime)

    R_6 = (R_t + alpha_prime * R_c) / (1 + alpha_prime)

    gamma_6 = c_p6 / (c_p6 - R_6)

    T_t16_T_t5_ratio = T_0 * tau_f / (T_t4 * tau_tH * tau_tL)

    tau_M = c_pt / c_p6 * ((1 + alpha_prime * (c_pc / c_pt) * T_t16_T_t5_ratio) / (1 + alpha_prime))

    PHI = ((1 + alpha_prime) / (1 / (sqrt(phi(M_5, gamma_t))) + alpha_prime * sqrt(
        (R_c * gamma_t * T_t16_T_t5_ratio) / (R_t * gamma_c * phi(M_16, gamma_c))))) ** 2 * (R_6 * gamma_t * tau_M) / (
                  R_t * gamma_6)

    M_6 = sqrt(2 * PHI / (1 - 2 * gamma_6 * PHI + sqrt(1 - 2 * (gamma_6 + 1) * PHI)))

    A_16_A_5_ratio = alpha_prime * sqrt(T_t16_T_t5_ratio) / (M_16 / M_5 * sqrt(
        gamma_c * R_t / (gamma_t * R_c) * (1 + (gamma_c - 1) / 2 * M_16 ** 2) / (1 + (gamma_t - 1) / 2 * M_5 ** 2)))

    pi_M = ((1 + alpha_prime) * sqrt(tau_M)) / (1 + A_16_A_5_ratio) * MFP(M_5, gamma_t, R_t) / MFP(M_6, gamma_6, R_6)

    '''PÓS QUEIMADOR'''

    if has_after_burner:
        tau_lambdaAB = c_pAB * T_t7 / (c_p6 * T_0)

        f_AB = (tau_lambdaAB - c_pc / c_pt * tau_lambda * tau_tL * tau_tH * tau_M) / (
                eta_AB * h_PR / (c_p6 * T_0) - tau_lambdaAB)

    else:
        tau_lambdaAB = 1
        f_AB = 0

    '''TUBEIRA (7-8-9)'''

    tau_n = 1

    P_t9_P_9_ratio = P_0_P_9_ratio * pi_r * pi_d * pi_f * pi_cH * pi_b * pi_tH * pi_tL * pi_M * pi_AB * pi_n

    M_9 = sqrt(2 / (gamma_t - 1) * (P_t9_P_9_ratio ** ((gamma_t - 1) / gamma_t) - 1))

    T_t9_T_0_ratio = c_pc / c_pt * tau_lambda * tau_tH * tau_tL * tau_M * tau_lambdaAB * tau_n

    T_9_T_0_ratio = T_t9_T_0_ratio / P_t9_P_9_ratio ** ((gamma_AB - 1) / gamma_AB)

    V_9_a_0_ratio = M_9 * sqrt(gamma_AB * R_t / (gamma_c * R_c) * T_9_T_0_ratio)

    a_0 = sqrt(gamma_c * R_c * T_0)

    '''EMPUXO ESPECÍFICO'''

    ISP = a_0 * (V_9_a_0_ratio - M_0 + R_t * T_9_T_0_ratio / (R_c * V_9_a_0_ratio) * (1 - P_0_P_9_ratio) / gamma_c)

    '''CONSUMO ESPECÍFICO'''

    f = f_b / (1 + alpha) + f_AB

    S = f / ISP

    '''EFICIENCIAS'''

    eta_t = a_0 ** 2 * (V_9_a_0_ratio ** 2 - M_0 ** 2) / (2 * f * h_PR)

    eta_p = (2 * M_0 * ISP / a_0) / (V_9_a_0_ratio ** 2 - M_0 ** 2)

    eta_o = eta_p * eta_t

    return ISP, S, alpha


if __name__ == '__main__':
    print(calculate(20, 2, 2))
