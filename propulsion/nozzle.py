from math import pi, sqrt

gamma = 1.4

D = 1.

A_e = pi * (D / 2) ** 2

A_t = A_e / 40
p_atm = 1602


P_c = 759176  # bar

p_e_P_c_ratio = p_atm*0.282/P_c
def cf(p0):
    first_part = gamma * sqrt(2 / (gamma - 1) * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1)) * (
            1 - p_e_P_c_ratio ** ((gamma - 1) / gamma)))

    second_part = (p_e_P_c_ratio - p0 / P_c) * A_e / A_t

    return first_part + second_part
