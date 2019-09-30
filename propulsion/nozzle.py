from math import pi, sqrt

gamma = 1.4

D = 1.3

A_e = pi * (D / 2) ** 2

A_t = A_e / 40

p_e_P_c_ratio = 9.54E-4

P_c = 300  # bar


def cf(p0):
    first_part = gamma * sqrt(2 / (gamma - 1) * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1)) * (
            1 - p_e_P_c_ratio ** ((gamma - 1) / gamma)))

    second_part = (p_e_P_c_ratio - p0 / P_c) * A_e / A_t

    return first_part + second_part
