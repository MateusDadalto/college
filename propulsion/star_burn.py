from numpy import pi, tan, sin, cos, rad2deg, arcsin, arccos, arctan

n = 7
R_c = 0.495  # m
l_p = 0.280  # m
c_p = 0.006  # m
alpha = pi / 6  # rad
f = 0.010  # m
epsilon = 0.875  # rad
L = 3.956  # m

'''Q 1'''

theta = pi / n

beta = epsilon * theta + (pi / 2 - alpha)

gamma = (pi / 2 + alpha) / 2

print('Q 1\n')

print(f'theta = {rad2deg(theta):.2f}', f'beta = {rad2deg(beta):.2f}', f'gamma = {rad2deg(gamma):.2f}', sep='\n')

'''Q 2 '''

rho = 1710  # kg/m³

y1 = c_p * tan(gamma)

y = 0

first_term = (l_p + f + y) ** 2 * (1 - epsilon) * theta / 2

second_term = (f + y) ** 2 * beta / 2

third_term = l_p ** 2 / 2 * sin(epsilon * theta) * (cos(epsilon * theta) - sin(epsilon * theta) / tan(alpha))

fourth_term = sin(alpha) / 2 * (c_p - y / tan(gamma)) ** 2

fifth_term = (f + y) / 2 * (2 * l_p * sin(epsilon * theta) / sin(alpha) - (f + y) / tan(alpha))

A_p = 2 * n * (first_term + second_term + third_term + fourth_term + fifth_term)

A_total = pi * R_c ** 2

mass = (A_total - A_p) * L * rho

print('\nQ 2\n')

print(f'A_p = {A_p:.2f} m²', f'Mass: {mass:.2f} kg', sep='\n')

''' Q 3 '''

y2 = R_c - f - l_p

y3 = max([l_p * sin(epsilon * theta) / cos(alpha) - f, y2])

# y4 = sqrt(l_p ** 2 + R_c ** 2 - 2 * l_p * R_c * cos(epsilon * theta)) - f

phi = arccos((R_c ** 2 - l_p ** 2 - (f + y3) ** 2) / (2 * l_p * (f + y3)))

delta = arcsin((f + y3) / R_c * sin(phi))

nu = pi / 2 - alpha - arcsin(l_p * sin(epsilon * theta) / (f + y3))

psi = (beta - delta - nu)

first_term_sleeve = (1 - epsilon) * theta * R_c ** 2 / 2

second_term_sleeve = l_p ** 2 / 2 * sin(beta) / cos(alpha) * sin(epsilon * theta)

third_term_sleeve = (f + y3) * l_p * sin(epsilon * theta) / cos(alpha) * sin(nu) / 2

fourth_term_sleeve = (f + y3) ** 2 * (beta - phi - nu) / 2

fifth_term_sleeve = - (f + y3) / 2 * l_p * sin(delta) / sin(psi) * sin(beta - phi - nu)

sixth_term_sleeve = delta / 2 * R_c ** 2

seventh_term_sleeve = - l_p ** 2 / 2 * sin(beta - nu) / sin(psi) * sin(delta)

A_p_sleeve = 2 * n * (
        first_term_sleeve + second_term_sleeve + third_term_sleeve + fourth_term_sleeve +
        fifth_term_sleeve + sixth_term_sleeve + seventh_term_sleeve)

mass_sleeve = (A_total - A_p_sleeve) * L * rho

print('\nQ 3\n')

print(f'y_3 = {y3:.2f} m', f'A_p_sleeve = {A_p_sleeve:.2f} m²', f'Sleeve Mass: {mass_sleeve:.2f} kg', sep='\n')

''' Q 4 '''

alpha_4 = arctan(1 / ((1 - epsilon) * theta + beta))

print('\nQ 4\n')

print(f'alpha = {rad2deg(alpha_4):.2f} graus')
