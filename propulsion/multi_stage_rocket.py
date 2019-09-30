import numpy as np

eps1 = 5.66
eps2 = 5.25
eps3 = 4.24
eps4 = 5.82

M_f1 = 21830
M_f2 = 8490
M_f3 = 2440
M_f4 = 270

M_i1 = 50730
M_i2 = 15630
M_i3 = 6810
M_i4 = 1090

delta_v1 = 1489.15
delta_v2 = 1137.6
delta_v3 = 2315.54
delta_v4 = 2757

c1 = 1766
c2 = 1864
c3 = 2256
c4 = 1976

A = M_f1/M_i1
B = M_f2/M_i2
C = M_f3/M_i3
D = M_f4/M_i4

Mu = 100

linha1 = np.array([1, -eps1, -eps2, -eps3, -eps4])

linha2 = np.array([A, -1, -eps2, -eps3, -eps4])

linha3 = np.array([0, 0, -(B * eps2 - 1) / (B - 1), -eps3, -eps4])

linha4 = np.array([0, 0, 0, -(C * eps3 - 1) / (C - 1), -eps4])

linha5 = np.array([0, 0, 0, 0, -(D * eps4 - 1) / (D - 1)])

matriz = np.array([linha1, linha2, linha3, linha4, linha5])

print(matriz)

mu_array = [Mu, Mu, Mu, Mu, Mu]

print(np.linalg.solve(matriz, mu_array))
