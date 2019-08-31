import numpy as np
from math import sin, cos, pi, sqrt


def retornaMatrizRotacao(i, omega, capomega):
    omegarad = omega * pi / 180
    irad = i * pi / 180
    capomegarad = capomega * pi / 180

    rzomega = np.array([[cos(-omegarad), sin(-omegarad), 0],
                    [-sin(-omegarad), cos(-omegarad), 0],
                    [0, 0, 1]])
    rxi = np.array([[1, 0, 0],
                        [0, cos(-irad), sin(-irad)],
                        [0, -sin(-irad), cos(-irad)]])
    rzcapomega = np.array([[cos(-capomegarad), sin(-capomegarad), 0],
                           [-sin(-capomegarad), cos(-capomegarad), 0],
                           [0, 0, 1]])

    r1 = rzcapomega.dot(rxi)
    R = r1.dot(rzomega)
    return R

#faz a mesma coisa que o método acima, foi criado só para teste.
def retornaMatrizRotacao2(i, omega, capomega):
    omegarad = omega * pi / 180
    irad = i * pi / 180
    capomegarad = capomega * pi / 180

    rzomega = np.array([[cos(omegarad), sin(omegarad), 0],
                    [-sin(omegarad), cos(omegarad), 0],
                    [0, 0, 1]])
    rxi = np.array([[1, 0, 0],
                        [0, cos(irad), sin(irad)],
                        [0, -sin(irad), cos(irad)]])
    rzcapomega = np.array([[cos(capomegarad), sin(capomegarad), 0],
                           [-sin(capomegarad), cos(capomegarad), 0],
                           [0, 0, 1]])

    r1 = rzomega.dot(rxi)
    R = r1.dot(rzcapomega)
    return R.transpose()


def calculaU(M,e):

    u = M ##valor inicial de u
    f1 = u - e*sin(u) - M
    f2 = 1 - e*cos(u)

    while abs(f1) > 0.0001:
        u = u - f1/f2
        f1 = u - e * sin(u) - M
        f2 = 1 - e * cos(u)

    return u


mi = 398600
rt = 6378
a = 6959
i = 45
omega = 60
capomega = 30
e = 0.01617
deltat = 36835

## CALCULAR U

n = sqrt(mi/a**3)
M = deltat*n
u = 29.54*pi/180
print("{0:.30f}".format(n))

## CALCULAR n E r

r = a*(1 - e*cos(u))

## CALCULAR x, y, xdot e ydot

x = a*(cos(u)-e)
y = a*sin(u)*pow(1-pow(e, 2), 0.5)

xdot = -(n*a**2)/r * sin(u)
ydot = (n*a**2)/r * cos(u)*pow(1-pow(e, 2), 0.5)

print(x, y, xdot, ydot, '\n', sep='\n')
## MONTAR VETORES x E xdot

vecx = np.array([x, y, 0])
vecxdot = np.array([xdot, ydot, 0])

## CALCULAR MATRIZ DE ROTAÇAO

R = retornaMatrizRotacao(i, omega, capomega)
R2 = retornaMatrizRotacao2(i, omega, capomega)
print("R = {}\n\nR2 = {}".format(R, R2))
det = np.linalg.det(R)
## ROTACIONAR x E xdot
X = R.dot(vecx)
Xdot = R.dot(vecxdot)

mod1 = np.linalg.norm(Xdot)
mod2 = np.linalg.norm(vecxdot)
print(mod1,mod2)
print('X = {}'.format(X), "Xdot = {}".format(Xdot), sep='\n')

