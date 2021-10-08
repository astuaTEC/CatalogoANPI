import numpy as np
from sympy import *
from scipy.optimize import fminbound

def cota_traz_cubico(func, n):
    x = Symbol('x')
    f = sympify(func)

    fnum = lambdify(x, f)
    f4 = f.diff(x, 4)

    xv = np.array([1, 1.5, 1.75, 2.15, 2.4, 3], dtype=float)
    dist = np.zeros(n-1)
    #calculo de h
    for i in range(n-1):
        dist[i] = xv[i+1] - xv[i]

    h = max(dist)

    faux = -1*f4
    faux_n = lambdify(x , faux)
    a = xv[0]
    b = xv[-1]
    x_max = fminbound(faux_n, a, b)

    print(x_max)
    fn4 = lambdify(x, f4)

    cota = (5*h**4/384)*fn4(x_max)

    print(cota)


f = 'exp(x/2)'

cota_traz_cubico(f, 6)