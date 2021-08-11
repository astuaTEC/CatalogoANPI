###############################
### METODO DE LA BISECCION ####
###                        ####
###############################

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def biseccion(f, a, b, tol, iterMax):
    x = sp.Symbol('x')
    f1 = sp.sympify(f)
    er = []
    k = 0
    if f1.subs(x, a)*f1.subs(x, b) < 0: #se verifica el teorema de Bolzano
        while k < iterMax:
            k += 1
            xk = (a+b)/2
            if f1.subs(x, a)*f1.subs(x, xk) < 0: #se verifica el teorema de Bolzano
                #se cumple en el primer intervalo
                b = xk
            else:
                a = xk
            error = abs(f1.subs(x, xk))
            er.append(sp.N(error))
            if error < tol:
                break

        fig, graf = plt.subplots()
        ejeX = np.arange(1, k+1, 1)
        graf.plot(ejeX, er)
        plt.show()

        return [xk, error]
    else:
        xk = "NA"
        k = "NA"
        error = "NA"
        print("El intervalo seleccionado no cumple las condiciones del teorema de Bolzano")
        return [xk, k, error]

f = 'x^2-3'
a = 0
b = 2
tol = 10**-5
iterMax = 1000

y = biseccion(f, a, b, tol, iterMax)
print(y)
