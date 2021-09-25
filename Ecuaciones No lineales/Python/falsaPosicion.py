#Metodo de la falsa posicion

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def falsaPosicion(f, a, b, tol, iterMax):
    x = sp.Symbol('x')
    f1 = sp.sympify(f)
    print(f1)
    er = []
    k = 0
    xk = a
    xk_1 = b

    if f1.subs(x, a)*f1.subs(x, b) < 0: #se verifica el teorema de Bolzano
        n = sp.N(f1.subs(x, xk_1)) * (xk_1 - xk)
        d = sp.N(f1.subs(x, xk_1)) - sp.N(f1.subs(x, xk))
        xk = xk_1
        xk_1 = xk_1 - n/d #se usa la formula de la secante

        while k < iterMax:
            k += 1
            if abs(d) > tol:
                if f1.subs(x, a)*f1.subs(x, xk_1) < 0: #se verifica el teorema de Bolzano
                    #se cumple en el primer intervalo
                    b = xk_1
                    n = sp.N(f1.subs(x, xk_1)) * (xk_1 - a)
                    d = sp.N(f1.subs(x, xk_1)) - sp.N(f1.subs(x, a))
                    xk = xk_1
                    xk_1 = xk_1 - n/d #se usa la formula de la secante
                    
                    error = abs((xk_1 - xk)/xk_1)
                    er.append(sp.N(error))
                elif f1.subs(x, b)*f1.subs(x, xk_1) < 0:
                    a = xk_1
                    n = sp.N(f1.subs(x, xk_1)) * (xk_1 - b)
                    d = sp.N(f1.subs(x, xk_1)) - sp.N(f1.subs(x, b))
                    xk = xk_1
                    xk_1 = xk_1 - n/d #se usa la formula de la secante
                    
                    error = abs((xk_1 - xk)/xk_1)
                    er.append(sp.N(error))
                else:
                    break
            if error < tol:
                break

        fig, graf = plt.subplots()
        ejeX = np.arange(1, k+1, 1)
        graf.plot(ejeX, er)
        plt.show()

        return [xk_1, k, error]
    else:
        xk = "NA"
        k = "NA"
        error = "NA"
        print("El intervalo seleccionado no cumple las condiciones del teorema de Bolzano")
        return [xk, k, error]

f = 'log(x)- exp(-x) - cos(x)'
a = 1
b = 2
tol = 10**-5
iterMax = 1000

y = falsaPosicion(f, a, b, tol, iterMax)
print(y)