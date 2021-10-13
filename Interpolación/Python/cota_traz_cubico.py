import numpy as np
from sympy import *
from scipy.optimize import fminbound
from sympy.solvers.diophantine.diophantine import length

"""
Función para calcular la cota
de error del trazador cúbico
Entradas:
    func: Función a calcularle la cota del trazador
    n: numero de puntos a utilizar en el cálculo
    xv: vector utilizado en el trazador cúbico
Salidas:
    cota: El valor numérico de la cota obtenida
"""
def cota_traz_cubico(func, n, xv):
    x = Symbol('x')
    f = sympify(func)

    #se deriva 4 veces
    f4 = f.diff(x, 4)

    dist = np.zeros(n-1)
    
    #calculo de h
    for i in range(n-1):
        dist[i] = xv[i+1] - xv[i]

    h = max(dist) #es la máxima distancia

    # Para sacar un máximo de la función en
    # un intervalo
    # es equivalente a sacarle el mínimo
    # a la función negada en ese intervalo
    # max(f) = min(-1*f)
    faux = -1*f4
    faux_n = lambdify(x , faux)
    
    #se obtienen los extremos
    a = xv[0]
    b = xv[-1]
    
    #se obtiene el máximo en ese intervalo
    x_max = fminbound(faux_n, a, b)

    fn4 = lambdify(x, f4)

    #se calcula la cota con la fórmula
    cota = (5*h**4/384)*fn4(x_max)

    return cota


if __name__ == "__main__":
    
    #se define la función y el vector de puntos a utilizar
    f = 'exp(x/2)'
    xv = np.array([1, 1.5, 1.75, 2.15, 2.4, 3], dtype=float)
    n = len(xv)

    #se llama a la función
    cota = cota_traz_cubico(f, n, xv)

    print("La cota de error es: ", cota)