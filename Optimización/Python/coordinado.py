import numpy as np
import sympy as sp
from scipy import optimize
import matplotlib.pyplot as plt

"""
Metodo para convertir una serie de variables
ingresadas a variables simbolicas
variables: las varianbles en un string
            en formato 'x y z'
Retorna: una lista con las variables simbolicas
"""
def variables_simbolicas(variables):
    n = len(variables)
    tam = np.arange(0, n, 2)
    v_var = []
    for i in tam:
        v_var.append(sp.Symbol(variables[i]))
    return v_var

"""
Metodo para calcular el gradiente de una funcion dada
y retornarlo en su forma simbolica
f: la funcion a calcularle su gradiente
v_var: variables simbolicas para calcular el gradiente
Retorna: el gradiente en forma simbolica
         dentro de una lista
"""
def gradiente(f, v_var):
    n = len(v_var)
    g = []
    for i in np.arange(0, n, 1):
        g.append(sp.diff(f, v_var[i]))
    return g


"""
Metodo numérico descenso coordinado para encontrar 
al menos un minimo de una función f multivariable
a partir de un valor inical (x0, y0, ...).
Además esta función grafica las iteraciones vs el error.
funcion: la funcion a aplicarle el método de optimizacion
variables: Representan todas las variables simbolicas
           con las que se va a trabajar.
puntoInicial: punto de referencia inicial para comenzar
              a aplicar el metodo
tol: valor de la tolerancia de resultado aceptable
iterMax: cantidad máxima de iteraciones que se pueden realizar
Retorno: un vector con la aproximación y el error.
"""
def coordinado(funcion, variables, puntoInicial, iterMax, tol):
    #funciones
    f_s = sp.sympify(funcion) #simbolica

    #gradirente (simbolico)
    g = gradiente(f_s, variables)

    #punto inicial
    xk = puntoInicial

    #se accede al numero de variables
    numVariables = len(variables)

    #vector para guardar el error
    er = []

    k = 0
    while k < iterMax:
        k += 1
        n = 0
        indexSuperior = numVariables - 1
        
        #por cada variable se "fijan" las demas
        # y se calcula el nuevo punto xk
        while n < numVariables:
            #se usa el metodo de Gauss-Seidel
            simbolica = f_s.subs(variables[indexSuperior], 
                        xk[indexSuperior])

            funcion_num = sp.lambdify(variables[n], simbolica)

            punto = optimize.fmin(funcion_num, puntoInicial[n])

            xk[n] = punto[0]

            n += 1

            indexSuperior -= 1
        
        #se utiliza para evaluar el gradiente
        #de acuerdo al numero de variables existente
        data ={}
        g_evaluado = []
        for i in range(0, numVariables):
            data[variables[i]] = xk[i]

        for i in range(0, numVariables):
            g_evaluado.append(float(g[i].subs(data).evalf()))

        #se calcula el error y se agrega
        #al vector de errores
        error = np.linalg.norm(g_evaluado)
        er.append(error)
        if error < tol:
            break

    fig, graf = plt.subplots() #se crea la gráfica
    ejeX = np.arange(1, k+1, 1) #se crea el eje X (son las iteraciones)
    graf.plot(ejeX, er) #se grafican los datos
    plt.show() #se muestra la gráfica
    return [xk, error]


if __name__ == "__main__":
    #se definen las variables
    variables = 'x y'
    #se crean las variables simbólicas
    v_var = variables_simbolicas(variables)
    #se define la función a trabajar
    f = '(x-2)^2+(y+3)^2 + x*y'
    #se define un punto inicial
    puntoInicial = [1, 1]
    #se definen las iteraciones máximas y la tolerancia
    iterMax = 9
    tol = 10e-6
    #se llama a la función
    coordinado(f, v_var, puntoInicial, iterMax, tol)
