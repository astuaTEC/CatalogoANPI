#################################
### METODO DE NEWTON RAPHSON ####
###                          ####
#################################


#Función auxiliar para llamar a la función de newton_raphson
#   En esta función se definen los atributos que la funcion
#   newton_raphson necesita para calcular la aproximación

function newton_raphson_prueba
  clc; clear; #se limpia la consola
  format long; #se configura el formato long para el resultado numérico
  
  f = 'exp(x) + x-2'; # se define la función a utilizar.
  x0 = 5; #se define el punto incial a utilizar
  tol = 10^-9; #se define la tolerancia máxima aceptada
  iterMax = 100; #se definen las iteraciones máximas antes de parar
  [xk error] = newton_raphson(f, x0, tol, iterMax) #se llama a la función y se muestra el resultado
  
end

# Metodo numérico newton_raphson para encontrar al menos un cero de una función f a partir de un valor inicial dado
# Además esta función grafica las iteraciones vs el error.
# f: la funcion a aplicarle el método
# x0: el punto inicial para comenzar el método
# tol: valor de la tolerancia de resultado aceptable
# iterMax: cantidad máxima de iteraciones que se pueden realizar
#
# Retorno: un vector con la aproximación y el error.

function [xk error] = newton_raphson(f, x0, tol, iterMax)
  %cargar el paquete symbolic
  pkg load symbolic
  
  f1 = matlabFunction(sym(f)); #se pasa la función escrita a una función que octave pueda operar
  df1 = matlabFunction(diff(sym(f))); #se le calcula la derivada a la función ingresada y se pasa a una función operable en octave
  er = []; #se crea un vector para los errores
  iter = []; #se crea un vector con las iteraciones
  error = tol + 1; #se define un valor inicial para el error
  k = 0; #se inicializa la variable de la iteración
  xk = x0; #se asigna el valor inicial a xk
  while (error > tol && k < iterMax)
    iter = [iter; k]; #se agrega k al vector de las iteraciones
    k = k + 1; #se suma una iteración
    n = f1(xk); #se calcula el numerador de la fórmula de newton raphson
    d = df1(xk); #se calcula el denominador de la fórmula de newton raphson
    
    if d > tol #se verifica que el denominador esté dentro de la tolerancia
      xk = xk - n/d; #se calcula el nuevo valor de la aproximación
      error = abs(f1(xk)); #se calcula el error
      er = [er; error]; #se inserta el error en el vector de errores
    else
      break
    endif
  endwhile
  
  plot(iter, er) #se grafica las iteraciones vs el error
  grid on #se activa la cuadrícula en la gráfica
  
end
