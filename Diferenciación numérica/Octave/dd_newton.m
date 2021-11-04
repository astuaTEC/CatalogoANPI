
##function dd_newton_aux
##  clc; clear; %se limpia la consola
##  format long; %se configura el formato long para el resultado numérico
##  warning off; %se desactivan los mensajes de advertencia
##  
##  xv = [-2 0 1];
##  yv = [0 1 -1];
##  
##  polinomioInterpolacion = dd_newton(xv, yv)
##end


function polinomioInterpolacion = dd_newton(xv, yv)
   pkg load symbolic;
   syms x;
    
   n = length(xv);
    
   if length(xv) != length(yv)
     display("Error, el tamano de los vectores xv y yv es diferente.");
     return
   endif
  
   # Se calcula el primer termino del polinomio de interpolacion.
   polinomioInterpolacion = yv(1);

   variable = 1;

   # Cantidad de operaciones a realizar por cada iteracion.
   contador = n - 1;
   
   for i=2:n
     # Se van almacenando las variables (x-x0)(x-x1)...(x-xn-1).
     variable = variable * (x-xv(i-1));
     
     # Lista para almacenar los multiplicadores.
     multiplicadores = [];
     
     for j=1:contador
      numerador = yv(j+1) - yv(j);
      denominador = xv(j+i - 1) - xv(j);
      multiplicadores = [multiplicadores numerador/denominador];
     endfor
     
     # Se reduce en uno la cantidad de operaciones a realizar para 
     # la siguiente iteracion.
     contador = contador - 1;

     # Se le agrega el nuevo termino al polinomio de interpolacion.
     polinomioInterpolacion = polinomioInterpolacion + multiplicadores(1)*variable;

     # Se actualiza la lista en y.
     yv = multiplicadores;
     
   endfor

   polinomioInterpolacion = expand(polinomioInterpolacion);
   
end