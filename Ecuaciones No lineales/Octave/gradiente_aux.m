############################################
### METODO GRADIENTE CONJUGADO NO LINEAL####
###                                     ####
############################################


function gradiente_aux
  clc; clear; #se limpia la consola
  format long; #se configura el formato long para el resultado numérico
  warning off;
  
  f = '(x - 2)^4 + (x - 2*y)^2'; # se define la función a utilizar.
  x0 = [0, 3];
  variables = ['x', 'y'];
  tol = 10^-3;
  iterMax = 13;
  
  [xk error] = gradiente(f, x0, variables, iterMax, tol)
  
  
  
end


function [xk error] = gradiente(f, x0, variables, iterMax, tol)
  %cargar el paquete symbolic
  pkg load symbolic
  
  f = sym(f); #se pasa la función escrita a una funcion simbolica
  f1 = matlabFunction(f); #Funci�n f en formato del lenguaje M
  er = []; #se crea un vector para los errores
  iter = []; #se crea un vector con las iteraciones
  error = tol + 1; #se define un valor inicial para el error
  xk = x0';
  
  
  numVariables = length(variables);
  symbolicas = []; #lista de variables simbolicas
  
  grad = gradient(f);
  
  for n=1: numVariables
    symbolicas = [symbolicas, sym(variables(n))];
  endfor
  
  gk = [double(subs(grad, symbolicas, xk'))];
  dk = -gk;
  
  iterMax2 = 1000; #iteraciones maximas para el valor de alpha_k
  
  error = norm(double(subs(grad, symbolicas, xk')));
  er = [er, error];
  iter = [iter, 0];
  for k = 0: iterMax
    k = k + 1;
    iter = [iter, k];
    delta = 0.5; #se define delta
    alpha_k = 1; #se define alpha_k
    
    for i = 1: iterMax2
      vectorEval = xk + alpha_k*dk;

      izq = double(subs(f, symbolicas, vectorEval') - subs(f, symbolicas, xk'));
      der = delta*alpha_k*(gk)'*dk;
      if izq <= der
        break
      else
        alpha_k = alpha_k/2;
      endif
    endfor
    
    xk = xk + (alpha_k*dk);
    error = norm([double(subs(grad, symbolicas, xk'))]);
    er = [er, error];
    if error < tol
      break
    endif
    gk_anterior = gk;
    gk = double(subs(grad, symbolicas, xk'));
    
    beta_k = (norm(gk)^2)/(norm(gk_anterior)^2);
    dk = -gk + beta_k*dk;   
  endfor
  
  
  plot(iter, er)
  
end