%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METODO GRADIENTE CONJUGADO NO LINEAL%
%                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Funcion auxiliar para llamar a la funcion de gradiente
%   En esta función se definen los atributos que la funcion
%   gradiente necesita para calcular la minimizacion
function gradiente_aux
  clc; clear; %se limpia la consola
  format long; %se configura el formato long para el resultado numÃ©rico
  warning off; %se desactivan los mensajes de advertencia
  
  f = '(x - 2)^4 + (x - 2*y)^2'; % se define la funciÃ³n a utilizar.
  x0 = [0, 3]; %se define el punto inicial
  variables = ['x', 'y']; %se definen las variables a utilizar
  tol = 10^-3; %se define la tolerancia
  iterMax = 13; %se define el maximo de iteraciones
  
  %se llama a la funcion
  [xk error] = gradiente(f, x0, variables, iterMax, tol)
  
end


% Metodo numérico gradiente conjugado no lineal
% para encontrar al menos un punto minimo en una
% función f multivariable a partir de un valor inicial dado
% Además esta función grafica las iteraciones vs el error.
% f: la funcion a aplicarle el método
% x0: el punto inicial para comenzar el método
% variables: las variables con las que se trabaja la funcion f
% tol: valor de la tolerancia de resultado aceptable
% iterMax: cantidad máxima de iteraciones que se pueden realizar
%
% Retorno: un vector con la aproximación y el error.
function [xk error] = gradiente(f, x0, variables, iterMax, tol)
  %cargar el paquete symbolic
  pkg load symbolic
  
  f = sym(f); %se pasa la funciÃ³n escrita a una funcion simbolica
  f1 = matlabFunction(f); %Funciï¿½n f en formato del lenguaje M
  er = []; %se crea un vector para los errores
  iter = []; %se crea un vector con las iteraciones
  error = tol + 1; %se define un valor inicial para el error
  
  xk = x0'; %se define el punto para la iteracion inicial
            % importante recalcar que es la transpuesta
            % del punto x0
            
  numVariables = length(variables); %se define el numero de
                                    %variables a trabajar
  symbolicas = []; %lista de variables simbolicas
  
  grad = gradient(f); %se calcula el gradiente
  
  %se convierte a cada una de las variables
  %ingresadas a variables simbolicas
  for n=1: numVariables
    symbolicas = [symbolicas, sym(variables(n))];
  endfor
  
   %se evalua el gradiente en el punto inicial
   % el resultado es un vector columna
  gk = [double(subs(grad, symbolicas, xk'))];
  
  %se define el vector dk
  % este es un vector fila
  dk = -gk;
  
  iterMax2 = 1000; %iteraciones maximas para el valor de alpha_k
  
  %se calcula el error inicial y se agrega al
  %vector de errores
  error = norm(double(subs(grad, symbolicas, xk')));
  er = [er, error];
  
  %se agrega la primera iteracion al vector
  % de iteraciones
  iter = [iter, 0];
  for k = 0: iterMax
    k = k + 1;
    iter = [iter, k];
    delta = 0.5; %se define delta
    alpha_k = 1; %se define alpha_k
    
    %primero se debe calcular el valor de alpha_k
    %para poder calcular el valor de la siguiente
    %iteracion
    %esto se hace mediante la siguiente validacion
    for i = 1: iterMax2
      
      %se define un vector para 
      %evaluar la funcion ingresada
      vectorEval = xk + alpha_k*dk;
      
      %se tiene una ecuacion en la que se tiene que evaluar 
      %que el lado izquierdo sea menor o igual que el lado derecho
      
      %se calcula el lado izquierdo
      izq = double(subs(f, symbolicas, vectorEval') 
            - subs(f, symbolicas, xk'));
            
      %se calcula el lado derecho
      der = delta*alpha_k*(gk)'*dk;
      if izq <= der
        break
      else
        alpha_k = alpha_k/2; %se define un valor mas
                             %pequeno para alpha_k
      endif
    endfor
    
    xk = xk + (alpha_k*dk); %se calcula la nueva aproximacion
    
    %se calcula el error y se agrega al
    % vector de errores
    error = norm([double(subs(grad, symbolicas, xk'))]);
    er = [er, error];
    if error < tol
      break
    endif
    
    % se calcula el valor de beta
    % esto se hace por medio del grediente
    % evaluado en la iteracion anterior
    % y el gradiente evaluado en el punto nuevo 
    gk_anterior = gk;
    gk = double(subs(grad, symbolicas, xk'));
    
    beta_k = (norm(gk)^2)/(norm(gk_anterior)^2);
    
    %se calcula el nuevo valor de dk
    dk = -gk + beta_k*dk;   
  endfor
  
  %se grafica las iteraciones vs el error
  plot(iter, er)
  
end