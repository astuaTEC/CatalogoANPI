%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Metodo de punto fijo%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Metodo de punto fijo, el cual calcula almenos un creo segun  los parametros  
% de entrada
function punto_fijo 
  clc; clear; %Limpia consola
  format long; %Permite retornar un valor con muchas decimales
  f = 'log(2*x + 1)'; %Funcion entrada
  x0 = 1.5; %Valor inicial
  tol=10^-9; % toleracia de erroe
  iterMax=150; %Maxima interacion
  [xk error] = punto_fijo_aux(f, x0, tol, iterMax)% Funcion Axiliar
end

function [xk error] = punto_fijo_aux(f, x0, tol, iterMax)
  %Se inicalisan las distintas variables
  er = [];
  error = tol;
  iter = [];
  xk = x0;
  f1 = matlabFunction(sym(f));%Se convierte el string de F en una funcion simbolica
  for(k=0:iterMax)%For para empezar las iteraciones
    iter = [iter; k];
    k +=1; %Se aumenta en uno la condicion de parada
    xkaux = xk;
    xk = f1(xk); %Se calcula un nuvo valor para la iteracion 
    error = abs((xk-xkaux)/xk);% Se calcula el error
    er = [er; error]; %Se guardan valores para graficar el error 
    if(error < tol) %En caso de llega a una solucion valida se rompe la iteracion
      break;
    endif
  endfor
  %Mensajes para graficar
  plot(iter, er)
  grid on
end


