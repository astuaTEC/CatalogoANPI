function euler_aux
  pkg load symbolic;
  %% Ejemplo del M���todo de Euler
  clc; clear; close all
  intervalo = [0 5];
  num_pt = 11;
  funcion = 'y - x^2 + 1';
  
  [xv, yv] = euler(funcion, intervalo, num_pt)
  
   
endfunction


function [xv, yv] = euler(funcion, intervalo, n)
  % Simbolico
  syms x;
  f_simbolica = sym(funcion);
  
  f = matlabFunction(f_simbolica);
  
  a = intervalo(1);
  b = intervalo(2);
  
  h=(b-a)/(n-1);
  xv=a:h:b;
  yv=[0.5];
  for n=1:n-1  
    yv(n+1)=yv(n)+h*f(xv(n),yv(n));
  end
  
  hold on
  plot(xv, yv, 'r');
  
  %Soluci���n anal���tica
  y_s=@(x) (x+1).^2-0.5*exp(x);
  %Graficar la solucion
  
  x_g=a:0.0001:b;
  y_g=y_s(x_g);
  plot(x_g,y_g,'b')
  
endfunction 



