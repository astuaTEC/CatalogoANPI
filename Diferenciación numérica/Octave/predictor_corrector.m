%Predictor corrector
function predictor_corrector_aux
  pkg load symbolic;

  clc; clear; close all
  warning off;
  
  intervalo = [0 2];
  pasoh = 11;
  funcion = 'y - x^2 + 1';
  yinicial = 0.5;
  
  [xv, yv, polinomioInterpolacion] = predictor_corrector(funcion, intervalo, pasoh, yinicial)
  
endfunction

% f: Funcion del tipo string
% Entradas:
%	   intervalo: Lista con los valores frontera en x con la forma: [xi, xf]
%    pasoh: tama�o de paso 
%    yinicial: Valor inicial de la solucion para y(x=0) 

%Salidas:
% 	 xv: coordenadas x evaluadas
% 	 yv: coodernadas y obtenidas
% 	 polinomioInterpolacion: polinomio de interpolacion en formato simbolico
function [xv, yv, polinomioInterpolacion] = predictor_corrector(funcion, intervalo, pasoh, yinicial)
  % Simbolico
  syms x;
  f_simbolica = sym(funcion);
  
  n = pasoh;
  
  f = matlabFunction(f_simbolica);
  
  a = intervalo(1);
  b = intervalo(2);
  
  h=(b-a)/(n-1);
  xv=a:h:b;
  yv=[yinicial];
  
  %se calculan las y con las formulas
  for i=1:n-1 
    zv = yv(i)+h*f(xv(i),yv(i));
    yv(i+1)= yv(i) + h/2*(f(xv(i),yv(i)) + f(xv(i+1),zv));
  end
  
  polinomioInterpolacion = dd_newton(xv, yv);
  poli_n = matlabFunction(polinomioInterpolacion);
  
  hold on
  ezplot(poli_n, intervalo);
  stem(xv, yv, 'b');
  
endfunction 


% Funcion para calcular el polinomio de interpolacion, mediante el metodo
% de diferencias divididas de Newton.
% Entradas:       xv = vector con los puntos en x.
%                 yv = vector con los puntos en y.
%Salidas:        polinomioInterpolacion = polinomio de interpolacion
%                tras ser calculado.
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
