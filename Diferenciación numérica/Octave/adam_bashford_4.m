

function adam_bashford_4_aux
  pkg load symbolic;
  clc; clear; close all
  warning off;
  
  intervalo = [2 4];
  pasoh = 11;
  funcion = '1 + (x-y)^2';
  yinicial = 1;
  
  [xv, yv, polinomioInterpolacion] = adam_bashford_4(funcion, intervalo, pasoh, yinicial)
  
endfunction

% f: Funcion del tipo string
% Entradas:
%	   intervalo: Lista con los valores frontera en x con la forma: [xi, xf]
%    pasoh: tamaño de paso 
%    yo: Valor inicial de la solucion para y(x=0) 

%Salidas:
% 	 xv: coordenadas x evaluadas
% 	 yv: coodernadas y obtenidas
% 	 polinomioInterpolacion: polinomio de interpolacion en formato simbolico
function [xv, yv, polinomioInterpolacion] = adam_bashford_4(funcion, intervalo, pasosh, y0)
  % Simbolico
  syms x;
  syms y;
  
  f_simbolica = sym(funcion);
  
  n = pasosh;
  
  f = matlabFunction(f_simbolica);
  
  a = intervalo(1);
  b = intervalo(2);
  
  h = (b-a)/(n-1);
  xv = a:h:b;
  
  %Se calulan las y con las formulas del metodo
  yv(1) = y0;
  yv(2) = yv(1) + h*f(xv(1), yv(1));
  yv(3) = yv(2) + h/2*(3*f(xv(2), yv(2)) - f(xv(1), yv(1)));
  yv(4) = yv(3) + h/12*(23*f(xv(3), yv(3)) - 16*f(xv(2), yv(2)) + 5*f(xv(1), yv(1)));
  
  for i=4:n-1
      yv(i+1) = yv(i) + h/24*(55*f(xv(i), yv(i)) 
      - 59*f(xv(i-1), yv(i-1)) + 37*f(xv(i-2), yv(i-2)) - 9*f(xv(i-3), yv(i-3)));
  endfor 
  
  polinomioInterpolacion = dd_newton(xv, yv);
  poli_n = matlabFunction(polinomioInterpolacion);
  
  hold on
  ezplot(poli_n, intervalo);
  stem(xv, yv, 'b');
  
endfunction 

% Funcion o
%
%
function polinomioInterpolacion = dd_newton(xv, yv)
   pkg load symbolic;
   syms x;
    
   n = length(xv);
    
   if length(xv) != length(yv)
     display("Error, el tamano de los vectores xv y yv es diferente.");
     return
   endif
  
   % Se calcula el primer termino del polinomio de interpolacion.
   polinomioInterpolacion = yv(1);

   variable = 1;

   % Cantidad de operaciones a realizar por cada iteracion.
   contador = n - 1;
   
   for i=2:n
     % Se van almacenando las variables (x-x0)(x-x1)...(x-xn-1).
     variable = variable * (x-xv(i-1));
     
     % Lista para almacenar los multiplicadores.
     multiplicadores = [];
     
     for j=1:contador
      numerador = yv(j+1) - yv(j);
      denominador = xv(j+i - 1) - xv(j);
      multiplicadores = [multiplicadores numerador/denominador];
     endfor
     
     % Se reduce en uno la cantidad de operaciones a realizar para 
     % la siguiente iteracion.
     contador = contador - 1;

     % Se le agrega el nuevo termino al polinomio de interpolacion.
     polinomioInterpolacion = polinomioInterpolacion + multiplicadores(1)*variable;

     % Se actualiza la lista en y.
     yv = multiplicadores;
     
   endfor

   polinomioInterpolacion = expand(polinomioInterpolacion);
   
end

