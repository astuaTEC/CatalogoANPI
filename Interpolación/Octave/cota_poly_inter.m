
function cota_poly_interpolacion
  pkg load symbolic;
  warning off;
  clc;
  clear;
  f = "sin(pi*x/2)";
  vector = [-1 0 1 2];
  n = length(vector);
  cota_de_error = cota_poly_inter(f,n,vector)
end

%Calcula la cota de error de polinomio de interpolacion a partir de la funcion 
%a la cual se le calcula el polinomio de interpolacion y un conjunto de valores
%
%Entradas:
%     f: Funcion a la cual se le va a calcular la cota de error
%                                          (con x como variable)
%     vector: Conjunto de puntos con los cuales se calcula la cota de error
%     n: cantidad de datos dentro del vactor
%
%Salida:
% Cota de error maxima del polinomio de interpolacion
%

function [cota_de_error] = cota_poly_inter(funcion,n,vector)
  % Simbolico
  syms x;
  f_simbolica = sym(funcion);
  a = vector(1);
  b = vector(n);
  

  %Se calcula el valor de alphamax 
  
  %Para sacar un máximo de la función en un intervalo es equivalente a sacarle 
  %el mínimo a la función negada en ese intervalo max(f) = min(-1*f)
  
  f_derivada_n = diff(f_simbolica,n);
  f_derivada_n_evaluable = matlabFunction(f_derivada_n);
  f_aux = -1 * abs(f_derivada_n); 
  f_aux_num = matlabFunction(f_aux);
  x_max = fminbnd(f_aux_num,a,b);
  alpha_max = abs(f_derivada_n_evaluable(x_max));
  
  %Se calcula el valor maximo del funcion del tipo (x-x1)(x-x2)...(x-xn)

  funcionMult = generarMultilplicadores(vector, n);
  funcionMult_evaluable = matlabFunction(funcionMult);
  funcionMult_aux = -1 * abs(funcionMult);
  funcionMult_aux_evaluable = matlabFunction(funcionMult_aux);
  funcionMult_xmax = fminbnd(funcionMult_aux_evaluable, a, b);
  funcionMult_max = abs(funcionMult_evaluable(funcionMult_xmax));
  
  
  cota_de_error = (1/(factorial(n)))*alpha_max*funcionMult_max;
end


% Funcion que genera una funcion simbolica del tipo (x-x1)(x-x2)...(x-xn) 
% a partir de un conjunto de vector
% Entradas: 
%      Conjunto soporte
% Salidas: 
%      Funcion simbolica del tipo (x-x1)(x-x2)...(x-xn)
function symFMul = generarMultilplicadores(vector)
  n = length(vector);
  symStr = "1";
  for i=1:n
    charS = num2str(vector(i));
    symStr = strcat(symStr, "*(x-",charS,")");
  end
  symFMul = sym(symStr);
end