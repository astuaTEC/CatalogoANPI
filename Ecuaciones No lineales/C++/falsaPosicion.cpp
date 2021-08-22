/*
Comandos para la compilacion
 [1] $ g++ -std=c++14 -W -o fp.out falsaPosicion.cpp -I/usr/include/python3.8 -lpython3.8 -lcln -lginac
 [2] $ ./fp.out
*/

#include <iostream>
#include "matplotlibcpp.h"
#include <math.h>
#include <cmath>
#include <vector>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;
namespace plt = matplotlibcpp;

symbol x("x"); //se inicializa la variable simbólica
ex func; //se inicializa la funcion matemática a trabajar (GiNaC)

/**
 * Funcion para evaluar la funcion matemática a trabajar
 * @param x_: el valor numérico a sustituir en la funcion
 * @return el valor de evaluar la función matemática
 */
double f(double x_){
    ex xk = x_;
    ex result = evalf(subs(func, x==xk));
    return ex_to<numeric>(result).to_double();
}

/**
 * Grafica el error en funcion de la cantidad de iteraciones
 * @param x: El set de valores del eje x (iteraciones)
 * @param y: El set de valores del eje y (error)
 */ 
void plot(vector<double> x, vector<double> y){
    plt::named_plot("Error |f(xk)|", x, y);
    plt::title("Error falsa posición");
    plt::legend();
    plt::show();
}

/**
 * Metodo numérico falsa posicion para encontrar al menos un cero de una 
 * función ingresada f_
 * @param f_: la funcion a aplicarle el método
 * @param a: punto menor del conjunto a evaluar [a, b]
 * @param b: punto mayor del conjunto a evaluar [a, b]
 * @param tol: valor de la tolerancia de resultado aceptable
 * @param iterMax: cantidad máxima de iteraciones que se pueden realizar
 * @return un vector con el la aproximación a la solución de f(x) y un error asociado
 */
vector<double> falsa_posicion(string f_, double a, double b, double tol, int iterMax){
    // se definen aspectos de funcionamiento
    // de la librería GiNaC
    symtab table;
    table["x"] = x;
    parser reader(table);

    func = reader(f_); // se pasa de la función en texto a una función que se 
                        //puede evaluar
    vector<double> er, iter, resultado; //er: vector de errores, 
                                        //iter: vector de iteraciones, 
                                        //resultado: vector para retornar el resultado

    int k = 0; //se inicializa las iteración en cero
    double xk = a; //se actualiza el valor de xk con el límite inferior
    double xk_1 = b; //se actualiza el valor de xk_1 con el límite superior
    double error = 0; // se inicializa el error en cero

    if(f(a)*f(b) < 0){ //se verfica el teorema de Bolzano
        double n = f(xk_1) * (xk_1 - xk); //numerador formula de la secante
        double d = f(xk_1) - f(xk); //denominador de la formula de la secante
        xk = xk_1;
        xk_1 = xk_1 - n/d;
        
        while (k < iterMax){
            k++;
            iter.push_back(k); //se agrega la iteración al vector
            if(abs(d) > tol){ //condicion para la tolerancia en el denominador
                if(f(a)*f(xk_1) < 0){ //se verifica el teorema de Bolzano
                    //se cumple en el primer intervalo
                    b = xk_1;
                    n = f(xk_1) * (xk_1 - a); //numerador formula de la secante
                    d = f(xk_1) - f(a); //denominador de la formula de la secante
                    xk = xk_1; //se actualiza el valor de la iteración anterior
                    xk_1 = xk_1 - n/d; //se calcula el nuevo valor de la aproximación

                    error = abs((xk_1 - xk)/xk_1); //se calcula el error
                    er.push_back(error); //se inserta el error en el vector de errores
                    //se verifica que el error cumpla con la tolerancia
                    if(error < tol){ 
                        break;
                    }
                }
                else if(f(b)*f(xk_1) < 0){ //se verifica el teorema de Bolzano
                    //se cumple en el segundo intervalo
                    a = xk_1;
                    n = f(xk_1) * (xk_1 - b); //numerador formula de la secante
                    d = f(xk_1) - f(b); //denominador de la formula de la secante
                    xk = xk_1; //se actualiza el valor de la iteración anterior
                    xk_1 = xk_1 - n/d; //se calcula el nuevo valor de la aproximación

                    error = abs((xk_1 - xk)/xk_1); //se calcula el error
                    er.push_back(error); //se inserta el error en el vector de errores

                    //se verifica que el error cumpla con la tolerancia
                    if(error < tol){ 
                        break;
                    }
                }
                else{
                    break;
                }
            }
            else{
                break;
            }
        }
        resultado.push_back(xk_1); //se inserta el valor de la aproximacion
        resultado.push_back(error); // se inserta el valor del error asociado
        plot(iter, er); // se grafica el error vs iteraciones
        return resultado; // se retorna el resultado
    }
    else{
        cout << "No se cumple el teorema de Bolzano" << endl;
        resultado.push_back(0);
        resultado.push_back(0);

        return resultado;
    }

}


int main(){
    double a = 1; //se define el límite inferior del rango
    double b = 2; //se define el límite superior del rango
    string f = "log(x) - exp(-x) - cos(x)"; //se define la función a calcularle 
                                            //al menos un cero
    double tol = 10e-5; //se define la tolerancia
    int iterMax = 100; //se definen las iteraciones máximas

    vector<double> res = falsa_posicion(f, a, b, tol, iterMax); //se llama al método 
                                            //y el resultado se guarda en una variable
    
    cout << "Aproximación:  " + to_string(res.at(0)) << endl;
    cout << "Error:  " + to_string(res.at(1)) << endl;

    return 0;
}


