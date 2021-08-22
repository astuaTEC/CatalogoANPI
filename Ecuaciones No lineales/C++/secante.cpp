/*
Comandos para la compilacion
 [1] $ g++ -std=c++14 -W -o sec.out secante.cpp -I/usr/include/python3.8 -lpython3.8 -lcln -lginac
 [2] $ ./sec.out
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
    plt::title("Error Secante");
    plt::legend();
    plt::show();
}

/**
 * Metodo numérico de la secante iterativo para encontrar al menos un cero 
 * de una función
 * @param f_: la funcion a aplicarle el método
 * @param x0: un punto inicial que cumpla que esté dentro de los rangos [a, b]
 * @param x1: un segundo punto que cumpla esté dentro de los rangos [a, b]
 * @param tol: valor de la tolerancia de resultado aceptable
 * @param iterMax: cantidad máxima de iteraciones que se pueden realizar
 * @return un vector con el la aproximación a la solución de f(x) y un error asociado
 */
vector<double> secante(string f_, double x0, double x1, double tol, int iterMax){
    // se definen aspectos de funcionamiento
    // de la librería GiNaC
    symtab table;
    table["x"] = x;
    parser reader(table);

    func = reader(f_); // se pasa de la función en texto a una función que 
                       // se puede evaluar
    double error = tol + 1; // se define una tolerancia

    vector<double> er, iter, resultado; //er: vector de errores, 
                                        //iter: vector de iteraciones,
                                        // resultado: vector para retornar el resultado

    int k = 0; //se inicializa el las iteraciones en cero
    double xk = x1; //se le asigna el valor inicial a la variable de la aproximacion
    double xk_1 = x0; //se le asigna un valor a la iteración siguiente.
    double n, d; //se define el numerador y el denominador
    while (error > tol && k < iterMax){
        k += 1; //se suma uno a las iteraciones
        iter.push_back(k); //se inserta een el vector de las iterciones
        n = f(xk)*(xk - xk_1); //se calcula el numerador de la fórmula
        d = f(xk) - f(xk_1); // se calcula el denominador de la fórmula

        if(abs(d) > tol){ // se verifica que  el denominador no se indefina
            // se actualizan valores
            xk_1 = xk; 
            xk = xk - n/d;
            error = abs(f(xk));
            er.push_back(error);
        }
        else{
            break;
        }
    }

    resultado.push_back(xk); //se inserta el valor de la aproximacion
    resultado.push_back(error); // se inserta el valor del error asociado
    plot(iter, er);
    return resultado;
}


int main(){
    //se definen valores para llamar al método
    string f = "exp(-x^2) - x"; 
    double x0 = 0;
    double x1 = 1;
    double tol = 10e-9;
    int iterMax = 100;

    vector<double> res = secante(f, x0, x1, tol, iterMax); //se llama a la función
    cout << "Aproximación:  " + to_string(res.at(0)) << endl;
    cout << res.at(1) << endl;
    
    return 0;
}