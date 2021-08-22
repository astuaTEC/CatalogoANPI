/*
Comandos para la compilacion
 [1] $ g++ -std=c++14 -W -o pf.out puntoFijo.cpp -I/usr/include/python3.8 -lpython3.8 -lcln -lginac
 [2] $ ./pf.out
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
    plt::title("Error falsa posicion");
    plt::legend();
    plt::show();
}

/**
 * Metodo numérico punto fijo iterativo para encontrar al menos un cero de una función
 * g(x) que previamente ha cumplido con la existencia y unicidad
 * @param g: la funcion a aplicarle el método
 * @param x0: un punto inicial que cumpla los requisitos del punto fijo
 * @param tol: valor de la tolerancia de resultado aceptable
 * @param iterMax: cantidad máxima de iteraciones que se pueden realizar
 * @return un vector con el la aproximación a la solución de g(x) y un error asociado
 */
vector<double> puntoFijo(string g, double x0, double tol, int iterMax){
    // se definen aspectos de funcionamiento
    // de la librería GiNaC
    symtab table;
    table["x"] = x;
    parser reader(table);

    func = reader(g); // se pasa de la función en texto a una función que se puede evaluar
    vector<double> er, iter, resultado; //er: vector de errores, iter: vector de iteraciones, resultado: vector para retornar el resultado

    int k = 0; //se inicializa el las iteraciones en cero
    double xk = x0; //se le asigna el valor inicial a la variable de la aproximacion
    double error = 0; // se inicializa el error en cero

    for(int i = 0; i < iterMax; i++){
        k +=1; // se aumenta en 1 las iteraciones
        iter.push_back(k); //se introduce dentro del vector de las iteraciones
        xk = f(xk); // se evalúa la funcion
        error = abs(f(xk)); // se calcula el error
        er.push_back(error); //se introduce en el vector de errores

        if(error < tol){ // evalúa si ya existe un error aceptable
            break;
        }
    }

    resultado.push_back(xk); //se inserta el valor de la aproximacion
    resultado.push_back(error); // se inserta el valor del error asociado
    plot(iter, er); // se grafica el error vs iteraciones
    return resultado; // se retorna el resultado
}   

int main(){

    double x0 = 1.5; //se define el punto incial para el algoritmo
    string g = "log(2*x + 1)"; //se define la función a encontrarle al menos un cero
    double tol = 10e-9; //se define la tolerancia máxima
    int iterMax = 100; //se definen las iteraciones máximas

    vector<double> res = puntoFijo(g, x0, tol, iterMax); //se llama a la función y 
                                          //el resultado se guarda en una variable
    
    cout << "Aproximación:  " + to_string(res.at(0)) << endl;
    cout << "Error:  " + to_string(res.at(1)) << endl;
    return 0;
}