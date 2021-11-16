#include "matplotlibcpp.h"
#include <armadillo>
#include <math.h>

/*
Comandos para la compilacion
 [1] $ g++ -std=c++14 -W -o pi.out potencia_inversa.cpp -I/usr/include/python3.8 -lpython3.8 -O2 -DARMA_DONT_USE_WRAPPER -larmadillo -llapack -lblas
 [2] $ ./pi.out
*/

using namespace std;
using namespace arma;

namespace plt = matplotlibcpp;

/**
 * Grafica el error en funcion de la cantidad de iteraciones
 * @param x: El set de valores del eje x (iteraciones)
 * @param y: El set de valores del eje y (error)
 */ 
void plot(vector<double> x, vector<double> y){
    plt::named_plot("Error |xk_n - xk|", x, y);
    plt::title("Error potencia inversa");
    plt::legend();
    plt::show();
}

/**
 * Funcion para aproximar proxima el modulo del valor propio de 
 * menor magnitud deuna matriz A y el vector propio normalizadoasociado 
 * a dicho valor propio.
 * @param A: Matriz de tamanio mxm
 * @param x0: vector inicial
 * @param iterMax: iteraciones maximas
 * @param tol: tolerancia maxima
 * @return valor propio y vector propio respectivo
 */
tuple<double, vec> potencia_inversa(mat A, vec x0, int iterMax, double tol){

    vec xk = x0;

    vec yk, xk_n;

    double error, ck;

    vector<double> er, iter;
    
    // se calcula la aproximacion
    // mediante la formula de la
    // potencia inversa
    for (int i = 0; i < iterMax; i++)
    {
        //en lugar de calcular la inversa
        // se resuelve el sistema de ecuaciones
        yk = solve(A, xk);

        ck = norm(yk, "inf");

        xk_n = (1.0/ck)*yk;

        error = norm(xk_n - xk);

        xk = xk_n;

        iter.push_back(i);
        er.push_back(error);

        if(error < tol){
            break;
        }

    }

    // se grafica
    plot(iter, er);
    return tuple<double, vec>{ck, xk}; 
}

int main(){

    mat A = {{3, -1, 0},
            { -1, 2, -1},
            { 0, -1, 3}};

    vec x0 = {1, 1, 1};

    int iterMax = 9;

    double tol = 10e-10;

    tuple<double, vec> res = potencia_inversa(A, x0, iterMax, tol);

    cout << "Valor propio mayor: " + to_string(get<0>(res)) << endl;

    get<1>(res).raw_print("Vector propio respectivo: ");

    return 0;
}