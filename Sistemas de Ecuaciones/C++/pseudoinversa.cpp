// [1] $ g++ -std=c++14 -W -o psi.out pseudoinversa.cpp -I/usr/include/python3.8 -lpython3.8 -O2 -DARMA_DONT_USE_WRAPPER -larmadillo -llapack -lblas
// [2] $ ./psi.out
// Nota: Asegurarse de tener armadillo bien instalado.

#include <iostream>
#include "matplotlibcpp.h"
#include <armadillo>

using namespace std;
using namespace arma;
namespace plt = matplotlibcpp;


/**
 * Grafica el error en funcion de la cantidad de iteraciones
 * @param x: El set de valores del eje x (iteraciones)
 * @param y: El set de valores del eje y (error)
 */ 
void plot(vector<double> x, vector<double> y){
    plt::named_plot("Error |f(xk)|", x, y);
    plt::title("Error pseudoinversa");
    plt::legend();
    plt::show();
}

void pseudoinversa(mat A, vec b, double tol, int iterMax=15){
    int n = A.n_rows;
    int m = A.n_cols;
    double alpha = eig_sym(A*A.t()).max();

    vector<double> er, iter;

    mat I(n, m);
    I.eye();
    mat Xk = (1/alpha)*A.t();
    vec xk = Xk*b;

    double error = 1;
    int k = 0;

    vec xk_n;
    
    iter.push_back(k);
    while (k < iterMax){
        Xk = Xk*(2*I-A*Xk);
        xk_n = Xk*b;

        error = norm(xk_n - xk)/norm(xk_n);
        er.push_back(error);
        xk = xk_n;

        if (error < tol)
            break;
        
        iter.push_back(k);
        k++;
    }

    // Mostrar los resultados
    xk.print("xk: \n");
    cout<<"Error del sistema: "<< norm(A*xk-b)<<endl;

    plot(iter, er);
}


int main()
{
    mat A = {{ 1, 2, 4},
             { 2,-1, 1},
             { 1, 0, 1}};
             
    vec b = {4, 3, 9};
    pseudoinversa(A,b,10e-5,1000);
    return 0;
}   