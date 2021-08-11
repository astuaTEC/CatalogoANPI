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

symbol x("x");
ex func;

double f(double x_){
    ex xk = x_;
    ex result = evalf(subs(func, x==xk));
    return ex_to<numeric>(result).to_double();
}

void plot(vector<double> x, vector<double> y){
    plt::named_plot("Error |f(xk)|", x, y);
    plt::title("Error Secante");
    plt::legend();
    plt::show();
}

vector<double> secante(string f_, double x0, double x1, double tol, int iterMax){
    symtab table;
    table["x"] = x;
    parser reader(table);
    func = reader(f_);
    double error = tol + 1;

    vector<double> er, iter, resultado;

    int k = 0;
    double xk = x1;
    double xk_1 = x0;
    double n, d;
    while (error > tol && k < iterMax){
        k += 1;
        iter.push_back(k);
        n = f(xk)*(xk - xk_1);
        d = f(xk) - f(xk_1);

        if(abs(d) > tol){
            xk_1 = xk;
            xk = xk - n/d;
            error = abs(f(xk));
            er.push_back(error);
        }
        else{
            break;
        }
    }

    resultado.push_back(xk);
    resultado.push_back(error);
    plot(iter, er);
    return resultado;
}


int main(){
    string f = "exp(-x^2) - x";
    double x0 = 0;
    double x1 = 1;
    double tol = 10e-9;
    int iterMax = 100;

    vector<double> res = secante(f, x0, x1, tol, iterMax);
    cout << "Aproximación:  " + to_string(res.at(0)) << endl;
    cout << res.at(1) << endl;
    
    return 0;
}