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

symbol x("x");
ex func;

double f(double x_){
    ex xk = x_;
    ex result = evalf(subs(func, x==xk));
    return ex_to<numeric>(result).to_double();
}

void plot(vector<double> x, vector<double> y){
    plt::named_plot("Error |f(xk)|", x, y);
    plt::title("Error falsa posicion");
    plt::legend();
    plt::show();
}

vector<double> puntoFijo(string g, double x0, double tol, int iterMax){
    symtab table;
    table["x"] = x;
    parser reader(table);
    func = reader(g);
    vector<double> er, iter, resultado;

    int k = 0;
    double xk = x0;
    double error = 0;
    cout << func << endl;

    for(int i = 0; i < iterMax; i++){
        k +=1;
        iter.push_back(k);
        xk = f(xk);
        error = abs(f(xk));
        er.push_back(error);

        if(error < tol){
            k +=1;
            iter.push_back(k);
            break;
        }
    }

    resultado.push_back(xk);
    resultado.push_back(error);
    plot(iter, er);
    return resultado;
}   

int main(){

    double x0 = 1.5;
    string g = "log(2*x + 1)";
    double tol = 10e-9;
    int iterMax = 100;

    vector<double> res = puntoFijo(g, x0, tol, iterMax);
    cout << "Aproximación:  " + to_string(res.at(0)) << endl;
    cout << res.at(1) << endl;
    return 0;
}