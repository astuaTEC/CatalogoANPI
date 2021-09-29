/*
Comandos para la compilacion
 [1] $ g++ -std=c++14 -W -o th.out thomas.cpp -O2 -DARMA_DONT_USE_WRAPPER -larmadillo -llapack -lblas
 [2] $ ./th.out
*/

#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

/**
 * @brief Método para verificar si una matriz
 * es tridiagonal
 * @param A una matriz cuadrada
 * @param n el tamanio de la matriz
 * @return true si cumple, false si no
*/
bool verificarTridiagonal(mat A, int n){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double elem = A(i, j);
            if ((i == j) || (i-1 == j) || (i+1 == j)){
                if (elem == 0){
                    return false;
                }
            }
            else{
                if (elem != 0){
                    return false;
                }
            }
        }
    }
    return true;
}

/**
 * @brief Método para obtener los vectores
 * a, b y c de la matriz tridiagonal
 * @param A la matriz tridiagonal
 * @param n el tamanio de la matriz
 * @return un vector con los vectores a, b, c
*/
vector<colvec> obtener_abc(mat A, int n){
    vec a(n);
    a.zeros();
    vec b(n);
    b.zeros();
    vec c(n);
    c.zeros();

    vector<vec> result;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if(j == i){
                b(i) = A(i, j);
            }
            else if(i+1 == j){
                c(i) = A(i, j);
            }
            else if(i-1 == j){
                a(i) = A(i, j);
            }
            else{
                if (A(i, j) != 0){
                    cout << "La matriz no es tridiagonal" << endl;
                    throw exception();
                }      
            }
        } 
    }

    result.push_back(a);
    result.push_back(b);
    result.push_back(c);

    return result; 
}

/**
 * @brief Método auxiliar de thomas para trabajar
 * con los vectores a, b, c y d
 * @param a vector a
 * @param b vector b
 * @param c vector c
 * @param d vector d
 * @param n tamanio de los vectores
 * @return la solución del sistema xk
*/
vec thomas_aux(vec a, vec b, vec c, vec d, int n){
    vec p(n);
    p.zeros();
    vec q(n);
    q.zeros();
    vec xk(n);
    xk.zeros();

    p(0) = c(0)/b(0); //Primer coeficiente

    for (int i = 1; i < n-1; i++)
    {
        if( (b(i)- p(i-1)*a(i)) == 0){ //Comprueba que el divisor no sea 0
            cout << "No se puede dividir entre 0" << endl;
            throw exception();
        }
        else{
            p(i)= c(i)/(b(i)- p(i-1)*a(i)); //Calcula los nuevos coeficientes
        }
    }

    q(0) = d(0)/b(0);

    for (int i = 1; i < n; i++)
    {
        if ((b(i) - p(i-1)*a(i)) == 0) { //Comprueba que el divisor no sea 0
            cout << "No se puede dividir entre 0" << endl;
            throw exception();
        }
        else{
            q[i] = (d(i) - q(i-1)*a(i))/(b(i) - p(i-1)*a(i));
        }
    }
    
    xk(n-1) = q(n-1);
    
    for (int i = n-2; i >= 0; i--)
    {
        xk(i) = q(i) - p(i)*xk(i+1);
    }

    return xk; 
}

/**
 * @brief Función para resolver sistemas de ecuaciones
 * mediante el método directo de Thomas
 * @param A una matriz cuadrada
 * @param b vector de términos independientes
 * @return la solución del sistema
*/
vec thomas(mat A, vec b){
    int n = A.n_rows;
    if(!verificarTridiagonal(A, n)){
        cout << "La matriz debe ser tridiagonal" << endl;
        throw exception();
    }

    vector<vec> ABC = obtener_abc(A, n);

    return thomas_aux(ABC.at(0), ABC.at(1), ABC.at(2), b, n);
}

int main(){

    mat A = {{5, 1, 0, 0},
             {1, 5, 1, 0},
             {0, 1, 5, 1},
             {0, 0, 1, 5}};

    vec b = {-12, -14, -14, -12};

    vec sol = thomas(A, b);

    // Presentar la solución del Sistema.
    sol.print("x: ");

    return 0;
}