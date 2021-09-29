/*
Comandos para la compilacion
 [1] $ g++ -std=c++14 -W -o flu.out fact_lu.cpp -O2 -DARMA_DONT_USE_WRAPPER -larmadillo -llapack -lblas
 [2] $ ./flu.out
*/

#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

/**
 * @brief  Método que realiza la sustitución hacia adelante 
 * en un sistema de ecuaciones
 * @param A una matriz cuadrada
 * @param b vector de términos independientes
 * @param n tamanio de la matriz
 */
vec sust_adelante(mat A, vec b, int n){
    vec soluciones(n);
    soluciones.zeros();

     for (int i = 0; i < n; i++){
        double sumatoria = 0;
        for (int j = 0; j < i; j++){
            sumatoria += A(i,j) * soluciones(j);
        }
        soluciones(i) = (b(i) - sumatoria) * (1/A(i,i));
    } 

    return soluciones;
}

/**
 * @brief  Método que realiza la sustitución hacia atrás 
 * en un sistema de ecuaciones
 * @param A una matriz cuadrada
 * @param b vector de términos independientes
 * @param n tamanio de la matriz
 * @return vector con la solución
 */
vec sust_atras(mat A, vec b, int n){
    vec soluciones(n);
    soluciones.zeros();

    for (int i = n-1; i >= 0; i--){
        double sumatoria = 0;
        for (int j = i+1; j < n; j++){
            sumatoria += A(i,j) * soluciones(j);
        }
        soluciones(i) = (1/A(i,i)) * (b(i) - sumatoria);
    } 

    return soluciones;
}

/**
 * @brief Verifica si la matriz tiene una única
 * factorización LU
 * @param A matriz cuadrada
 * @param n el tamanio de la matriz
 * @return true si cumple, false si no
 */
bool verificar_fact_lu(mat A, int n){
    double tol = datum::eps;
    for (int i = 0; i < n; i++)
    {
        mat Ak = A.submat(0,0,i,i);
        double deter = det(A);
        if (abs(deter) < tol){
            return false;
        }
    }
    return true;
    
}

/**
 * @brief Método para obtener la factorización LU
 * de una determinada matriz
 * @param A una matriz cuadrada
 * @param n el tamanio de la matriz
 * @return un vector con la matriz L y U
 */
vector<mat> obtener_lu(mat A, int n)
{
    mat U = A;
    mat L = eye(n, n);
    vector<mat> result;

    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double m_ik = U(i, k) / U(k, k);
            L(i, k) = m_ik; //se agrega el multiplicador

            for (int j = k; j < n; j++)
            {
                U(i, j) = U(i, j) - (m_ik * U(k, j));
            }
        }
    }

    result.push_back(U);
    result.push_back(L);

    return result;
}

/**
 * @brief Función para resolver sistemas de ecuaciones
 * mediante el método directo de factorización LU
 * @param A una matriz cuadrada
 * @param b vector de términos independientes
 * @return la solución del sistema
*/
vec fact_lu(mat A, vec b){
    int n = A.n_rows;
    if(!verificar_fact_lu(A, n)){ //se verifica el teorema de convergencia
        return "No cumple con los requisitos de la fact LU";
    }

    vector<mat> LU = obtener_lu(A, n);

    vec y = sust_adelante(LU.at(1), b, n); //resuelve Ly = b
    vec x = sust_atras(LU.at(0), y, n); //Resuelve Ux = y

    return x;
}

int main(){

    mat A = {{ 4, -2, 1},
            { 20, -7, 12},
            { -8, 13, 17}};

    vec b = {11, 70, 17};

    vec sol = fact_lu(A, b);

    // Presentar la solución del Sistema.
    sol.print("x: ");

    return 0;
}