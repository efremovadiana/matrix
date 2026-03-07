#ifndef MATRIX_GAUSS_H_INCLUDED
#define MATRIX_GAUSS_H_INCLUDED

#include "matrix.h"

/* Решение СЛАУ методом Гаусса, A -  матрица коэффициентов (квадратная), B - мтарица правых частей
 * tol - допуск для определения вырожденности
 * return - матрица X решений (нужно освободить через matrix_free) или NULL при ошибке
 */

matrix *matrix_solve_gauss(const matrix *A, const matrix *B, double tol);


 /* Вычисление невязки ||A*X - B||
  * A матрица коэффициентов, B матрица правых частей, X найденное решение
  */

double matrix_check_solution(const matrix *A, const matrix *B, const matrix *X);

#endif // MATRIX_GAUSS_H_INCLUDED

