#ifndef MATRIX_GAUSS_H_INCLUDED
#define MATRIX_GAUSS_H_INCLUDED

#include "matrix.h"

// решает СЛУ A*X = B, где A - матрица коэффициентов, квадратная невырожденная, B - матрица правых частей
matrix *matrix_solve_gauss(const matrix *A, const matrix *B);

// проверка решения
double matrix_check_solution(const matrix *A, const matrix *B, const matrix *X);

#endif // MATRIX_GAUSS_H_INCLUDED
