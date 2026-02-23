#ifndef MATRIX_EXP_H_INCLUDED
#define MATRIX_EXP_H_INCLUDED

#include "matrix.h"

// Вычисление матричной экпоненты, где A - указатель на исходную квадратную матрицу, exp - точность вычислений
matrix *matrix_exp(const matrix *A, double exp);

#endif
