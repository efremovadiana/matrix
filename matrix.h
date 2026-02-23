#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <stddef.h>        // Для типа size_t (беззнаковый целый тип)
#include <stdio.h>          // Для FILE* (нужно для вывода)


typedef struct matrix matrix;

// новая матрица, w - количество столбцов, h - количество строк
matrix *matrix_alloc(size_t w, size_t h);

void matrix_free(matrix *m);

// копирует матрицу
matrix *matrix_copy(const matrix *m);

// получение указателя на элемент
double *matrix_ptr(matrix *m, size_t i, size_t j);

//константный доступ
const double *matrix_cptr(const matrix *m, size_t i, size_t j);

// количество строк матрицы, m - указатель на матрицу
size_t matrix_rows(const matrix *m);

// количество столбцов матрицы, m - указатель на матрицу
size_t matrix_cols(const matrix *m);

matrix *matrix_alloc_id(size_t n);

matrix *matrix_assign(matrix *m1, const matrix *m2);

// сложение матриц
int matrix_add(matrix *m1, const matrix *m2);

// умножение на число
void matrix_smul(matrix *m, double d);

// деление на число
void matrix_sdiv(matrix *m, double d);

// вычисляет норму матрицы
double matrix_norm(const matrix *m);

// умножение двух матриц, res - матрица результата, a, b - матрицы
int matrix_multiply(matrix *res, const matrix *a, const matrix *b);

// вывод
void matrix_print(const matrix *m, const char *format);

#endif // MATRIX_H_INCLUDED
