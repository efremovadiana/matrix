#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// функция умножения матриц
static int gauss_matrix_multiply(matrix *res, const matrix *a, const matrix *b) {
    if (res == NULL || a == NULL || b == NULL) return -1;

    size_t a_rows = matrix_rows(a);
    size_t a_cols = matrix_cols(a);
    size_t b_rows = matrix_rows(b);
    size_t b_cols = matrix_cols(b);
    size_t res_rows = matrix_rows(res);
    size_t res_cols = matrix_cols(res);

    if (a_cols != b_rows) return -1;
    if (res_rows != a_rows) return -1;
    if (res_cols != b_cols) return -1;

    for (size_t i = 0; i < a_rows; i++) {
        for (size_t j = 0; j < b_cols; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < a_cols; k++) {
                sum += *matrix_cptr(a, i, k) * *matrix_cptr(b, k, j);
            }
            *matrix_ptr(res, i, j) = sum;
        }
    }

    return 0;
}

// Поиск строки с максимальным элементом в столбе
static size_t find_pivot_row(const matrix *A, size_t col, size_t start) {
    size_t pivot = start;
    double max_val = fabs(*matrix_cptr(A, start, col));

    size_t rows = matrix_rows(A);
    for (size_t i = start + 1; i < rows; i++) {
        double val = fabs(*matrix_cptr(A, i, col));
        if (val > max_val) {
            max_val = val;
            pivot = i;
        }
    }

    return pivot;
}

// Перестановка двух строк матрицы
static void swap_rows(matrix *m, size_t r1, size_t r2) {
    if (r1 == r2) return;

    size_t cols = matrix_cols(m);
    for (size_t j = 0; j < cols; j++) {
        double tmp = *matrix_ptr(m, r1, j);
        *matrix_ptr(m, r1, j) = *matrix_ptr(m, r2, j);
        *matrix_ptr(m, r2, j) = tmp;
    }
}


matrix *matrix_solve_gauss(const matrix *A, const matrix *B) {

    if (A == NULL || B == NULL) {
        fprintf(stderr, "Ошибка: NULL указатель\n");
        return NULL;
    }

    size_t n = matrix_rows(A);


    if (n != matrix_cols(A)) {
        fprintf(stderr, "Ошибка: матрица A должна быть квадратной\n");
        return NULL;
    }


    if (n != matrix_rows(B)) {
        fprintf(stderr, "Ошибка: разное количество строк в A и B\n");
        return NULL;
    }

    size_t m = matrix_cols(B);


    matrix *A_copy = matrix_copy(A);
    matrix *X = matrix_copy(B);

    if (A_copy == NULL || X == NULL) {
        matrix_free(A_copy);
        matrix_free(X);
        return NULL;
    }

    for (size_t k = 0; k < n; k++) {

        size_t pivot = find_pivot_row(A_copy, k, k);
        double pivot_val = *matrix_cptr(A_copy, pivot, k);

        if (fabs(pivot_val) < 1e-15) {
            fprintf(stderr, "Ошибка: матрица вырождена\n");
            matrix_free(A_copy);
            matrix_free(X);
            return NULL;
        }

        if (pivot != k) {
            swap_rows(A_copy, k, pivot);
            swap_rows(X, k, pivot);
            pivot_val = *matrix_cptr(A_copy, k, k);
        }


        for (size_t j = k; j < n; j++) {
            *matrix_ptr(A_copy, k, j) /= pivot_val;
        }
        for (size_t j = 0; j < m; j++) {
            *matrix_ptr(X, k, j) /= pivot_val;
        }


        for (size_t i = k + 1; i < n; i++) {
            double factor = *matrix_cptr(A_copy, i, k);
            if (fabs(factor) < 1e-15) continue;

            for (size_t j = k; j < n; j++) {
                *matrix_ptr(A_copy, i, j) -= factor * *matrix_cptr(A_copy, k, j);
            }
            for (size_t j = 0; j < m; j++) {
                *matrix_ptr(X, i, j) -= factor * *matrix_cptr(X, k, j);
            }
        }
    }


    for (size_t k = n; k-- > 0;) {
        for (size_t i = 0; i < k; i++) {
            double factor = *matrix_cptr(A_copy, i, k);
            if (fabs(factor) < 1e-15) continue;

            for (size_t j = 0; j < m; j++) {
                *matrix_ptr(X, i, j) -= factor * *matrix_cptr(X, k, j);
            }
        }
    }

    matrix_free(A_copy);
    return X;
}

double matrix_check_solution(const matrix *A, const matrix *B, const matrix *X) {

    if (A == NULL || B == NULL || X == NULL) {
        return -1.0;
    }

    if (matrix_rows(A) != matrix_rows(B) ||
        matrix_cols(A) != matrix_rows(X) ||
        matrix_cols(B) != matrix_cols(X)) {
        return -1.0;
    }

    matrix *AX = matrix_alloc(matrix_cols(X), matrix_rows(A));
    if (AX == NULL) {
        return -1.0;
    }

    if (gauss_matrix_multiply(AX, A, X) != 0) {
        matrix_free(AX);
        return -1.0;
    }

    double max_diff = 0.0;
    size_t ax_rows = matrix_rows(AX);
    size_t ax_cols = matrix_cols(AX);

    for (size_t i = 0; i < ax_rows; i++) {
        for (size_t j = 0; j < ax_cols; j++) {
            double diff = fabs(*matrix_cptr(AX, i, j) - *matrix_cptr(B, i, j));
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
    }

    matrix_free(AX);
    return max_diff;
}
