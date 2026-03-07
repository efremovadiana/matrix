#include "matrix.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

matrix *matrix_assign(matrix *m1, const matrix *m2) {
    if (!m1 || !m2) return NULL;

    size_t rows = matrix_rows(m1);
    size_t cols = matrix_cols(m1);

    if (rows != matrix_rows(m2) || cols != matrix_cols(m2)) {
        return NULL;
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            *matrix_ptr(m1, i, j) = *matrix_cptr(m2, i, j);
        }
    }

    return m1;
}

int matrix_add(matrix *m1, const matrix *m2) {
    if (!m1 || !m2) return -1;

    size_t rows = matrix_rows(m1);
    size_t cols = matrix_cols(m1);

    if (rows != matrix_rows(m2) || cols != matrix_cols(m2)) {
        return -1;
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            *matrix_ptr(m1, i, j) += *matrix_cptr(m2, i, j);
        }
    }

    return 0;
}

int matrix_sub(matrix *m1, const matrix *m2) {
    if (!m1 || !m2) return -1;

    size_t rows = matrix_rows(m1);
    size_t cols = matrix_cols(m1);

    if (rows != matrix_rows(m2) || cols != matrix_cols(m2)) {
        return -1;
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            *matrix_ptr(m1, i, j) -= *matrix_cptr(m2, i, j);
        }
    }

    return 0;
}

void matrix_smul(matrix *m, double d) {
    if (!m) return;

    size_t rows = matrix_rows(m);
    size_t cols = matrix_cols(m);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            *matrix_ptr(m, i, j) *= d;
        }
    }
}

void matrix_sdiv(matrix *m, double d) {
    if (!m || d == 0.0) return;
    matrix_smul(m, 1.0 / d);
}

matrix *matrix_alloc_id(size_t n) {
    if (n == 0) return NULL;

    matrix *m = matrix_alloc(n, n);
    if (!m) return NULL;

    for (size_t i = 0; i < n; i++) {
        *matrix_ptr(m, i, i) = 1.0;
    }

    return m;
}

double matrix_norm(const matrix *m) {
    if (!m) return -1.0;

    size_t rows = matrix_rows(m);
    size_t cols = matrix_cols(m);
    double max = 0.0;

    for (size_t i = 0; i < rows; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < cols; j++) {
            sum += fabs(*matrix_cptr(m, i, j));
        }
        if (sum > max) max = sum;
    }
    return max;
}

int matrix_mul(matrix *res, const matrix *a, const matrix *b) {
    if (!res || !a || !b) return -1;

    size_t a_rows = matrix_rows(a);
    size_t a_cols = matrix_cols(a);
    size_t b_rows = matrix_rows(b);
    size_t b_cols = matrix_cols(b);

    if (a_cols != b_rows ||
        matrix_rows(res) != a_rows ||
        matrix_cols(res) != b_cols) {
        return -1;
    }

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

int matrix_mul_assign(matrix *a, const matrix *b) {
    if (!a || !b) return -1;

    size_t a_rows = matrix_rows(a);
    size_t a_cols = matrix_cols(a);
    size_t b_rows = matrix_rows(b);
    size_t b_cols = matrix_cols(b);

    if (a_cols != b_rows) return -1;

    matrix *temp = matrix_alloc(b_cols, a_rows);
    if (!temp) return -1;

    if (matrix_mul(temp, a, b) != 0) {
        matrix_free(temp);
        return -1;
    }

    matrix_assign(a, temp);
    matrix_free(temp);
    return 0;
}

int matrix_transpose(matrix *m) {
    if (!m) return -1;

    size_t n = matrix_rows(m);
    if (n != matrix_cols(m)) return -1;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            double tmp = *matrix_ptr(m, i, j);
            *matrix_ptr(m, i, j) = *matrix_ptr(m, j, i);
            *matrix_ptr(m, j, i) = tmp;
        }
    }
    return 0;
}

matrix *matrix_transpose_new(const matrix *m) {
    if (!m) return NULL;

    size_t rows = matrix_rows(m);
    size_t cols = matrix_cols(m);

    matrix *res = matrix_alloc(rows, cols);
    if (!res) return NULL;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            *matrix_ptr(res, j, i) = *matrix_cptr(m, i, j);
        }
    }

    return res;
}

void matrix_swap_rows(matrix *m, size_t r1, size_t r2) {
    if (!m) return;

    size_t rows = matrix_rows(m);
    size_t cols = matrix_cols(m);

    if (r1 >= rows || r2 >= rows || r1 == r2) return;

    for (size_t j = 0; j < cols; j++) {
        double tmp = *matrix_ptr(m, r1, j);
        *matrix_ptr(m, r1, j) = *matrix_ptr(m, r2, j);
        *matrix_ptr(m, r2, j) = tmp;
    }
}

void matrix_swap_cols(matrix *m, size_t c1, size_t c2) {
    if (!m) return;

    size_t rows = matrix_rows(m);
    size_t cols = matrix_cols(m);

    if (c1 >= cols || c2 >= cols || c1 == c2) return;

    for (size_t i = 0; i < rows; i++) {
        double tmp = *matrix_ptr(m, i, c1);
        *matrix_ptr(m, i, c1) = *matrix_ptr(m, i, c2);
        *matrix_ptr(m, i, c2) = tmp;
    }
}

void matrix_row_mul(matrix *m, size_t row, double d) {
    if (!m) return;

    size_t rows = matrix_rows(m);
    size_t cols = matrix_cols(m);

    if (row >= rows) return;

    for (size_t j = 0; j < cols; j++) {
        *matrix_ptr(m, row, j) *= d;
    }
}

void matrix_row_add(matrix *m, size_t to, size_t from, double coeff) {
    if (!m) return;

    size_t rows = matrix_rows(m);
    size_t cols = matrix_cols(m);

    if (to >= rows || from >= rows) return;

    for (size_t j = 0; j < cols; j++) {
        *matrix_ptr(m, to, j) += coeff * *matrix_cptr(m, from, j);
    }
}

void matrix_print(const matrix *m, const char *format) {
    if (!m) {
        printf("NULL\n");
        return;
    }

    size_t rows = matrix_rows(m);
    size_t cols = matrix_cols(m);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            printf(format, *matrix_cptr(m, i, j));
        }
        printf("\n");
    }
}
