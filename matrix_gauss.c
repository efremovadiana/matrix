#include "matrix.h"
#include <math.h>

/* Поиск главного элемента, А - матрица, col - индекс столбца,
 *start - начальна строка, max_val - указатель для сохранения максимального значения
 */
static size_t find_pivot(const matrix *A, size_t col, size_t start, double *max_val) {
    size_t pivot = start;
    *max_val = fabs(*matrix_cptr(A, start, col));

    double scale = 0.0;
    for (size_t j = 0; j < matrix_cols(A); j++) {
        scale += fabs(*matrix_cptr(A, start, j));
    }
    if (scale == 0.0) scale = 1.0;

    double scaled_max = *max_val / scale;

    size_t rows = matrix_rows(A);
    for (size_t i = start + 1; i < rows; i++) {
        double val = fabs(*matrix_cptr(A, i, col));

        scale = 0.0;
        for (size_t j = 0; j < matrix_cols(A); j++) {
            scale += fabs(*matrix_cptr(A, i, j));
        }
        if (scale == 0.0) scale = 1.0;

        double scaled = val / scale;
        if (scaled > scaled_max) {
            scaled_max = scaled;
            *max_val = val;
            pivot = i;
        }
    }

    return pivot;
}

matrix *matrix_solve_gauss(const matrix *A, const matrix *B, double tol) {

    if (!A || !B) return NULL;

    size_t n = matrix_rows(A);
    size_t rhs_count = matrix_cols(B);

    if (n != matrix_cols(A) || n != matrix_rows(B)) {
        return NULL;
    }

    if (tol <= 0.0) {
        double normA = matrix_norm(A);
        tol = (normA > 0.0) ? normA * 1e-15 : 1e-15;
    }

    matrix *A_copy = matrix_copy(A);
    matrix *X = matrix_copy(B);

    if (!A_copy || !X) {
        matrix_free(A_copy);
        matrix_free(X);
        return NULL;
    }

    for (size_t k = 0; k < n; k++) {
        double pivot_val;
        size_t pivot = find_pivot(A_copy, k, k, &pivot_val);

        if (pivot_val <= tol) {
            matrix_free(A_copy);
            matrix_free(X);
            return NULL;
        }

        if (pivot != k) {
            matrix_swap_rows(A_copy, k, pivot);
            matrix_swap_rows(X, k, pivot);
        }

        double inv_pivot = 1.0 / pivot_val;
        matrix_row_mul(A_copy, k, inv_pivot);
        matrix_row_mul(X, k, inv_pivot);


        for (size_t i = k + 1; i < n; i++) {
            double factor = *matrix_cptr(A_copy, i, k);
            if (fabs(factor) <= tol) continue;

            matrix_row_add(A_copy, i, k, -factor);
            matrix_row_add(X, i, k, -factor);
        }
    }


    for (size_t k = n; k-- > 0;) {
        for (size_t i = 0; i < k; i++) {
            double factor = *matrix_cptr(A_copy, i, k);
            if (fabs(factor) <= tol) continue;


            for (size_t j = 0; j < rhs_count; j++) {
                *matrix_ptr(X, i, j) -= factor * *matrix_cptr(X, k, j);
            }
        }
    }

    matrix_free(A_copy);
    return X;
}


double matrix_check_solution(const matrix *A, const matrix *B, const matrix *X) {
    if (!A || !B || !X) return -1.0;

    size_t n = matrix_rows(A);
    size_t m = matrix_cols(B);

    if (matrix_cols(A) != matrix_rows(X) ||
        n != matrix_rows(B) ||
        m != matrix_cols(X)) {
        return -1.0;
    }

    matrix *AX = matrix_alloc(m, n);
    if (!AX) return -1.0;

    if (matrix_mul(AX, A, X) != 0) {
        matrix_free(AX);
        return -1.0;
    }

    double max_diff = 0.0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            double diff = fabs(*matrix_cptr(AX, i, j) - *matrix_cptr(B, i, j));
            if (diff > max_diff) max_diff = diff;
        }
    }

    matrix_free(AX);
    return max_diff;
}
