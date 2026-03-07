#include "matrix.h"
#include <math.h>

matrix *matrix_exp(const matrix *A, double eps) {
    if (!A || eps <= 0.0) return NULL;

    size_t n = matrix_rows(A);
    if (n != matrix_cols(A)) return NULL;

    matrix *result = matrix_alloc_id(n);
    if (!result) return NULL;

    matrix *term = matrix_copy(A);
    if (!term) {
        matrix_free(result);
        return NULL;
    }

    matrix_add(result, term);

    matrix *temp = NULL;
    unsigned int k = 2;

    while (1) {
        temp = matrix_alloc(n, n);
        if (!temp) break;

        if (matrix_mul(temp, term, A) != 0) {
            matrix_free(temp);
            break;
        }

        matrix_assign(term, temp);
        matrix_free(temp);

        matrix_sdiv(term, (double)k);

        if (matrix_norm(term) < eps) break;

        matrix_add(result, term);
        k++;
    }

    matrix_free(term);
    return result;
}
