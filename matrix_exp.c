#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

matrix *matrix_exp(const matrix *A, double eps) {

    if (A == NULL) {
        return NULL;
    }

    size_t n = matrix_rows(A);

    if (n != matrix_cols(A)) {
        fprintf(stderr, "Матрица должна быть квадратной\n");
        return NULL;
    }

    matrix *result = matrix_alloc_id(n);
    if (result == NULL) {
        return NULL;
    }

    matrix *term = matrix_copy(A);
    if (term == NULL) {
        matrix_free(result);
        return NULL;
    }


    matrix_add(result, term);

    matrix *temp = NULL;
    unsigned int k = 2;

    while (1) {

        temp = matrix_alloc(n, n);
        if (temp == NULL) {
            break;
        }

        if (matrix_multiply(temp, term, A) != 0) {
            matrix_free(temp);
            break;
        }

        matrix_assign(term, temp);
        matrix_free(temp);

        matrix_sdiv(term, (double)k);

        double norm = matrix_norm(term);
        if (norm < eps) {
            break;
        }

        matrix_add(result, term);

        k++;
    }

    matrix_free(term);

    return result;
}
