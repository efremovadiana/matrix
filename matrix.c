#include "matrix.h"
#include <stdlib.h>
#include <string.h>

struct matrix {
    double *data;
    size_t w;
    size_t h;
};

matrix *matrix_alloc(size_t w, size_t h) {
    if (w == 0 || h == 0) return NULL;

    matrix *m = malloc(sizeof(matrix));
    if (!m) return NULL;

    m->data = calloc(w * h, sizeof(double));
    if (!m->data) {
        free(m);
        return NULL;
    }

    m->w = w;
    m->h = h;
    return m;
}

matrix *matrix_copy(const matrix *m) {
    if (!m) return NULL;

    matrix *new = matrix_alloc(m->w, m->h);
    if (!new) return NULL;

    memcpy(new->data, m->data, m->w * m->h * sizeof(double));
    return new;
}

void matrix_free(matrix *m) {
    if (m) {
        free(m->data);
        free(m);
    }
}

double *matrix_ptr(matrix *m, size_t i, size_t j) {
    return m->data + m->w * i + j;
}

const double *matrix_cptr(const matrix *m, size_t i, size_t j) {
    return m->data + m->w * i + j;
}

size_t matrix_rows(const matrix *m) {
    return m ? m->h : 0;
}

size_t matrix_cols(const matrix *m) {
    return m ? m->w : 0;
}
