#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

struct matrix {
    double *data;
    size_t w;
    size_t h;
};

matrix *matrix_alloc(size_t w, size_t h){


        if (w == 0 || h == 0){
            return NULL;
        }

    matrix *m = (matrix *)malloc(sizeof(matrix));

    if (m == NULL){
        return NULL;
    }

    m -> data = (double *)calloc(w * h, sizeof(double));

    if ( m -> data == NULL){
        free(m);
        return NULL;
    }

    m -> w = w;
    m -> h = h;

    return m;

}

matrix *matrix_copy(const matrix *m){

    if (m == NULL){
        return NULL;
    }

    matrix *new_m = matrix_alloc(m -> w, m -> h);
    if (new_m == NULL){
        return NULL;
    }

    size_t total_bytes = m->w * m->h * sizeof(double);
    memcpy(new_m -> data, m -> data, total_bytes);

    return new_m;
}

void matrix_free(matrix *m){

    if (m != NULL){

        free(m -> data);
        free(m);
    }
}

double *matrix_ptr(matrix *m, size_t i, size_t j){

    return m -> data + m -> w * i + j;
}

const double *matrix_cptr(const matrix *m, size_t i, size_t j){

    return m -> data + m -> w * i + j;
}

size_t matrix_rows(const matrix *m) {
    return (m != NULL) ? m->h : 0;
}

size_t matrix_cols(const matrix *m) {
    return (m != NULL) ? m->w : 0;
}

matrix *matrix_alloc_id(size_t n) {

    if (n == 0) {
        return NULL;
    }

    matrix *m = matrix_alloc(n, n);
    if (m == NULL) {
        return NULL;
    }

     for (size_t i = 0; i < n; i++) {

        *matrix_ptr(m, i, i) = 1.0;
    }

    return m;
}

matrix *matrix_assign(matrix *m1, const matrix *m2) {

    if (m1 == NULL || m2 == NULL) {
        return NULL;
    }

    if (m1->w != m2->w || m1->h != m2->h) {
        return NULL;
    }

    size_t total_bytes = m1->w * m1->h * sizeof(double);
    memcpy(m1->data, m2->data, total_bytes);

    return m1;
}

int matrix_add(matrix *m1, const matrix *m2) {

    if (m1 == NULL || m2 == NULL) {
        return -1;
    }

    if (m1->w != m2->w || m1->h != m2->h) {
        return -1;
    }


    size_t total = m1->w * m1->h;
    for (size_t k = 0; k < total; k++) {
        m1->data[k] += m2->data[k];
    }

    return 0;
}

void matrix_smul(matrix *m, double d) {
    if (m == NULL) {
        return;
    }

    size_t total = m->w * m->h;
    for (size_t k = 0; k < total; k++) {
        m->data[k] *= d;
    }
}

void matrix_sdiv(matrix *m, double d) {

    if (m == NULL || d == 0.0) {
        return;
    }

    matrix_smul(m, 1.0 / d);
}


double matrix_norm(const matrix *m) {
    if (m == NULL) return -1.0;

    double max_sum = 0.0;

    for (size_t i = 0; i < m->h; i++) {
        double row_sum = 0.0;

        for (size_t j = 0; j < m->w; j++) {
            row_sum += fabs(*matrix_cptr(m, i, j));
        }

        if (row_sum > max_sum) {
            max_sum = row_sum;
        }
    }

    return max_sum;
}

int matrix_multiply(matrix *res, const matrix *a, const matrix *b) {

    if (res == NULL || a == NULL || b == NULL) {
        return -1;
    }

    if (a->w != b->h) {
        return -1;
    }
    if (res->h != a->h || res->w != b->w) {
        return -1;
    }

    for (size_t i = 0; i < a->h; i++) {
        for (size_t j = 0; j < b->w; j++) {
            double sum = 0.0;

            for (size_t k = 0; k < a->w; k++) {
                sum += *matrix_cptr(a, i, k) * *matrix_cptr(b, k, j);
            }

            *matrix_ptr(res, i, j) = sum;
        }
    }

    return 0;
}

void matrix_print(const matrix *m, const char *format) {
    if (m == NULL) {
        printf("NULL матрица\n");
        return;
    }

    for (size_t i = 0; i < m->h; i++) {
        for (size_t j = 0; j < m->w; j++) {
            printf(format, *matrix_cptr(m, i, j));
        }
        printf("\n");
    }
}


