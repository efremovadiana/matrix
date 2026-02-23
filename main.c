#include "matrix.h"
#include "matrix_exp.h"
#include "matrix_gauss.h"
#include <stdio.h>
#include <math.h>
#include <locale.h>


void test_matrix_exp(void) {

    printf("Тест 1: Матричная экспонента\n");

    matrix *A = matrix_alloc_id(3);
    *matrix_ptr(A, 0, 0) = 1.0;
    *matrix_ptr(A, 1, 1) = 2.0;
    *matrix_ptr(A, 2, 2) = -1.0;

    printf("\nИсходная матрица A:\n");
    matrix_print(A, "%8.3f");

    double eps = 1e-10;
    printf("\nВычисление exp(A) с точностью %g...\n", eps);

    matrix *expA = matrix_exp(A, eps);

    if (expA != NULL) {
        printf("\nРезультат exp(A):\n");
        matrix_print(expA, "%12.6f");

        printf("\nОжидаемый результат:\n");
        printf("  %12.6f %12.6f %12.6f\n", exp(1.0), 0.0, 0.0);
        printf("  %12.6f %12.6f %12.6f\n", 0.0, exp(2.0), 0.0);
        printf("  %12.6f %12.6f %12.6f\n", 0.0, 0.0, exp(-1.0));

        matrix_free(expA);
    } else {
        printf("\nОШИБКА: не удалось вычислить экспоненту\n");
    }

    matrix_free(A);
    printf("\nТест 1 завершён.\n");
}


void test_gauss_simple(void) {

    printf("Тест 2: Метод Гаусса (одна правая часть)\n");

    matrix *A = matrix_alloc(3, 3);
    *matrix_ptr(A, 0, 0) = 2.0;  *matrix_ptr(A, 0, 1) = 1.0;  *matrix_ptr(A, 0, 2) = -1.0;
    *matrix_ptr(A, 1, 0) = 1.0;  *matrix_ptr(A, 1, 1) = 3.0;  *matrix_ptr(A, 1, 2) = 2.0;
    *matrix_ptr(A, 2, 0) = 1.0;  *matrix_ptr(A, 2, 1) = 0.0;  *matrix_ptr(A, 2, 2) = 0.0;

    matrix *B = matrix_alloc(1, 3);
    *matrix_ptr(B, 0, 0) = 8.0;
    *matrix_ptr(B, 1, 0) = 9.0;
    *matrix_ptr(B, 2, 0) = 3.0;

    printf("\nМатрица A:\n");
    matrix_print(A, "%6.2f");
    printf("\nВектор B:\n");
    matrix_print(B, "%6.2f");

    matrix *X = matrix_solve_gauss(A, B);

    if (X != NULL) {
        printf("\nРешение X:\n");
        matrix_print(X, "%10.6f");

        double norm = matrix_check_solution(A, B, X);
        printf("\nНевязка ||AX - B|| = %g\n", norm);

        printf("\nОжидаемое решение: x = 3, y = 2, z = 0\n");

        matrix_free(X);
    } else {
        printf("\nОшибка: не удалось решить систему\n");
    }

    matrix_free(A);
    matrix_free(B);
    printf("\nТест 2 завершён.\n");
}


void test_gauss_multi(void) {

    printf("Тест 3: Метод Гаусса (несколько правых частей)\n");

    matrix *A = matrix_alloc(2, 2);
    *matrix_ptr(A, 0, 0) = 1.0;  *matrix_ptr(A, 0, 1) = 2.0;
    *matrix_ptr(A, 1, 0) = 3.0;  *matrix_ptr(A, 1, 1) = 4.0;

    matrix *B = matrix_alloc(2, 2);
    *matrix_ptr(B, 0, 0) = 5.0;   *matrix_ptr(B, 0, 1) = 4.0;
    *matrix_ptr(B, 1, 0) = 11.0;  *matrix_ptr(B, 1, 1) = 10.0;

    printf("\nМатрица A:\n");
    matrix_print(A, "%6.2f");
    printf("\nМатрица B:\n");
    matrix_print(B, "%6.2f");

    matrix *X = matrix_solve_gauss(A, B);

    if (X != NULL) {
        printf("\nРешение X (по столбцам):\n");
        matrix_print(X, "%10.6f");

        double norm = matrix_check_solution(A, B, X);
        printf("\nНевязка ||AX - B|| = %g\n", norm);

        matrix_free(X);
    } else {
        printf("\nОшибка: не удалось решить систему\n");
    }

    matrix_free(A);
    matrix_free(B);
    printf("\nТест 3 завершён.\n");
}


void test_singular(void) {

    printf("Тест 4: Вырожденная матрица\n");

    matrix *A = matrix_alloc(2, 2);
    *matrix_ptr(A, 0, 0) = 1.0;  *matrix_ptr(A, 0, 1) = 2.0;
    *matrix_ptr(A, 1, 0) = 2.0;  *matrix_ptr(A, 1, 1) = 4.0;

    matrix *B = matrix_alloc(1, 2);
    *matrix_ptr(B, 0, 0) = 5.0;
    *matrix_ptr(B, 1, 0) = 10.0;

    printf("\nМатрица A (вырожденная):\n");
    matrix_print(A, "%6.2f");


    matrix *X = matrix_solve_gauss(A, B);

    if (X == NULL) {
        printf("\nОжидаемая ошибка: матрица вырождена - решение невозможно\n");
    } else {
        printf("\nОшибка: получено решение для вырожденной матрицы\n");
        matrix_free(X);
    }

    matrix_free(A);
    matrix_free(B);
    printf("\nТест 4 завершён.\n");
}


void test_memory_leaks(void) {

    printf("Тест 5: Проверка управления памятью\n");

    printf("Создание и освобождение нескольких матриц\n");

    matrix *m1 = matrix_alloc(3, 3);
    matrix *m2 = matrix_alloc(4, 4);
    matrix *m3 = matrix_copy(m1);

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            *matrix_ptr(m1, i, j) = i * 3 + j;
        }
    }

    printf("m1 (3x3):\n");
    matrix_print(m1, "%5.0f");

    matrix_free(m1);
    matrix_free(m2);
    matrix_free(m3);

    printf("\nВсе матрицы освобождены.\n");
    printf("\nТест 5 завершён.\n");
}

int main(void) {

    setlocale(LC_ALL, "Russian");

    printf("Решение двух задач:\n");
    printf("1. Матричная экспонента\n");
    printf("2. Метод Гаусса для СЛАУ\n");

    test_matrix_exp();
    test_gauss_simple();
    test_gauss_multi();
    test_singular();
    test_memory_leaks();

    printf("Все тесты завершены\n");
    printf("\nДля проверки утечек памяти:\n");
    printf("valgrind --leak-check=full ./matrix_tasks\n");

    return 0;
}
