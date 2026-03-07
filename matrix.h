#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <stddef.h>

typedef struct matrix matrix;

/* Создаёт новую матрицу, w - ширина, h - высота
 * return: новая матрица или NULL (нужно освободить через matrix_free)
 */
matrix *matrix_alloc(size_t w, size_t h);

/*
 * Создаёт копию матрицы
 * m - исходная матрица
 * return: копия или NULL (нужно освободить через matrix_free)
 */
matrix *matrix_copy(const matrix *m);

/*
 * Освобождает память, занятую матрицей
 * m матрица
 */
void matrix_free(matrix *m);

/*
 * Доступ к элементу (без проверки границ)
 * m матрица, i строка, j столбец
 */
double *matrix_ptr(matrix *m, size_t i, size_t j);

/*
 * Константный доступ к элементу
 * m матрица, i строка, j столбец
 */
const double *matrix_cptr(const matrix *m, size_t i, size_t j);

/*
 * Получить количество строк
 * m матрица
 */
size_t matrix_rows(const matrix *m);

/*
 * Получить количество столбцов
 * m матрица
 */
size_t matrix_cols(const matrix *m);


/*
 * Создаёт единичную матрицу
 * n размер
 * return: новая единичная матрица (нужно освободить через matrix_free)
 */
matrix *matrix_alloc_id(size_t n);

/*
 * Копирует содержимое m2 в m1
 * m1 - целевая матрица, m2 - исходная матрица
 */
matrix *matrix_assign(matrix *m1, const matrix *m2);

/*
 * Поэлементное сложение m1 += m2
 * m1, m2 - матрицы
 */
int matrix_add(matrix *m1, const matrix *m2);

/*
 * Поэлементное вычитание m1 -= m2
 * m1, m2 - матрицы
 */
int matrix_sub(matrix *m1, const matrix *m2);

/*
 * Умножение всех элементов на число
 * m матрица, d множитель
 */
void matrix_smul(matrix *m, double d);

/*
 * Деление всех элементов на число
 * m матрица, d делитель
 */
void matrix_sdiv(matrix *m, double d);



/*
 * Норма матрицы (максимум суммы модулей по строкам)
 * m матрица
 */
double matrix_norm(const matrix *m);

/*
 * Умножение матриц res = a * b
 * res матрица для результата
 * a, b - матрицы
 */
int matrix_mul(matrix *res, const matrix *a, const matrix *b);

/*
 * Умножение с присваиванием a *= b
 * a левая матрица, b правая матрица
 */
int matrix_mul_assign(matrix *a, const matrix *b);

/*
 * Транспонирование квадратной матрицы на месте
 * m матрица
 */
int matrix_transpose(matrix *m);

/*
 * Создать транспонированную копию
 * m исходная матрица
 * return: новая матрица (нужно освободить через matrix_free)
 */
matrix *matrix_transpose_new(const matrix *m);

/*
 * Перестановка строк
 * m матрица, r1 первая строка, r2 вторая строка
 */
void matrix_swap_rows(matrix *m, size_t r1, size_t r2);

/*
 * Перестановка столбцов
 * m матрица, c1 первый столбец, c2 второй столбец
 */
void matrix_swap_cols(matrix *m, size_t c1, size_t c2);

/*
 * Умножение строки на число
 * m матрица, row индекс строки, d множитель
 */
void matrix_row_mul(matrix *m, size_t row, double d);

/*
 * Сложение строк: row_to += row_from * coeff
 * m матрица, to строка-приёмник, from строка-источник, coeff коэффициент
 */
void matrix_row_add(matrix *m, size_t to, size_t from, double coeff);



/*
 * Вычислить матричную экспоненту
 * A исходная матрица (квадратная), eps точность
 * return: новая матрица exp(A) (нужно освободить через matrix_free) или NULL при ошибке
 */
matrix *matrix_exp(const matrix *A, double eps);

/*
 * Решение СЛАУ методом Гаусса
 * A матрица коэффициентов (квадратная), B матрица правых частей, tol допуск для определения вырожденности (0 - автоматический)
 * return: матрица X решений (нужно освободить через matrix_free) или NULL при ошибке
 */
matrix *matrix_solve_gauss(const matrix *A, const matrix *B, double tol);

/*
 * Вычислить невязку ||A*X - B||
 * A матрица коэффициентов, B матрица правых частей, X найденное решение
 */
double matrix_check_solution(const matrix *A, const matrix *B, const matrix *X);



/*
 * Вывод матрицы на экран
 * m матрица, format формат вывода
 */
void matrix_print(const matrix *m, const char *format);


#endif // MATRIX_H_INCLUDED
