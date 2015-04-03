#include "ranlib/ranlib.h"

#include <lapacke.h>
#include <lapacke_utils.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define $(i, j) ((i) * 3 + (j))
#define $$(i, j) ((i) * 3 + (j) + 9)

#define N_EQUATIONS 243
#define SHORTED_N_EQUATIONS 122
#define N_VARS 18

// Выделяет память для матрицы разсера w * l
inline lapack_complex_double * create_matrix(size_t w, size_t l) {
    return (lapack_complex_double *)calloc(w * l, sizeof(lapack_complex_double));
}


// Заполняет матрицу случайными комплексными числами
void fill_matrix(lapack_complex_double *matrix) {
    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j)  {
            // Генерируются числа в отрезке [-5.0,5.0]
            // (ranf() возврашает числа в отрезке [0,1])
            int c = 5;
            double a = (2.0 *ranf() - 1.0) * c;
            double b = (2.0 *ranf() - 1.0) * c;
            matrix[$(i, j)] = lapack_make_complex_double(a, b);
        }
    }
}

void fill_matrix_except_one_element(lapack_complex_double *matrix) {
    int i, j;

    int not_fill_i = rand() % 3;
    int not_fill_j = rand() % 3;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j)  {
            if (i == not_fill_i && j == not_fill_j) {
                int c = 5;
                double a = (2.0 *ranf() - 1.0) * c;
                double b = (2.0 *ranf() - 1.0) * c;
                matrix[$(i, j)] = lapack_make_complex_double(a, b);
            }
        }
    }
}


// Просто распечатывает матрицу
inline void  print_matrix(lapack_complex_double *matrix, int w, int l) {
    int i, j;

    for (i = 0; i < w; ++i) {
        for (j = 0; j < l; ++j) {
            printf("(%.3g, %.3g)  ", lapack_complex_double_real(matrix[i * l + j]),
                                     lapack_complex_double_imag(matrix[i * l + j]));
        }
        printf("\n");
    }
}


/**
* Если рассмотреть  матрицы X_tAX_t^{-1}, то получается, если не учитывать \epsilon
* (так как в итоговые уравнения все равно они войдут), то после необходимого перемножения
* мы получим только лишь две различные матрицы A и A' (в А' "крутятся" 2 и 3 строка, 2 и 3 столбец)
*/


// Является перестановкой индексов, которая оставляет на месте индекс 0, но меняет местам 1 и 2
inline int map_index_x(int i) {
    switch (i) {
        case 0: return 0;
        case 1: return 2;
        case 2: return 1;
    }
    return i;
}


/**
* Если рассмотреть  матрицы Y_tAY_t^{-1}, то получается, если не учитывать \epsilon
* (так как в итоговые уравнения все равно они войдут), то после необходимого перемножения
* мы получим только лишь две различные матрицы A и A' (в А' "крутятся" 1 и 3 строка, 1 и 3 столбец)
*/


// Является перестановкой индексов, которая оставляет на месте индекс 1, но меняет местам 0 и 2
inline int map_index_y(int i) {
    switch (i) {
        case 0: return 2;
        case 1: return 1;
        case 2: return 0;
    }
    return i;
}


/**
* Если рассмотреть  матрицы Z_tAZ_t^{-1}, то получается, если не учитывать \epsilon
* (так как в итоговые уравнения все равно они войдут), то после необходимого перемножения
* мы получим только лишь две различные матрицы A' и A'' (в А' "крутятся" 1 и 3 строка, 1 и 3 столбец, а
* в А'' происходит циклический сдвиг по столбцам и строкам).
*
*/

// Является перестановкой индексов, которая оставляет на месте индекс 0, но меняет местам 1 и 2
inline int map_index_z1(int i) {
    return map_index_x(i);
}


// Является циклическим сдвигом, который совершают матрицы Z_t
inline int map_index_z2(int i) {
    switch (i) {
        case 0: return 2;
        case 1: return 0;
        case 2: return 1;
    }
    return i;
}


// Генерирует матрицу системы при свободныз переменных из матриц C и N
void fill_row_cn(lapack_complex_double *row,
                 lapack_complex_double *A,
                 lapack_complex_double *B,
                 lapack_complex_double *K,
                 lapack_complex_double *M,
                 int i, int j, int k, int l, int r, int s) {

    int row_size = ((3 * 3) * 2);
    memset(row, 0, sizeof(lapack_complex_double) * row_size);

    int rcx, scx, icx, jcx, kcx, lcx;

    rcx = map_index_x(r);
    scx = map_index_x(s);
    icx = map_index_x(i);
    jcx = map_index_x(j);
    kcx = map_index_x(k);
    lcx = map_index_x(l);

    row[$(r, s)] += A[$(i, j)] * B[$(k, l)];
    row[$(rcx, scx)] += A[$(icx, jcx)] * B[$(kcx, lcx)];

    int rny, sny, iny, jny, kny, lny;

    rny = map_index_y(r);
    sny = map_index_y(s);
    iny = map_index_y(i);
    jny = map_index_y(j);
    kny = map_index_y(k);
    lny = map_index_y(l);

    row[$$(r, s)] += K[$(i, j)] * M[$(k, l)];
    row[$$(rny, sny)] += K[$(iny, jny)] * M[$(kny, lny)];

    int rnz, snz, inz, jnz, knz, lnz;

    rnz = map_index_z2(r);
    snz = map_index_z2(s);
    inz = map_index_z2(i);
    jnz = map_index_z2(j);
    knz = map_index_z2(k);
    lnz = map_index_z2(l);

    r = map_index_z1(r);
    s = map_index_z1(s);
    i = map_index_z1(i);
    j = map_index_z1(j);
    k = map_index_z1(k);
    l = map_index_z1(l);

    row[$$(r, s)] += K[$(i, j)] * M[$(k, l)];
    row[$$(rnz, snz)] += K[$(inz, jnz)] * M[$(knz, lnz)];
}


// Генерирует матрицу системы при свободныз переменных из матриц B и M
void fill_row_bm(lapack_complex_double *row,
        lapack_complex_double *A,
        lapack_complex_double *C,
        lapack_complex_double *K,
        lapack_complex_double *N,
        int i, int j, int k, int l, int r, int s) {

    int row_size = ((3 * 3) * 2);
    memset(row, 0, sizeof(lapack_complex_double) * row_size);

    int rbx, sbx, ibx, jbx, kbx, lbx;

    rbx = map_index_x(r);
    sbx = map_index_x(s);
    ibx = map_index_x(i);
    jbx = map_index_x(j);
    kbx = map_index_x(k);
    lbx = map_index_x(l);

    row[$(k, l)] += A[$(i, j)] * C[$(r, s)];
    row[$(kbx, lbx)] += A[$(ibx, jbx)] * C[$(rbx, sbx)];

    int rmy, smy, imy, jmy, kmy, lmy;

    rmy = map_index_y(r);
    smy = map_index_y(s);
    imy = map_index_y(i);
    jmy = map_index_y(j);
    kmy = map_index_y(k);
    lmy = map_index_y(l);

    row[$$(k, l)] += K[$(i, j)] * N[$(r, s)];
    row[$$(kmy, lmy)] += K[$(imy, jmy)] * N[$(rmy, smy)];

    int rmz, smz, imz, jmz, kmz, lmz;

    rmz = map_index_z2(r);
    smz = map_index_z2(s);
    imz = map_index_z2(i);
    jmz = map_index_z2(j);
    kmz = map_index_z2(k);
    lmz = map_index_z2(l);

    r = map_index_z1(r);
    s = map_index_z1(s);
    i = map_index_z1(i);
    j = map_index_z1(j);
    k = map_index_z1(k);
    l = map_index_z1(l);

    row[$$(k, l)] += K[$(i, j)] * N[$(r, s)];
    row[$$(kmz, lmz)] += K[$(imz, jmz)] * N[$(rmz, smz)];
}


// Генерирует матрицу системы при свободныз переменных из матриц A и K
void fill_row_ak(lapack_complex_double *row,
        lapack_complex_double *B,
        lapack_complex_double *C,
        lapack_complex_double *M,
        lapack_complex_double *N,
        int i, int j, int k, int l, int r, int s) {

    int row_size = ((3 * 3) * 2);
    memset(row, 0, sizeof(lapack_complex_double) * row_size);

    int rax, sax, iax, jax, kax, lax;

    rax = map_index_x(r);
    sax = map_index_x(s);
    iax = map_index_x(i);
    jax = map_index_x(j);
    kax = map_index_x(k);
    lax = map_index_x(l);

    row[$(i, j)] += B[$(k, l)] * C[$(r, s)];
    row[$(iax, jax)] += B[$(kax, lax)] * C[$(rax, sax)];

    int rky, sky, iky, jky, kky, lky;

    rky = map_index_y(r);
    sky = map_index_y(s);
    iky = map_index_y(i);
    jky = map_index_y(j);
    kky = map_index_y(k);
    lky = map_index_y(l);

    row[$$(i, j)] += M[$(k, l)] * N[$(r, s)];
    row[$$(iky, jky)] += M[$(kky, lky)] * N[$(rky, sky)];

    int rkz, skz, ikz, jkz, kkz, lkz;

    rkz = map_index_z2(r);
    skz = map_index_z2(s);
    ikz = map_index_z2(i);
    jkz = map_index_z2(j);
    kkz = map_index_z2(k);
    lkz = map_index_z2(l);

    r = map_index_z1(r);
    s = map_index_z1(s);
    i = map_index_z1(i);
    j = map_index_z1(j);
    k = map_index_z1(k);
    l = map_index_z1(l);

    row[$$(i, j)] += M[$(k, l)] * N[$(r, s)];
    row[$$(ikz, jkz)] += M[$(kkz, lkz)] * N[$(rkz, skz)];
}

// Читает из файла те наборы индексов, которые остаются в итоговой системе
void read_sets(int **sets) {
    int i;
    FILE *f;

    if((f = fopen("filter_sets.txt", "r")) < 0) {
        perror("File reading error:");
        exit(1);
    }

    for (i = 0; i < SHORTED_N_EQUATIONS; ++i) {
        fscanf(f, "%d %d %d %d %d %d", &sets[i][0], &sets[i][1], &sets[i][2], &sets[i][3], &sets[i][4], &sets[i][5]);
    }

    for (i = 0; i < SHORTED_N_EQUATIONS; ++i) {
        sets[i][0]--;
        sets[i][1]--;
        sets[i][2]--;
        sets[i][3]--;
        sets[i][4]--;
        sets[i][5]--;
    }

    fclose(f);
}

void free_matrix(lapack_complex_double *matrix) {
    free(matrix);
}


void free_sets(int **sets) {
    int i;

    for(i = 0; i < SHORTED_N_EQUATIONS; i++) {
        free(sets[i]);
    }

    free(sets);
}


// Создает вектор столбец правой части уравнения
lapack_complex_double create_lhs_vector(int i, int j, int k, int l, int r, int s) {
    if (j == k && l == r && s == i) {
        if (i == j && k == l && r == s) {
            return lapack_make_complex_double(0, 0);
        } else {
            return lapack_make_complex_double(1.0 / 3.0, 0);
        }
    } else {
        if (i == j && k == l && r == s) {
            return lapack_make_complex_double(- 1.0 / 3.0, 0);
        } else {
            return lapack_make_complex_double(0, 0);
        }
    }
}


// Вычисляется невязка
double compute_residual(lapack_complex_double *main_vector) {
    double residual = 0.0;
    int i;
    for (i = N_VARS; i < SHORTED_N_EQUATIONS; ++i) {
        residual += pow(lapack_complex_float_real(main_vector[i]), 2) +
                pow(lapack_complex_float_imag(main_vector[i]), 2);
    }

    return residual;
}


void fill_matrix_new(lapack_complex_double *dest, lapack_complex_double *src) {
    memcpy(dest, src, sizeof(lapack_complex_double) * 9);
}


double get_max_abs_of_matrix(lapack_complex_double *matrix) {
    int i, j;

    double max_abs = -1;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            double elem_max = sqrt(pow(lapack_complex_double_real(matrix[$(i,j)]), 2.0) +
                                   pow(lapack_complex_double_imag(matrix[$(i,j)]), 2.0));
            if (elem_max > max_abs) {
                max_abs = elem_max;
            }
        }
    }

    return max_abs;
}


void scalar_mult(lapack_complex_double *matrix, double scalar) {
    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            matrix[$(i, j)] *= scalar;
        }
    }
}


int main(void) {
    int i;

    lapack_complex_double *A, *B, *C, *K, *M, *N;
    lapack_complex_double *MAIN_MATRIX, *MAIN_VECTOR;
    int **sets = (int **)calloc(SHORTED_N_EQUATIONS, sizeof(int *));

    for(i = 0; i < SHORTED_N_EQUATIONS; i++) {
        sets[i] = (int *)calloc(6, sizeof(int));
    }

    read_sets(sets);
    srand(time(0));
    setall(time(0), time(0));

    A = create_matrix(3, 3);
    B = create_matrix(3, 3);
    C = create_matrix(3, 3);
    K = create_matrix(3, 3);
    M = create_matrix(3, 3);
    N = create_matrix(3, 3);

    fill_matrix(A);
    fill_matrix(B);
    fill_matrix(C);
    fill_matrix(K);
    fill_matrix(M);
    fill_matrix(N);

    MAIN_MATRIX = (lapack_complex_double *)calloc(SHORTED_N_EQUATIONS * 18, sizeof(lapack_complex_double));
    MAIN_VECTOR = (lapack_complex_double *)calloc(SHORTED_N_EQUATIONS, sizeof(lapack_complex_double));

    void (*generate_func[3])(lapack_complex_double *,
                          lapack_complex_double *,
                          lapack_complex_double *,
                          lapack_complex_double *,
                          lapack_complex_double *,
                          int, int, int, int, int, int);

    generate_func[0] = fill_row_cn;
    generate_func[1] = fill_row_ak;
    generate_func[2] = fill_row_bm;

    int iterate = 0;
    int main_count = 200;
    double prev_residual = 0.0;
    double cur_residual = 0.0;

    double max_abs_a, max_abs_b, max_abs_c;
    double max_abs_k, max_abs_m, max_abs_n;

    lapack_complex_double *matrix1, *matrix2, *matrix3, *matrix4;

    // MAIN LOOP START
    do {
        int matrix_index = iterate % 3;

        if (matrix_index == 0) {
            matrix1 = A;
            matrix2 = B;
            matrix3 = K;
            matrix4 = M;
        } else if (matrix_index == 1) {
            matrix1 = B;
            matrix2 = C;
            matrix3 = M;
            matrix4 = N;
        } else {
            matrix1 = A;
            matrix2 = C;
            matrix3 = K;
            matrix4 = N;
        }

        for (i = 0; i < SHORTED_N_EQUATIONS; ++i) {
            generate_func[matrix_index](&MAIN_MATRIX[i * 18], matrix1, matrix2, matrix3, matrix4,
                    sets[i][0], sets[i][1], sets[i][2],
                    sets[i][3], sets[i][4], sets[i][5]);
            MAIN_VECTOR[i] = create_lhs_vector(sets[i][0], sets[i][1], sets[i][2],
                    sets[i][3], sets[i][4], sets[i][5]);
        }

        int info;
        info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', SHORTED_N_EQUATIONS, N_VARS, 1, MAIN_MATRIX, N_VARS, MAIN_VECTOR, 1);

        if (info > 0) {
            perror("The solution could not be computed (the matrix have not the full rank)");
            exit(1);
        }

        if (matrix_index == 0) {
            fill_matrix_new(C, &MAIN_VECTOR[0]);
            fill_matrix_new(N, &MAIN_VECTOR[9]);

            if (iterate == 0) {

                max_abs_a = get_max_abs_of_matrix(A);
                max_abs_b = get_max_abs_of_matrix(B);
                max_abs_c = get_max_abs_of_matrix(C);

                double k = pow(max_abs_a * max_abs_b * max_abs_c, 1 / 3.);

                scalar_mult(A, k / max_abs_a);
                scalar_mult(B, k / max_abs_b);
                scalar_mult(C, k / max_abs_c);

                max_abs_k = get_max_abs_of_matrix(K);
                max_abs_m = get_max_abs_of_matrix(M);
                max_abs_n = get_max_abs_of_matrix(N);

                k = pow(max_abs_k * max_abs_m * max_abs_n, 1 / 3.);

                scalar_mult(K, k / max_abs_k);
                scalar_mult(M, k / max_abs_m);
                scalar_mult(N, k / max_abs_n);
            }

        } else if (matrix_index == 1) {
            fill_matrix_new(A, &MAIN_VECTOR[0]);
            fill_matrix_new(K, &MAIN_VECTOR[9]);
        } else {
            fill_matrix_new(B, &MAIN_VECTOR[0]);
            fill_matrix_new(M, &MAIN_VECTOR[9]);
        }

        prev_residual = cur_residual;
        cur_residual = compute_residual(MAIN_VECTOR);

        iterate++;

        if (fabs(cur_residual - prev_residual) < 0.00000005) {
            iterate = 0;

            fill_matrix(A);
            fill_matrix(B);
            fill_matrix(K);
            fill_matrix(M);

            printf("\n\n################# NEW ITERATION #################\n");

            print_matrix(A, 3, 3);
            printf("\n");
            print_matrix(B, 3, 3);
            printf("\n");
            print_matrix(C, 3, 3);
            printf("\n");
            print_matrix(K, 3, 3);
            printf("\n");
            print_matrix(M, 3, 3);
            printf("\n");
            print_matrix(N, 3, 3);

            printf("RESIDUAL = %f\n", cur_residual);
            prev_residual = cur_residual = 0.0;
            main_count--;
        }
    } while (main_count != 0);
    // MAIN LOOP END

    free_matrix(A);
    free_matrix(B);
    free_matrix(C);
    free_matrix(K);
    free_matrix(M);
    free_matrix(N);
    free_matrix(MAIN_MATRIX);
    free_sets(sets);

    return 0;
}