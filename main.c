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
#define $$$(i, j) ((i) * 3 + (j) + 18)

#define N_EQUATIONS 243
#define SHORTED_N_EQUATIONS 122
#define N_VARS 27

#define I 0
#define J 1
#define R 2
#define S 3
#define K 4
#define L 5



inline lapack_complex_double * create_matrix(size_t w, size_t l) {
    return (lapack_complex_double *)calloc(w * l, sizeof(lapack_complex_double));
}

void fillin_matrix(lapack_complex_double *matrix) {
    int i, j;

    int c;
    double a, b;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j)  {
            c = rand() % 10 + 1;
            a = (2.0 *ranf() - 1.0) * c;
            b = (2.0 *ranf() - 1.0) * c;
            matrix[$(i, j)] = lapack_make_complex_double(a, b);
        }
    }
}

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


inline int map_index_x(int i) {
    switch (i) {
        case 0: return 0;
        case 1: return 2;
        case 2: return 1;
    }
    return i;
}


inline int map_index_y(int i) {
    switch (i) {
        case 0: return 2;
        case 1: return 1;
        case 2: return 0;
    }
    return i;
}


inline int map_index_z(int i) {
    switch (i) {
        case 0: return 1;
        case 1: return 0;
        case 2: return 2;
    }
    return i;
}


void fillin_row(lapack_complex_double *row,
                lapack_complex_double *A,
                lapack_complex_double *B,
                lapack_complex_double *C,
                lapack_complex_double *D,
                lapack_complex_double *E,
                lapack_complex_double *F,
                int indexes[]) {

    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

    int indexes_x[6], indexes_y[6], indexes_z[6];

    indexes_x[I] = map_index_x(indexes[I]);
    indexes_x[J] = map_index_x(indexes[J]);
    indexes_x[K] = map_index_x(indexes[K]);
    indexes_x[L] = map_index_x(indexes[L]);
    indexes_x[R] = map_index_x(indexes[R]);
    indexes_x[S] = map_index_x(indexes[S]);

    indexes_y[I] = map_index_y(indexes[I]);
    indexes_y[J] = map_index_y(indexes[J]);
    indexes_y[K] = map_index_y(indexes[K]);
    indexes_y[L] = map_index_y(indexes[L]);
    indexes_y[R] = map_index_y(indexes[R]);
    indexes_y[S] = map_index_y(indexes[S]);

    indexes_z[I] = map_index_z(indexes[I]);
    indexes_z[J] = map_index_z(indexes[J]);
    indexes_z[K] = map_index_z(indexes[K]);
    indexes_z[L] = map_index_z(indexes[L]);
    indexes_z[R] = map_index_z(indexes[R]);
    indexes_z[S] = map_index_z(indexes[S]);

    row[$(indexes[R], indexes[S])] += A[$(indexes[I], indexes[J])] * B[$(indexes[K], indexes[L])];
    row[$(indexes_x[R], indexes_x[S])] += A[$(indexes_x[I], indexes_x[J])] * B[$(indexes_x[K], indexes_x[L])];

    row[$$(indexes[R], indexes[S])] += C[$(indexes[I], indexes[J])] * D[$(indexes[K], indexes[L])];
    row[$$(indexes_y[R], indexes_y[S])] += C[$(indexes_y[I], indexes_y[J])] * D[$(indexes_y[K], indexes_y[L])];

    row[$$$(indexes[R], indexes[S])] += E[$(indexes[I], indexes[J])] * F[$(indexes[K], indexes[L])];
    row[$$$(indexes_z[R], indexes_z[S])] += E[$(indexes_z[I], indexes_z[J])] * F[$(indexes_z[K], indexes_z[L])];
}

int main(void) {

    return 0;
}