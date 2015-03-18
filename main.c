#include "ranlib/ranlib.h"

#include <lapacke.h>

#include <stdio.h>
#include <string.h>

#define $(i, j) ((i) * 3 + (j))


inline lapack_complex_double * create_matrix(void) {
    return (lapack_complex_double *)calloc(3 * 3, sizeof(lapack_complex_double));
}


void fill_matrix(lapack_complex_double *matrix) {
    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            matrix[$(i, j)] = lapack_make_complex_double(ranf(), ranf());
        }
    }
}


inline void print_matrix(lapack_complex_double *matrix) {
    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            printf("%.3g + %.3g * I\t", lapack_complex_double_real(matrix[$(i,j)]),
                                        lapack_complex_double_imag(matrix[$(i,j)]));
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


inline int map_index_z1(int i) {
    return map_index_x(i);
}


inline int map_index_z2(int i) {
    switch (i) {
        case 0: return 2;
        case 1: return 0;
        case 2: return 1;
    }
    return i;
}



void fill_row_cn(lapack_complex_double *row,
                 lapack_complex_double *A,
                 lapack_complex_double *B,
                 lapack_complex_double *K,
                 lapack_complex_double *M,
                 int i, int j, int k, int l, int r, int s) {
    int row_size = ((3 * 3) * 2);
    memset(row, 0, sizeof(lapack_complex_double) * row_size);


}


int main(void) {
    lapack_complex_double *A;

    A = create_matrix();

    fill_matrix(A);

    print_matrix(A);

    return 0;
}