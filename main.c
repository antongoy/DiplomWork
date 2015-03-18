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


int main(void) {
    lapack_complex_double *A;

    A = create_matrix();

    fill_matrix(A);

    return 0;
}