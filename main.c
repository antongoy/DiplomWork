#include "ranlib/ranlib.h"

#include <lapacke.h>
#include <lapacke_utils.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <gmpxx.h>

#define $(i, j) ((i) * 3 + (j))
#define $$(i, j) ((i) * 3 + (j) + 9)
#define $$$(i, j) ((i) * 3 + (j) + 18)

#define N_EQUATIONS 243
#define N_ATTEMPTS 1
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


lapack_complex_double fillin_rhs(int indexes[]) {
    if (indexes[J] == indexes[K] &&
        indexes[L] == indexes[R] &&
        indexes[S] == indexes[I]) {

        if (indexes[I] == indexes[J] &&
            indexes[K] == indexes[L] &&
            indexes[R] == indexes[S]) {

            return lapack_make_complex_double(0, 0);
        } else {
            return lapack_make_complex_double(1.0 / 3.0, 0);
        }

    } else {
        if (indexes[I] == indexes[J] &&
            indexes[K] == indexes[L] &&
            indexes[R] == indexes[S]) {

            return lapack_make_complex_double(-1.0 / 3.0, 0);
        } else {
            return lapack_make_complex_double(0, 0);
        }
    }
}


double compute_residual(lapack_complex_double *main_vector) {
    double residual = 0.0;
    int i;
    for (i = N_VARS; i < N_EQUATIONS; ++i) {
        residual += pow(lapack_complex_float_real(main_vector[i]), 2) +
                pow(lapack_complex_float_imag(main_vector[i]), 2);
    }

    return residual;
}

int **read_sets() {
    int i;
    FILE *f;

    int **sets = (int **)calloc(N_EQUATIONS, sizeof(int *));

    for(i = 0; i < N_EQUATIONS; i++) {
        sets[i] = (int *)calloc(6, sizeof(int));
    }

    if((f = fopen("sets.txt", "r")) < 0) {
        perror("File opening error:");
        exit(1);
    }

    for (i = 0; i < N_EQUATIONS; ++i) {
        int n = fscanf(f, "%d %d %d %d %d %d", &sets[i][0],
                &sets[i][1], &sets[i][2], &sets[i][3], &sets[i][4], &sets[i][5]);

        if (n < 0) {
            perror("File reading error:");
            exit(1);
        }
    }

    for (i = 0; i < N_EQUATIONS; ++i) {
        sets[i][0]--;
        sets[i][1]--;
        sets[i][2]--;
        sets[i][3]--;
        sets[i][4]--;
        sets[i][5]--;
    }

    fclose(f);

    return sets;
}


void free_matrix(lapack_complex_double *matrix) {
    free(matrix);
}


void free_sets(int **sets) {
    int i;

    for(i = 0; i < N_EQUATIONS; i++) {
        free(sets[i]);
    }

    free(sets);
}

int main(void) {
    int i;
    int internal_iteration = 0, external_iteration = 0;
    double prev_residual = 0.0, cur_residual = 0.0;
    lapack_complex_double *matrix1, *matrix2, *matrix3, *matrix4, *matrix5, *matrix6;

    srand(time(0));
    setall(time(0), time(0));

    lapack_complex_double *A, *B, *C, *Q, *O, *P, *U, *V, *W;
    lapack_complex_double *MATRIX, *VECTOR;

    int **sets = read_sets();

    A = create_matrix(3, 3);
    B = create_matrix(3, 3);
    C = create_matrix(3, 3);
    Q = create_matrix(3, 3);
    O = create_matrix(3, 3);
    P = create_matrix(3, 3);
    U = create_matrix(3, 3);
    V = create_matrix(3, 3);
    W = create_matrix(3, 3);

    fillin_matrix(A);
    fillin_matrix(B);
    fillin_matrix(C);
    fillin_matrix(Q);
    fillin_matrix(O);
    fillin_matrix(P);
    fillin_matrix(U);
    fillin_matrix(V);
    fillin_matrix(W);

    MATRIX = (lapack_complex_double *)calloc(N_EQUATIONS * N_VARS, sizeof(lapack_complex_double));
    VECTOR = (lapack_complex_double *)calloc(N_EQUATIONS, sizeof(lapack_complex_double));

    do {
        int matrix_index = internal_iteration % 3;

        if (matrix_index == 0) {
            matrix1 = A;
            matrix2 = B;
            matrix3 = Q;
            matrix4 = O;
            matrix5 = U;
            matrix6 = V;
        } else if (matrix_index == 1) {
            matrix1 = B;
            matrix2 = C;
            matrix3 = O;
            matrix4 = P;
            matrix5 = V;
            matrix6 = W;
        } else {
            matrix1 = A;
            matrix2 = C;
            matrix3 = O;
            matrix4 = Q;
            matrix5 = U;
            matrix6 = W;
        }

        for (i = 0; i < N_EQUATIONS; ++i) {
            fillin_row(&MATRIX[i * N_VARS], matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, sets[i]);
            VECTOR[i] = fillin_rhs(sets[i]);
        }

        int info;
        info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', N_EQUATIONS, N_VARS, 1, MATRIX, N_VARS, VECTOR, 1);

        if (info > 0) {
            perror("The solution could not be computed (the matrix have not the full rank)");
            exit(1);
        }

        if (matrix_index == 0) {
            fill_matrix_new(C, &VECTOR[0]);
            fill_matrix_new(P, &VECTOR[9]);
            fill_matrix_new(W, &VECTOR[18]);
        } else if (matrix_index == 1) {
            fill_matrix_new(A, &VECTOR[0]);
            fill_matrix_new(Q, &VECTOR[9]);
            fill_matrix_new(U, &VECTOR[18]);
        } else {
            fill_matrix_new(B, &VECTOR[0]);
            fill_matrix_new(O, &VECTOR[9]);
            fill_matrix_new(V, &VECTOR[18]);
        }

        external_iteration++;
    } while (external_iteration < N_ATTEMPTS);







    return 0;
}