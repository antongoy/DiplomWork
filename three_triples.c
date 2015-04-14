//
// Created by anton-goy on 13.04.15.
//

#include "ranlib/ranlib.h"
#include "utils.h"

#define N_VARS 27


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


void fill_row_cpw(lapack_complex_double *row,
                    lapack_complex_double *A,
                    lapack_complex_double *B,
                    lapack_complex_double *C,
                    lapack_complex_double *D,
                    lapack_complex_double *E,
                    lapack_complex_double *F,
                    int index[]) {

    int i;
    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

    int index_x[6], index_y[6], index_z[6];

    for (i = 0; i < 6; ++i) {
        index_x[i] = map_index_x(index[i]);
        index_y[i] = map_index_y(index[i]);
        index_z[i] = map_index_z(index[i]);
    }

    row[$(index[iR], index[iS])] += A[$(index[iI], index[iJ])] * B[$(index[iK], index[iL])];
    row[$(index_x[iR], index_x[iS])] += A[$(index_x[iI], index_x[iJ])] * B[$(index_x[iK], index_x[iL])];

    row[$$(index[iR], index[iS])] += C[$(index[iI], index[iJ])] * D[$(index[iK], index[iL])];
    row[$$(index_y[iR], index_y[iS])] += C[$(index_y[iI], index_y[iJ])] * D[$(index_y[iK], index_y[iL])];

    row[$$$(index[iR], index[iS])] += E[$(index[iI], index[iJ])] * F[$(index[iK], index[iL])];
    row[$$$(index_z[iR], index_z[iS])] += E[$(index_z[iI], index_z[iJ])] * F[$(index_z[iK], index_z[iL])];
}

void fill_row_aqu(lapack_complex_double *row,
                    lapack_complex_double *A,
                    lapack_complex_double *B,
                    lapack_complex_double *C,
                    lapack_complex_double *D,
                    lapack_complex_double *E,
                    lapack_complex_double *F,
                    int index[]) {

    int i;
    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

    int index_x[6], index_y[6], index_z[6];

    for (i = 0; i < 6; ++i) {
        index_x[i] = map_index_x(index[i]);
        index_y[i] = map_index_y(index[i]);
        index_z[i] = map_index_z(index[i]);
    }

    row[$(index[iI], index[iJ])] += A[$(index[iK], index[iL])] * B[$(index[iR], index[iS])];
    row[$(index_x[iI], index_x[iJ])] += A[$(index_x[iK], index_x[iL])] * B[$(index_x[iR], index_x[iS])];

    row[$$(index[iI], index[iJ])] += C[$(index[iK], index[iL])] * D[$(index[iR], index[iS])];
    row[$$(index_y[iI], index_y[iJ])] += C[$(index_y[iK], index_y[iL])] * D[$(index_y[iR], index_y[iS])];

    row[$$$(index[iI], index[iJ])] += E[$(index[iK], index[iL])] * F[$(index[iR], index[iS])];
    row[$$$(index_z[iI], index_z[iJ])] += E[$(index_z[iK], index_z[iL])] * F[$(index_z[iR], index_z[iS])];
}

void fill_row_bov(lapack_complex_double *row,
                    lapack_complex_double *A,
                    lapack_complex_double *B,
                    lapack_complex_double *C,
                    lapack_complex_double *D,
                    lapack_complex_double *E,
                    lapack_complex_double *F,
                    int index[]) {

    int i;
    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

    int index_x[6], index_y[6], index_z[6];

    for (i = 0; i < 6; ++i) {
        index_x[i] = map_index_x(index[i]);
        index_y[i] = map_index_y(index[i]);
        index_z[i] = map_index_z(index[i]);
    }

    row[$(index[iK], index[iL])] += A[$(index[iI], index[iJ])] * B[$(index[iR], index[iS])];
    row[$(index_x[iK], index_x[iL])] += A[$(index_x[iI], index_x[iJ])] * B[$(index_x[iR], index_x[iS])];

    row[$$(index[iK], index[iL])] += C[$(index[iI], index[iJ])] * D[$(index[iR], index[iS])];
    row[$$(index_y[iK], index_y[iL])] += C[$(index_y[iI], index_y[iJ])] * D[$(index_y[iR], index_y[iS])];

    row[$$$(index[iK], index[iL])] += A[$(index[iI], index[iJ])] * B[$(index[iR], index[iS])];
    row[$$$(index_z[iK], index_z[iL])] += E[$(index_z[iI], index_z[iJ])] * F[$(index_z[iR], index_z[iS])];
}

void print_all(lapack_complex_double *A,
               lapack_complex_double *B,
               lapack_complex_double *C,
               lapack_complex_double *Q,
               lapack_complex_double *O,
               lapack_complex_double *P,
               lapack_complex_double *U,
               lapack_complex_double *V,
               lapack_complex_double *W,
               double residual) {
    printf("\n\n################# NEW ITERATION #################\n");

    printf("------- MATRIX A -------\n");
    print_matrix(A, 3, 3);
    printf("\n------- MATRIX B -------\n");
    print_matrix(B, 3, 3);
    printf("\n------- MATRIX C -------\n");
    print_matrix(C, 3, 3);
    printf("\n------- MATRIX Q -------\n");
    print_matrix(Q, 3, 3);
    printf("\n------- MATRIX O -------\n");
    print_matrix(O, 3, 3);
    printf("\n------- MATRIX P -------\n");
    print_matrix(P, 3, 3);
    printf("\n------- MATRIX U -------\n");
    print_matrix(U, 3, 3);
    printf("\n------- MATRIX V -------\n");
    print_matrix(V, 3, 3);
    printf("\n------- MATRIX W -------\n");
    print_matrix(W, 3, 3);

    printf("RESIDUAL = %f\n", residual);

}

int main(void) {
    int i;
    int iterate = 0, main_count = 2000;
    double prev_residual = 0.0, cur_residual = 0.0;
    lapack_complex_double *matrix1, *matrix2, *matrix3, *matrix4, *matrix5, *matrix6;

    void (*generate_func[3])(lapack_complex_double *,
                             lapack_complex_double *,
                             lapack_complex_double *,
                             lapack_complex_double *,
                             lapack_complex_double *,
                             lapack_complex_double *,
                             lapack_complex_double *,
                             int []);

    generate_func[0] = fill_row_cpw;
    generate_func[1] = fill_row_aqu;
    generate_func[2] = fill_row_bov;

    srand(time(0));
    setall(time(0), time(0));

    lapack_complex_double *A, *B, *C, *Q, *O, *P, *U, *V, *W;
    lapack_complex_double *MAIN_MATRIX, *MAIN_VECTOR;

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

    fill_matrix(A, 10);
    fill_matrix(B, 10);
    fill_matrix(C, 10);
    fill_matrix(Q, 10);
    fill_matrix(O, 10);
    fill_matrix(P, 10);
    fill_matrix(U, 10);
    fill_matrix(V, 10);
    fill_matrix(W, 10);

    MAIN_MATRIX = (lapack_complex_double *)calloc(N_EQUATIONS * N_VARS, sizeof(lapack_complex_double));
    MAIN_VECTOR = (lapack_complex_double *)calloc(N_EQUATIONS, sizeof(lapack_complex_double));

    do {
        int matrix_index = iterate % 3;

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
            generate_func[matrix_index](&MAIN_MATRIX[i * N_VARS], matrix1, matrix2, matrix3,
                                        matrix4, matrix5, matrix6, sets[i]);
            MAIN_VECTOR[i] = create_lhs_vector(sets[i]);
        }

        int info;
        info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', N_EQUATIONS, N_VARS, 1, MAIN_MATRIX, N_VARS, MAIN_VECTOR, 1);

        if (info > 0) {
            perror("The solution could not be computed (the matrix have not the full rank)");
            printf("info = %d\n", info);

            print_all(A, B, C, Q, O, P, U, V, W, cur_residual);

            exit(1);
        }

        if (matrix_index == 0) {
            fill_matrix_new(C, &MAIN_VECTOR[0]);
            fill_matrix_new(P, &MAIN_VECTOR[9]);
            fill_matrix_new(W, &MAIN_VECTOR[18]);
        } else if (matrix_index == 1) {
            fill_matrix_new(A, &MAIN_VECTOR[0]);
            fill_matrix_new(Q, &MAIN_VECTOR[9]);
            fill_matrix_new(U, &MAIN_VECTOR[18]);
        } else {
            fill_matrix_new(B, &MAIN_VECTOR[0]);
            fill_matrix_new(O, &MAIN_VECTOR[9]);
            fill_matrix_new(V, &MAIN_VECTOR[18]);
        }

        normalization(A, B, C);
        normalization(Q, O, P);
        normalization(U, V, W);

        prev_residual = cur_residual;
        cur_residual = compute_residual(MAIN_VECTOR, N_VARS);

        iterate++;

        //printf("RESIDUAL = %g\n", cur_residual);

        if (iterate > 1000 && cur_residual > 1.1) {
            iterate = 0;

            print_all(A, B, C, Q, O, P, U, V, W, cur_residual);

            prev_residual = cur_residual = 0.0;

            fill_matrix(A, 10);
            fill_matrix(B, 10);
            fill_matrix(Q, 10);
            fill_matrix(O, 10);
            fill_matrix(P, 10);
            fill_matrix(U, 10);
            fill_matrix(V, 10);
            fill_matrix(W, 10);

            main_count--;
        }

    } while (main_count != 0);

    free_matrix(A);
    free_matrix(B);
    free_matrix(C);
    free_matrix(Q);
    free_matrix(O);
    free_matrix(P);
    free_matrix(U);
    free_matrix(V);
    free_matrix(W);
    free_matrix(MAIN_MATRIX);
    free_matrix(MAIN_VECTOR);
    free_sets(sets);

    return 0;
}


