#include "ranlib/ranlib.h"
#include "utils.h"

#define N_VARS 9


inline int map_index_x1(int i) {
    switch (i) {
        case 0: return 0;
        case 1: return 1;
        case 2: return 2;
    }
    return i;
}


inline int map_index_x2(int i) {
    switch (i) {
        case 0: return 1;
        case 1: return 0;
        case 2: return 2;
    }
    return i;
}


inline int map_index_y1(int i) {
    switch (i) {
        case 0: return 2;
        case 1: return 1;
        case 2: return 0;
    }
    return i;
}


inline int map_index_y2(int i) {
    switch (i) {
        case 0: return 0;
        case 1: return 2;
        case 2: return 1;
    }
    return i;
}


inline int map_index_z1(int i) {
    switch (i) {
        case 0: return 1;
        case 1: return 2;
        case 2: return 0;
    }
    return i;
}


inline int map_index_z2(int i) {
    switch (i) {
        case 0: return 2;
        case 1: return 0;
        case 2: return 1;
    }
    return i;
}


void fill_row_c(lapack_complex_double *row,
                 lapack_complex_double *A,
                 lapack_complex_double *B,
                 int index[]) {
    int i;

    int index_x1[6], index_x2[6], index_y1[6], index_y2[6], index_z1[6], index_z2[6];

    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

    for (i = 0; i < 6; ++i) {
        index_x1[i] = map_index_x1(index[i]);
        index_x2[i] = map_index_x2(index[i]);
        index_y1[i] = map_index_y1(index[i]);
        index_y2[i] = map_index_y2(index[i]);
        index_z1[i] = map_index_z1(index[i]);
        index_z2[i] = map_index_z2(index[i]);
    }

    row[$(index_x1[iR], index_x1[iS])] += A[$(index_x1[iI], index_x1[iJ])] * B[$(index_x1[iK], index_x1[iL])];
    row[$(index_x2[iR], index_x2[iS])] += A[$(index_x2[iI], index_x2[iJ])] * B[$(index_x2[iK], index_x2[iL])];


    row[$(index_y1[iR], index_y1[iS])] += A[$(index_y1[iI], index_y1[iJ])] * B[$(index_y1[iK], index_y1[iL])];
    row[$(index_y2[iR], index_y2[iS])] += A[$(index_y2[iI], index_y2[iJ])] * B[$(index_y2[iK], index_y2[iL])];


    row[$(index_z1[iR], index_z1[iS])] += A[$(index_z1[iI], index_z1[iJ])] * B[$(index_z1[iK], index_z1[iL])];
    row[$(index_z2[iR], index_z2[iS])] += A[$(index_z2[iI], index_z2[iJ])] * B[$(index_z2[iK], index_z2[iL])];
}


void fill_row_b(lapack_complex_double *row,
                lapack_complex_double *A,
                lapack_complex_double *C,
                int index[]) {
    int i;

    int index_x1[6], index_x2[6], index_y1[6], index_y2[6], index_z1[6], index_z2[6];

    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

    for (i = 0; i < 6; ++i) {
        index_x1[i] = map_index_x1(index[i]);
        index_x2[i] = map_index_x2(index[i]);
        index_y1[i] = map_index_y1(index[i]);
        index_y2[i] = map_index_y2(index[i]);
        index_z1[i] = map_index_z1(index[i]);
        index_z2[i] = map_index_z2(index[i]);
    }

    row[$(index_x1[iK], index_x1[iL])] += A[$(index_x1[iI], index_x1[iJ])] * C[$(index_x1[iR], index_x1[iS])];
    row[$(index_x2[iK], index_x2[iL])] += A[$(index_x2[iI], index_x2[iJ])] * C[$(index_x2[iR], index_x2[iS])];


    row[$(index_y1[iK], index_y1[iL])] += A[$(index_y1[iI], index_y1[iJ])] * C[$(index_y1[iR], index_y1[iS])];
    row[$(index_y2[iK], index_y2[iL])] += A[$(index_y2[iI], index_y2[iJ])] * C[$(index_y2[iR], index_y2[iL])];


    row[$(index_z1[iK], index_z1[iL])] += A[$(index_z1[iI], index_z1[iJ])] * C[$(index_z1[iR], index_z1[iS])];
    row[$(index_z2[iK], index_z2[iL])] += A[$(index_z2[iI], index_z2[iJ])] * C[$(index_z2[iR], index_z2[iS])];
}

void fill_row_a(lapack_complex_double *row,
                lapack_complex_double *B,
                lapack_complex_double *C,
                int index[]) {
    int i;

    int index_x1[6], index_x2[6], index_y1[6], index_y2[6], index_z1[6], index_z2[6];

    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

    for (i = 0; i < 6; ++i) {
        index_x1[i] = map_index_x1(index[i]);
        index_x2[i] = map_index_x2(index[i]);
        index_y1[i] = map_index_y1(index[i]);
        index_y2[i] = map_index_y2(index[i]);
        index_z1[i] = map_index_z1(index[i]);
        index_z2[i] = map_index_z2(index[i]);
    }

    row[$(index_x1[iI], index_x1[iJ])] += B[$(index_x1[iK], index_x1[iL])] * C[$(index_x1[iR], index_x1[iS])];
    row[$(index_x2[iI], index_x2[iJ])] += B[$(index_x2[iK], index_x2[iL])] * C[$(index_x2[iR], index_x2[iS])];


    row[$(index_y1[iI], index_y1[iJ])] += B[$(index_y1[iK], index_y1[iL])] * C[$(index_y1[iR], index_y1[iS])];
    row[$(index_y2[iI], index_y2[iJ])] += B[$(index_y2[iK], index_y2[iL])] * C[$(index_y2[iR], index_y2[iL])];


    row[$(index_z1[iI], index_z1[iJ])] += B[$(index_z1[iK], index_z1[iL])] * C[$(index_z1[iR], index_z1[iS])];
    row[$(index_z2[iI], index_z2[iJ])] += B[$(index_z2[iK], index_z2[iL])] * C[$(index_z2[iR], index_z2[iS])];
}


void print_all(lapack_complex_double *A,
               lapack_complex_double *B,
               lapack_complex_double *C,
               double residual) {
    printf("\n\n################# NEW ITERATION #################\n");

    printf("------- MATRIX A -------\n");
    print_matrix(A, 3, 3);
    printf("\n------- MATRIX B -------\n");
    print_matrix(B, 3, 3);
    printf("\n------- MATRIX C -------\n");
    print_matrix(C, 3, 3);

    printf("RESIDUAL = %f\n", residual);
}

int main(void) {

    int i;
    int iterate = 0;
    int main_count = 60000;

    double prev_residual = 0.0;
    double cur_residual = 0.0;

    lapack_complex_double *matrix1, *matrix2;
    lapack_complex_double *A, *B, *C;
    lapack_complex_double *MAIN_MATRIX, *MAIN_VECTOR;

    void (*generate_func[3])(lapack_complex_double *,
                             lapack_complex_double *,
                             lapack_complex_double *,
                             int[]);

    int **sets = read_sets();

    srand(time(0));
    setall(time(0), time(0));

    A = create_matrix(3, 3);
    B = create_matrix(3, 3);
    C = create_matrix(3, 3);

    fill_matrix(A, 10);
    fill_matrix(B, 10);
    fill_matrix(C, 10);

    MAIN_MATRIX = (lapack_complex_double *)calloc(N_EQUATIONS * N_VARS, sizeof(lapack_complex_double));
    MAIN_VECTOR = (lapack_complex_double *)calloc(N_EQUATIONS, sizeof(lapack_complex_double));

    generate_func[0] = fill_row_c;
    generate_func[1] = fill_row_a;
    generate_func[2] = fill_row_b;

    do {
        int matrix_index = iterate % 3;

        if (matrix_index == 0) {
            matrix1 = A;
            matrix2 = B;
        } else if (matrix_index == 1) {
            matrix1 = B;
            matrix2 = C;
        } else {
            matrix1 = A;
            matrix2 = C;
        }

        for (i = 0; i < N_EQUATIONS; ++i) {
            generate_func[matrix_index](&MAIN_MATRIX[i * N_VARS], matrix1, matrix2, sets[i]);
            MAIN_VECTOR[i] = create_lhs_vector(sets[i]);
        }

        int info;
        info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', N_EQUATIONS, N_VARS, 1, MAIN_MATRIX, N_VARS, MAIN_VECTOR, 1);

        if (info > 0) {
            print_all(A, B, C, cur_residual);
            perror("The solution could not be computed (the matrix have not the full rank)");
            exit(1);
        }

        if (matrix_index == 0) {
            fill_matrix_new(C, MAIN_VECTOR);
        } else if (matrix_index == 1) {
            fill_matrix_new(A, MAIN_VECTOR);
        } else {
            fill_matrix_new(B, MAIN_VECTOR);
        }

        normalization(A, B, C);

        prev_residual = cur_residual;
        cur_residual = compute_residual(MAIN_VECTOR, N_VARS);

        iterate++;

        if (iterate < 30000 || cur_residual > 1.91) {
            iterate = 0;

            print_all(A, B, C, cur_residual);

            fill_matrix(A, 10);
            fill_matrix(B, 10);
            fill_matrix(C, 10);

            prev_residual = cur_residual = 0.0;
            main_count--;
        }
    } while (main_count != 0);

    print_all(A, B, C, cur_residual);

    free_matrix(A);
    free_matrix(B);
    free_matrix(C);
    free_matrix(MAIN_MATRIX);
    free_matrix(MAIN_VECTOR);
    free_sets(sets);

    return 0;
}



