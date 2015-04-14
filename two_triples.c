#include "ranlib/ranlib.h"
#include "utils.h"

#define N_VARS 18


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

    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

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


void fill_row_bm(lapack_complex_double *row,
                 lapack_complex_double *A,
                 lapack_complex_double *C,
                 lapack_complex_double *K,
                 lapack_complex_double *N,
                 int i, int j, int k, int l, int r, int s) {

    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

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


void fill_row_ak(lapack_complex_double *row,
                 lapack_complex_double *B,
                 lapack_complex_double *C,
                 lapack_complex_double *M,
                 lapack_complex_double *N,
                 int i, int j, int k, int l, int r, int s) {

    memset(row, 0, sizeof(lapack_complex_double) * N_VARS);

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


void print_all(lapack_complex_double *A,
               lapack_complex_double *B,
               lapack_complex_double *C,
               lapack_complex_double *Q,
               lapack_complex_double *O,
               lapack_complex_double *P,
               double residual) {
    printf("\n\n################# NEW ITERATION #################\n");

    printf("------- MATRIX A -------\n");
    print_matrix(A, 3, 3);
    printf("\n------- MATRIX B -------\n");
    print_matrix(B, 3, 3);
    printf("\n------- MATRIX C -------\n");
    print_matrix(C, 3, 3);
    printf("\n------- MATRIX K -------\n");
    print_matrix(Q, 3, 3);
    printf("\n------- MATRIX M -------\n");
    print_matrix(O, 3, 3);
    printf("\n------- MATRIX N -------\n");
    print_matrix(P, 3, 3);

    printf("RESIDUAL = %f\n", residual);
}


int main(void) {
    int i;
    int iterate = 0;
    int main_count = 2000;

    double prev_residual = 0.0;
    double cur_residual = 0.0;

    lapack_complex_double *matrix1, *matrix2, *matrix3, *matrix4;
    lapack_complex_double *A, *B, *C, *K, *M, *N;
    lapack_complex_double *MAIN_MATRIX, *MAIN_VECTOR;

    void (*generate_func[3])(lapack_complex_double *,
                             lapack_complex_double *,
                             lapack_complex_double *,
                             lapack_complex_double *,
                             lapack_complex_double *,
                             int, int, int, int, int, int);

    int **sets = read_sets();

    srand(time(0));
    setall(time(0), time(0));

    A = create_matrix(3, 3);
    B = create_matrix(3, 3);
    C = create_matrix(3, 3);
    K = create_matrix(3, 3);
    M = create_matrix(3, 3);
    N = create_matrix(3, 3);

    fill_matrix(A, 5);
    fill_matrix(B, 5);
    fill_matrix(C, 5);
    fill_matrix(K, 5);
    fill_matrix(M, 5);
    fill_matrix(N, 5);

    MAIN_MATRIX = (lapack_complex_double *)calloc(N_EQUATIONS * N_VARS, sizeof(lapack_complex_double));
    MAIN_VECTOR = (lapack_complex_double *)calloc(N_EQUATIONS, sizeof(lapack_complex_double));

    generate_func[0] = fill_row_cn;
    generate_func[1] = fill_row_ak;
    generate_func[2] = fill_row_bm;

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

        for (i = 0; i < N_EQUATIONS; ++i) {
            generate_func[matrix_index](&MAIN_MATRIX[i * N_VARS], matrix1, matrix2, matrix3, matrix4,
                                        sets[i][0], sets[i][1], sets[i][2],
                                        sets[i][3], sets[i][4], sets[i][5]);
            MAIN_VECTOR[i] = create_lhs_vector(sets[i]);
        }

        int info;
        info = LAPACKE_zgels(LAPACK_ROW_MAJOR, 'N', N_EQUATIONS, N_VARS, 1, MAIN_MATRIX, N_VARS, MAIN_VECTOR, 1);

        if (info > 0) {
            print_all(A, B, C, K, M, N, cur_residual);
            perror("The solution could not be computed (the matrix have not the full rank)");
            exit(1);
        }

        if (matrix_index == 0) {
            fill_matrix_new(C, &MAIN_VECTOR[0]);
            fill_matrix_new(N, &MAIN_VECTOR[9]);
        } else if (matrix_index == 1) {
            fill_matrix_new(A, &MAIN_VECTOR[0]);
            fill_matrix_new(K, &MAIN_VECTOR[9]);
        } else {
            fill_matrix_new(B, &MAIN_VECTOR[0]);
            fill_matrix_new(M, &MAIN_VECTOR[9]);
        }

        normalization(A, B, C);
        normalization(K, M, N);

        prev_residual = cur_residual;
        cur_residual = compute_residual(MAIN_VECTOR, N_VARS);

        iterate++;

        if (fabs(cur_residual - prev_residual) < 0.00000005) {
            iterate = 0;

            print_all(A, B, C, K, M, N, cur_residual);

            fill_matrix(A, 20);
            fill_matrix(B, 20);
            fill_matrix(C, 20);
            fill_matrix(K, 20);
            fill_matrix(M, 20);
            fill_matrix(N, 20);

            prev_residual = cur_residual = 0.0;
            main_count--;
        }
    } while (main_count != 0);

    print_all(A, B, C, K, M, N, cur_residual);

    free_matrix(A);
    free_matrix(B);
    free_matrix(C);
    free_matrix(K);
    free_matrix(M);
    free_matrix(N);
    free_matrix(MAIN_MATRIX);
    free_matrix(MAIN_VECTOR);
    free_sets(sets);

    return 0;
}