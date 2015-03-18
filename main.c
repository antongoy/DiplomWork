#include "ranlib/ranlib.h"

#include <lapacke.h>

#include <stdio.h>
#include <string.h>
#include <time.h>

#define $(i, j) ((i) * 3 + (j))
#define $$(i, j) ((i) * 3 + (j) + 9)

#define N_EQUATIONS 243

inline lapack_complex_double * create_matrix(size_t w, size_t l) {
    return (lapack_complex_double *)calloc(w * l, sizeof(lapack_complex_double));
}


void fill_matrix(lapack_complex_double *matrix) {
    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            //matrix[$(i, j)] = lapack_make_complex_double(ranf(), ranf());
            matrix[$(i, j)] = lapack_make_complex_double(1, 1);
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

void read_sets(int **sets) {
    int i;
    FILE *f;

    if((f = fopen("sets.txt", "r")) < 0) {
        perror("File reading error:");
        exit(1);
    }

    for (i = 0; i < N_EQUATIONS; ++i) {
        fscanf(f, "%d %d %d %d %d %d", &sets[i][0], &sets[i][1], &sets[i][2], &sets[i][3], &sets[i][4], &sets[i][5]);
    }

    for (i = 0; i < N_EQUATIONS; ++i) {
        sets[i][0]--;
        sets[i][1]--;
        sets[i][2]--;
        sets[i][3]--;
        sets[i][4]--;
        sets[i][5]--;
    }
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

    lapack_complex_double *A, *B, *C, *K, *M, *N;
    lapack_complex_double *MAIN_MATRIX;
    int **sets = (int **)calloc(N_EQUATIONS, sizeof(int *));

    for(i = 0; i < N_EQUATIONS; i++) {
        sets[i] = (int *)calloc(6, sizeof(int));
    }

    read_sets(sets);
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

    MAIN_MATRIX = (lapack_complex_double *)calloc(N_EQUATIONS * 18, sizeof(lapack_complex_double));

    for (i = 0; i < N_EQUATIONS; ++i) {
        fill_row_cn(&MAIN_MATRIX[i * 18], A, B, K, M,
                    sets[i][0], sets[i][1], sets[i][2],
                    sets[i][3], sets[i][4], sets[i][5]);
    }

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