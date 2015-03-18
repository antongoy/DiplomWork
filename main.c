#include "ranlib/ranlib.h"

#include <lapacke.h>

#include <stdio.h>
#include <string.h>

#define $(i, j) ((i) * 3 + (j))
#define $$(i, j) ((i) * 3 + (j) + 9)

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

    int rcx, scx, icx, jcx, kcx, lcx;

    rcx = map_index_x(r);
    scx = map_index_x(s);
    icx = map_index_x(i);
    jcx = map_index_x(j);
    kcx = map_index_x(k);
    lcx = map_index_x(l);

    row[$(r, s)] = A[$(i, j)] * B[$(k, l)];
    row[$(rcx, scx)] = A[$(icx, jcx)] * B[$(kcx, lcx)];
    
    int rny, sny, iny, jny, kny, lny;

    rny = map_index_y(r);
    sny = map_index_y(s);
    iny = map_index_y(i);
    jny = map_index_y(j);
    kny = map_index_y(k);
    lny = map_index_y(l);
    
    row[$$(r, s)] = K[$(i, j)] * M[$(k, l)];
    row[$$(rny, sny)] = K[$(iny, jny)] * M[$(kny, lny)];

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

    row[$$(r, s)] = K[$(i, j)] * M[$(k, l)];
    row[$$(rnz, snz)] = K[$(inz, jnz)] * M[$(knz, lnz)];
    
    
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

    row[$(k, l)] = A[$(i, j)] * C[$(r, s)];
    row[$(kbx, lbx)] = A[$(ibx, jbx)] * C[$(rbx, sbx)];

    int rmy, smy, imy, jmy, kmy, lmy;

    rmy = map_index_y(r);
    smy = map_index_y(s);
    imy = map_index_y(i);
    jmy = map_index_y(j);
    kmy = map_index_y(k);
    lmy = map_index_y(l);

    row[$$(k, l)] = K[$(i, j)] * N[$(r, s)];
    row[$$(kmy, lmy)] = K[$(imy, jmy)] * N[$(rmy, smy)];

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

    row[$$(k, l)] = K[$(i, j)] * N[$(r, s)];
    row[$$(kmz, lmz)] = K[$(imz, jmz)] * N[$(rmz, smz)];
}


int main(void) {
    lapack_complex_double *A;

    A = create_matrix();

    fill_matrix(A);

    print_matrix(A);

    return 0;
}