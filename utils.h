//
// Created by anton-goy on 13.04.15.
//

#ifndef DIPLOMWORK_UTILS_H
#define DIPLOMWORK_UTILS_H

#include <lapacke.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define $(i, j) ((i) * 3 + (j))
#define $$(i, j) ((i) * 3 + (j) + 9)
#define $$$(i, j) ((i) * 3 + (j) + 18)

#define N_EQUATIONS 243

#define iI 0
#define iJ 1
#define iK 2
#define iL 3
#define iR 4
#define iS 5

inline lapack_complex_double * create_matrix(size_t, size_t);
void fill_matrix(lapack_complex_double *, int);
inline void  print_matrix(lapack_complex_double *, size_t, size_t);
int **read_sets();
void free_matrix(lapack_complex_double *);
void free_sets(int **);
lapack_complex_double create_lhs_vector(int *);
double compute_residual(lapack_complex_double *, size_t);
void fill_matrix_new(lapack_complex_double *, lapack_complex_double *);
double get_max_abs_of_matrix(lapack_complex_double *);
void scalar_multiplication(lapack_complex_double *, double);
void normalization(lapack_complex_double *, lapack_complex_double *, lapack_complex_double *);


#endif //DIPLOMWORK_UTILS_H
