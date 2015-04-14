//
// Created by anton-goy on 13.04.15.
//

#include "utils.h"
#include "ranlib/ranlib.h"


inline lapack_complex_double * create_matrix(size_t w, size_t l) {
    return (lapack_complex_double *)calloc(w * l, sizeof(lapack_complex_double));
}


void fill_matrix(lapack_complex_double *matrix, int max_abs_value) {
    int i, j;
    int c;
    double a, b;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j)  {
            c = rand() % max_abs_value + 1;
            a = (2.0 *ranf() - 1.0) * c;
            b = (2.0 *ranf() - 1.0) * c;
            matrix[$(i, j)] = lapack_make_complex_double(a, b);
        }
    }
}


inline void  print_matrix(lapack_complex_double *matrix, size_t w, size_t l) {
    int i, j;

    for (i = 0; i < w; ++i) {
        for (j = 0; j < l; ++j) {
            printf("(%.3g, %.3g)  ", lapack_complex_double_real(matrix[i * l + j]),
                   lapack_complex_double_imag(matrix[i * l + j]));
        }
        printf("\n");
    }
}


int** read_sets() {
    int i;
    FILE *f;

    int **sets = (int **)calloc(N_EQUATIONS, sizeof(int *));

    for(i = 0; i < N_EQUATIONS; i++) {
        sets[i] = (int *)calloc(6, sizeof(int));
    }

    if((f = fopen("sets.txt", "r")) < 0) {
        perror("File reading error:");
        exit(1);
    }

    for (i = 0; i < N_EQUATIONS; ++i) {
        int n = fscanf(f, "%d %d %d %d %d %d", &sets[i][0], &sets[i][1], &sets[i][2],
                       &sets[i][3], &sets[i][4], &sets[i][5]);
        if (n < 0) {
            perror("IO error");
            exit(-1);
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


lapack_complex_double create_lhs_vector(int *index) {
    int i, j, k, l, r, s;

    i = index[0];
    j = index[1];
    k = index[2];
    l = index[3];
    r = index[4];
    s = index[5];

    if (j == k && l == r && s == i) {
        if (i == j && k == l && r == s) {
            return lapack_make_complex_double(0, 0);
        } else {
            return lapack_make_complex_double(1.0 / 3.0, 0);
        }
    } else {
        if (i == j && k == l && r == s) {
            return lapack_make_complex_double(- 1.0 / 3.0, 0);
        } else {
            return lapack_make_complex_double(0, 0);
        }
    }
}


double compute_residual(lapack_complex_double *main_vector, size_t n_vars) {
    double residual = 0.0;
    int i;
    for (i = n_vars; i < N_EQUATIONS; ++i) {
        residual += pow(lapack_complex_float_real(main_vector[i]), 2) +
                    pow(lapack_complex_float_imag(main_vector[i]), 2);
    }

    return residual;
}


void fill_matrix_new(lapack_complex_double *dst, lapack_complex_double *src) {
    memcpy(dst, src, sizeof(lapack_complex_double) * 9);

    int i, j;
    double a = 0.0, b = 0.0;
    double alpha = -0.00001;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            a += pow(lapack_complex_double_real(dst[$(i, j)]), 2);
            b += pow(lapack_complex_double_imag(dst[$(i, j)]), 2);
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            dst[$(i, j)] += (alpha / (2 * 9)) * lapack_make_complex_double(a, b);
        }
    }

}


double get_max_abs_of_matrix(lapack_complex_double *matrix) {
    int i, j;
    double max_abs = -1;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            double elem_max = sqrt(pow(lapack_complex_double_real(matrix[$(i,j)]), 2.0) +
                                   pow(lapack_complex_double_imag(matrix[$(i,j)]), 2.0));
            if (elem_max > max_abs) {
                max_abs = elem_max;
            }
        }
    }
    return max_abs;
}


void scalar_multiplication(lapack_complex_double *matrix, double scalar) {
    int i, j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            matrix[$(i, j)] *= scalar;
        }
    }
}


void normalization(lapack_complex_double *A,
                   lapack_complex_double *B,
                   lapack_complex_double *C) {

    double max1, max2, max3;

    max1 = get_max_abs_of_matrix(A);
    max2 = get_max_abs_of_matrix(B);
    max3 = get_max_abs_of_matrix(C);

    double k = pow(max1 * max2 * max3, 1 / 3.);

    scalar_multiplication(A, k / max1);
    scalar_multiplication(B, k / max2);
    scalar_multiplication(C, k / max3);
}