#pragma once

#include <mkl_lapacke.h>

namespace lapacke {
    // "11. If you have a procedure with ten parameters, you probably missed some."
    //      -- A. J. Perlis
    template<typename T>
    int gesvd(int, char, char, int, int, T*, int, T*, T*, int, T*, int, T*);

    template<>
    int gesvd<float>(int matrix_layout, char jobu, char jobvt, int m, int n,
                     float *a, int lda, float *s, float *u, int ldu,
                     float *vt, int ldvt, float *superb) {
        return LAPACKE_sgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu,
                                vt, ldvt, superb);
    }

    template<>
    int gesvd<double>(int matrix_layout, char jobu, char jobvt, int m, int n,
                      double *a, int lda, double *s, double *u, int ldu,
                      double *vt, int ldvt, double *superb) {
        return LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu,
                                vt, ldvt, superb);
    }
}

