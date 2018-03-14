#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <chrono>
#include <random>
#include "nrutil.h"
void funcs( double x, double *afunc, int ma );
void svdfit( double *X, double *Y, double *Sig, int NData, double *A, int MA,
    double **U, double **V, double *W, double *ChiSq, void funcs(double x, double *afunc, int ma));
void svdcmp( double **A, int M, int N, double *W, double **V );
void svbksb( double **U, double *W, double **V, int M,
    int N, double *B, double *X );
void svdvar( double **V, int MA, double *W, double **CVM );
void fleg(double x, double pl[], int nl);

# define true ((int)1)
# define false ((int)0)
# define nmax ((int)1000)
# define mmax ((int)50)
# define tol ((double)1.0e-15)
