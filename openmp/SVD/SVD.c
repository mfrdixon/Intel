/*
 * Sample code for Least-Squares Monte Carlo implementing the Singular value decomposition, 
 * using MKL/VSL random number generation and compiler vectorisation based
 * #pragma omp simd directives, plus OpenMP parallelisation
 *
 * author: Matthew Dixon and Xin Zhang
 *
 * Date: Nov 30 2016
 */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mkl.h>
#include <mkl_vsl.h>

#include "mkl_vsl.h"
#include "mkl_lapacke.h"
#include "errcheck.inc"

#define SEED    1
#define BRNG    VSL_BRNG_SFMT19937
#define METHOD  VSL_RNG_METHOD_GAUSSIAN_ICDF
#define N       100000
#define NN      10
#define min(a,b) ((a)>(b)?(a):(b))
extern void print_cmatrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda );
extern void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda );
#define M_ 400
#define N_ 250
#define LDA M_
#define LDU M_
#define LDVT N_


int main(){
    
    double r[N]; // Random numbers
    VSLStreamStatePtr stream; // Random stream
    int i, errcode;
    float a=0.0,b=1.0;
    double tM,tD,tQ,tD2;
    double sM,sD;
    double sum, sum2;
    double n,s;
    double DeltaM,DeltaD;
    
    // Random number set up
    errcode = vslNewStream( &stream, BRNG,  SEED );
    CheckVslError( errcode );
    // Random number generation
    errcode = vdRngGaussian( METHOD, stream, N, r, a, b );
    CheckVslError( errcode );

    tM=a;
    tD=b;    
    sum=0.0;
    sum2=0.0;
    for(i=0;i<N;i++) {
        sum+=(double)r[i];
        sum2+=(double)r[i]*(double)r[i];
    }
    sM=sum/((double)N);
    sD=sum2/(double)N-(sM*sM);
    n=(double)N;
    tD2=tD*tD;
    s=((tQ-tD2)/n)-(2*(tQ-2*tD2)/(n*n))+((tQ-3*tD2)/(n*n*n));
    DeltaM=(tM-sM)/sqrt(tD/n);
    DeltaD=(tD-sD)/sqrt(s);

    printf("Sample of vdRngGaussian.\n");
    printf("-----------------------\n\n");
    printf("Parameters:\n");
    printf("    a=%.4f\n",a);
    printf("    b=%.4f\n\n",b);
    printf("Results (first NN of N):\n");
    printf("---------------------------\n");
    for(i=0;i<NN;i++) {
        printf("r[%d]=%.4f\n",i,r[i]);
    }
    printf("\n");
    if(fabs(DeltaM)>3.0 || fabs(DeltaD)>3.0) {
        printf("Error: sample moments (mean=%.2f, variance=%.2f) disagree with theory (mean=%.2f, variance=%.2f).\n",sM,sD,tM,tD);
        return 1;
    }
    else {
        printf("Sample moments (mean=%.2f, variance=%.2f) agree with theory (mean=%.2f, variance=%.2f).\n",sM,sD,tM,tD);
    }
    errcode = vslDeleteStream( &stream );

/***********************SVD*********************/
    //repackage stream into matrix M_*N_
    MKL_INT  m_=M_, n_=N_, lda_ = LDA, ldu_ = LDU,ldvt_=LDVT,info;
    float s_[M_];
    float superb[min(M_,N_)-1];
    MKL_Complex8 u_[LDU*M_],vt_[LDVT*N_];
    //populate matrix MKL_COmplex8 a[LDA*N_] by using random numbers
    MKL_Complex8 a_[LDA*N_];
    int j=0;
    for(j=0; j<LDA*N_; j++){
    	a_[j].real = r[j];
    	a_[j].imag = 0;
    }
    printf("LAPACKE_cgesvd(column-major, high level) runs....\n");
    // Singular value decomposition implementation
    info = LAPACKE_cgesvd( LAPACK_COL_MAJOR,'A','A',m_,n_,a_,lda_,s_,u_,ldu_,vt_,ldvt_,superb);    
    if(info>0)
      {
		printf("The algorithm failed to converge \n");
		exit(1);
      }
    // Print out the results
    print_rmatrix( "Singular values", 1, m_, s_, 1 );
    return 0;
}

/* Auxiliary routine: printing a complex matrix */
void print_cmatrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", a[i+j*lda].real, a[i+j*lda].imag );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
                printf( "\n" );
        }
}
