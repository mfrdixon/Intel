/*
 * Sample code implementing the random number set up 
 * and generation for Monte Carlo Method, using MKL/VSL
 * random number generation and compiler vectorisation, 
 * plus OpenMP parallelisation
 *
 * author: Matthew Dixon and Xin Zhang
 *
 * Date: Nov 30 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>
#include <mkl_vsl.h>
#include <memory.h>
#include <omp.h>

#define VECTOR_LENGTH 4
#define ROUND_UP(x, y) ((((x) + (y) - 1) / (y)) * y)

/* Every thread has its own RNG set up */
VSLStreamStatePtr stream;
double *variates;
#pragma omp threadprivate(stream, variates)

void generation(double,double,int,int,int);

int main(int argc, char **argv)
{
  double T=1.0, sigma=0.2, dt, start, end;
  int M, N, N2, N3, tid, num_t, result;
  long long skip;

  M  = 20;      /* number of timesteps */
  N  = 5120000; /* total number of MC samples */
  N2 = 512;     /* number of MC samples generated in a group */

  dt  = T / (double)M;

start = omp_get_wtime();

#pragma omp parallel private(num_t,tid,N3)		\
                     shared(sigma,dt,M,N,N2)	\
                     reduction(+:result)
 {
    num_t = omp_get_num_threads();
    tid   = omp_get_thread_num();

    printf("tid=%d, creating RNG generator and allocating memory \n",tid); 

    /* create RNG, then give each thread a unique skipahead */
    vslNewStream(&stream, VSL_BRNG_MRG32K3A,1337);
    skip = ((long long) (tid+1)) << 48;
    vslSkipAheadStream(stream,skip);

    variates = (double *)malloc(M*N2*sizeof(double));

    N3 = ROUND_UP(((tid+1)*N)/num_t,VECTOR_LENGTH)
       - ROUND_UP(( tid   *N)/num_t,VECTOR_LENGTH);

    generation(sigma,dt, M,N2, N3);
	result = result + tid;
  }

end = omp_get_wtime();

printf("Time is: %f\n",(end-start));

/* delete generator and storage */
#pragma omp parallel 
  {
    vslDeleteStream(&stream);
    free(variates);
  }
}

void generation(double sigma,double dt, int M, int N2, int N) {
  int n0;
  /* loop over all paths in groups of size N2 */
  for (n0=0; n0<N; n0+=N2) {
    /* may have to reduce size of final group */
    if (N2>N-n0) N2 = N-n0;

    /* generate required random numbers for this group */
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2,stream,M*N2,variates,0,sqrt(dt));
  }
}
