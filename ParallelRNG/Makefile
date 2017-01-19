LIB=-mkl -qopenmp
CFLAGS=-O3 -qopt-report=5 -qopt-report-phase:vec,openmp
CC=icc

KNL:
	$(CC) -xMIC-AVX512 $(CFLAGS) $(LIB) ParallelRNG.c -o ParallelRNG.MIC

KNC:
	$(CC) $(CFLAGS) $(LIB) ParallelRNG.c -o ParallelRNG.MIC

X86:
	$(CC) $(CFLAGS) $(LIB) ParallelRNG.c -o ParallelRNG.X86
