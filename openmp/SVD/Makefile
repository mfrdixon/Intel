LIB=-mkl -qopenmp
CFLAGS=-O3 -qopt-report=5 -qopt-report-phase:vec,openmp
CC=icc

KNL:
	$(CC) -xMIC-AVX512 $(CFLAGS) $(LIB) SVD.c -o SVD.MIC

KNC:
	$(CC) $(CFLAGS) $(LIB) SVD.c -o SVD.MIC

X86:
	$(CC) $(CFLAGS) $(LIB) SVD.c -o SVD.X86 
