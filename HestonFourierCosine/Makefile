LIB=-mkl -qopenmp
UT=utilities
CFLAGS=-O3 -I$(UT) -qopt-report=5 -qopt-report-phase:vec,openmp
CC=icc
SC=de4_0.cpp HestonFourierCosine.cpp utilities/dtime.cpp utilities/eval.cpp utilities/sort.cpp

KNL:
	$(CC) -xMIC-AVX512 $(CFLAGS) $(LIB) $(SC) -o HestonCalibration.MIC

KNC:
	$(CC) $(CFLAGS) $(LIB) $(SC) -o HestonCalibration.MIC

X86:
	$(CC) $(CFLAGS) $(LIB) $(SC) -o HestonCalibration.X86 
