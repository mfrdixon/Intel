LIBS=-mkl -qopenmp -lnlopt
UT=../utilities
INC=-I$(UT) -I$(UT)/nlopt/include
LIBPATH=$(UT)/nlopt/lib
CFLAGS=-O3 -qopt-report=5 -qopt-report-phase:vec,openmp
CC=icc
SC= main.cpp HestonFourierCosine.cpp $(UT)/de4_0.cpp $(UT)/eval.cpp $(UT)/sort.cpp

KNL:
	$(CC) -xMIC-AVX512 $(CFLAGS) $(SC) $(INC) -L$(LIBPATH) $(LIBS) -o HestonCalibration.MIC

KNC:
	$(CC) $(CFLAGS) $(SC) $(INC) -L$(LIBPATH) $(LIBS) -o HestonCalibration.MIC

X86:
	$(CC) $(CFLAGS) $(SC) $(INC) -L$(LIBPATH) $(LIBS) -o HestonCalibration.X86 
