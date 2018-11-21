# configure gcc with mpi
./configure CC=mpicc CXX=mpic++ CXXFLAGS="-m64 -O3 -g -Wall -std=c++11 -Wno-unused-arguments -Wno-unused-local-typedef" LIBS="-lboost_mpi -lboost_serialization -lboost_thread -lboost_system -lmpi"
