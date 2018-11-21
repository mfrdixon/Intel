# configure clang on mac
./configure --with-boost-include=/Users/peter/boost_1_59_0 --with-boost-lib=/Users/peter/boost_1_59_0/stage/lib CC=mpicc CXX=mpic++ CXXFLAGS="-m64 -O3 -g -Wall -std=c++11 -Qunused-arguments -Wno-unused-local-typedef" LIBS="-lboost_mpi -lboost_serialization -lboost_thread -lboost_system -lmpi"
