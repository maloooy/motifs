SRCS = main.cpp motiv.cpp io.cpp scan.cpp
CXXFLAGS = -std=c++17 -O3 -march=native -flto -fopenmp

all:
	g++ $(CXXFLAGS) $(SRCS) -o main

mpi:
	mpic++ $(CXXFLAGS) -DUSE_MPI $(SRCS) -o main