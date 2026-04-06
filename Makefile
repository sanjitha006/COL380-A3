# Makefile for Budgeted Maximum Clique — Sequential and MPI versions
CXX     = g++
MPICXX  = mpicxx
CFLAGS  = -O2 -std=c++17

all: compile_seq compile_mpi run_seq run_mpi

compile_seq: seq.cpp
	$(CXX) $(CFLAGS) -o clique_seq seq.cpp

compile_mpi: test.cpp
	$(MPICXX) $(CFLAGS) -o clique_mpi test.cpp

run_seq: clique_seq
	./clique_seq input1.txt output1_seq.txt

run_mpi: clique_mpi
	mpirun -n 4 ./clique_mpi input1.txt output1_mpi.txt

clean:
	rm -f clique_seq clique_mpi