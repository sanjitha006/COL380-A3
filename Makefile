# Makefile for Budgeted Maximum Clique — Sequential and MPI versions
CXX     = g++
MPICXX  = mpicxx
CFLAGS  = -O3 -std=c++17

all: compile_seq compile_mpi run_seq run_mpi

compile_seq: seq.cpp
	$(CXX) $(CFLAGS) -o clique_seq seq.cpp

compile_mpi: test.cpp
	$(MPICXX) $(CFLAGS) -o clique_mpi main.cpp

run_seq: clique_seq
	time ./clique_seq tainput.txt taoutput_seq.txt

run_mpi: clique_mpi
	time mpirun -x UCX_LOG_LEVEL=error -n 4 ./clique_mpi tainput.txt taoutput_mpi.txt
clean:
	rm -f clique_seq clique_mpi