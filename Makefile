CXX=mpicxx

OBJS = sparse_solve.o read_matrix.o find_x.o

all: sparse_solve

sparse_solve: $(OBJS)
	$(CXX) -O0 -g -o sparse_solve $(OBJS)  #-tracemode projections 

clean:
	rm -f *.o sparse_solve

sparse_solve.o: sparse_solve.cc 
	$(CXX) -c sparse_solve.cc read_matrix.cc find_x.cc

test: all
	mpiexec -n 2 ./sparse_solve nv_example.rcm
