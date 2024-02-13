#mpicxx -o sparse_solve sparse_solve.cc read_matrix.cc find_x.cc -Wall -Wextra -g
#mpicxx -o sparse_solve sparse_solve.cc read_matrix.cc find_x.cc -g
#mpicxx -O0  -fbacktrace -fmax-errors=10 -fbounds-check -o sparse_solve sparse_solve.cc read_matrix.cc find_x.cc -g
mpicxx -O3 -o sparse_solve sparse_solve.cc read_matrix.cc find_x.cc
#mpicxx -fopenmp -O3 -o sparse_solve sparse_solve.cc read_matrix.cc find_x.cc -g

echo -e "Compilation Done"

#mpiexec -np 4 ./sparse_solve nv_example.rcm nv_example.ccm
#mpiexec -np 4 ./sparse_solve nv_example.rcm
#mpiexec -np 8 ./sparse_solve matrixL10000.rcm matrixL10000.ccm
#mpiexec -np 16 ./sparse_solve matrixL10000.rcm
#mpiexec -np 1 ./sparse_solve TestMat2.rcm
#mpiexec -np 2 ./sparse_solve Testmatrix100.rcm
#mpiexec -np 4 ./sparse_solve ASmatrix10.rcm ASmatrix10.ccm
mpiexec -np 4 ./sparse_solve ASmatrix10.rcm 
#mpiexec -np 32 ./sparse_solve ASmatrix1000000.rcm
#mpiexec -np 4 ./sparse_solve Freelance1_test.rcm
#mpiexec -np 32 ./sparse_solve nlpkkt120_corr.rcm
#mpiexec -np 1 ./sparse_solve Geo_1438_corr.rcm
#mpiexec -np 8 ./sparse_solve Freescale1_corr.rcm
#mpiexec -np 4 ./sparse_solve Mars2020_FullAxi_90K_1_1_0_ordered.rcm
#mpiexec -np 1 ./sparse_solve EES_Full_3.1M_1_1_1_ordered.rcm
