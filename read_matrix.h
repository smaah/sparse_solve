
#include<iostream>
#pragma once

class matrix{

public:
  int m,n,nnz,myrank,size,maxDD; // GLOBAL row m ,column n, number of nonzeros nnz
  int *rowptr, *colidx,*rowidx,*colptr; // GLOBAL row pointer and column index
  int loc_m;   // local total row
  int *loc_m_size; // first row in each processor
  int *Ndpdnc; // number of dependency per row
  //int **Mapdpdnc; // dependency map for recieve
  int *MapSend; // dependency map for send
  int *SendTo;  // number of rows to send
  int *rowrank; // rank of the row
  int *x_done; // Array to trace of x found
  double *A;  // the matrix
  double *x;  // the unknown
  double *B;  // r.h.s
  char* file,*file2; // file name
  //int *dummy;
  matrix(char* fileIn,char* fileIn2){
    file = fileIn;
    file2 = fileIn2;
    m = 0;
    n = 0;
    nnz = 0;
    myrank = 0;
    size = 0;
    maxDD = 0;
  }
  ~matrix(){
    if(myrank==0) {
      std::cout<<"Deallocate Memory"<<std::endl;
    }
    delete[] rowptr;
    delete[] colidx;
    delete[] A;
    delete[] B;
    delete[] x;
    delete[] Ndpdnc;
    delete[] loc_m_size;
    delete[] SendTo;
    delete[] rowrank;
    delete[] x_done;
    
    delete [] MapSend;

  }
  void read_matrix();
  void read_matrix_ccs();
};


