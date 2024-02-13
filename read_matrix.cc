
#include<iostream>
#include<fstream>
#include<string>
#include"mpi.h"
#include"read_matrix.h"



void matrix::read_matrix()
{
  std::fstream rcm;
  rcm.open(file,std::ios::in);
  //
  rcm >> m >> n >> nnz;


  
  rowptr = new int[m+1];
  colidx = new int[nnz];
  A      = new double[nnz];
  x      = new double[m];
  B      = new double[m];
  Ndpdnc = new int[m];
  loc_m_size = new int[size+1];
  rowrank    = new int[m];
  x_done     = new int[m];


  for(int i=0;i<m;++i){
      rowptr[i]  = 0;
      x[i]       = 0.0;
      B[i]       = 0.0;
      Ndpdnc[i]  = 0.0;
      rowrank[i] = 0.0;
      x_done[i]  = 0;
    }
    
    rowptr[m] = 0;
    
    for(int i=0;i<nnz;++i){
      colidx[i] = 0;
      A[i]      = 0.0;
    }
  
  for(int i=0;i<m+1;++i){
    rcm >> rowptr[i];

    //initialization of arrays
    if(i<m) {
      x[i] = 0.0;
      B[i] = 1.0; // r.h.s is 1
      x_done[i] = 0;
    }
    
  }

  for(int i=0;i<nnz;++i){
    
    rcm >> colidx[i];
    
  }

  for(int i=0;i<nnz;++i){
    
    rcm >> A[i];
    
  }
 
 
}




void matrix::read_matrix_ccs()
{
  std::fstream rcm;
  rcm.open(file2,std::ios::in);
  //
  rcm >> m >> n >> nnz;

  
  
  colptr = new int[m+1];
  rowidx = new int[nnz];
  

  for(int i=0;i<m;++i){
      colptr[i]  = 0;      
    }
    
  colptr[m] = 0;
  
  for(int i=0;i<nnz;++i){
    rowidx[i] = 0;
  }
  
  for(int i=0;i<m+1;++i){
    rcm >> colptr[i];
    
  }

  for(int i=0;i<nnz;++i){
    
    rcm >> rowidx[i];
    
  }
 
}


