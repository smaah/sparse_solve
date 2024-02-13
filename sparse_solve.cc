
#include<iostream>
#include<fstream>
#include"read_matrix.h"
#include"find_x.h"
#include "mpi.h"
#include <chrono>
#include <ctime>
#include <string>
#include <cmath>

using namespace std::chrono;

#define debug 0
#define output 0
#define timing 0

int main(int argc, char *argv[]){
  
  int myrank,dummy,i,size,j;
  double time_elapsed_ms_dep(0.0),time_elapsed_ms_for(0.0),time_elapsed_ms_tot(0.0);
  char* fileIn;
  std::string filename = "OutL10000_vali_binary_opt_sm.dat";
  std::ofstream outdata;

  
  MPI_Init( &argc, &argv ); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);




  if(argc>2){
    fileIn = argv[2];
  }else{
    fileIn = new char[1];
    fileIn[0] = 'A';
  }



  matrix M(argv[1],fileIn);
  
  M.myrank = myrank;
  M.size   = size;
  dummy    = 0;
  
  if(myrank==0){
    
    if(argc>2){
      M.read_matrix();
      M.read_matrix_ccs();
    }else{
      M.read_matrix();
    }
  }

  MPI_Bcast(&M.m,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&M.n,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&M.nnz,1,MPI_INT,0,MPI_COMM_WORLD);
  
  
  if(myrank>0){
    M.rowptr   = new int[M.m+1];
    M.colidx   = new int[M.nnz];
    M.colptr   = new int[M.m+1];
    M.rowidx   = new int[M.nnz]; 
    M.A        = new double[M.nnz];
    M.x        = new double[M.m];
    M.B        = new double[M.m];
    M.Ndpdnc   = new int[M.m];
    M.loc_m_size = new int[M.size+1];
    M.rowrank    = new int[M.m];
    M.x_done     = new int[M.m];
    
    //initialization
    for(i=0;i<M.m;++i){
      M.rowptr[i]  = 0;
      M.colptr[i]  = 0;
      M.x[i]       = 0.0;
      M.B[i]       = 0.0;
      M.Ndpdnc[i]  = 0.0;
      M.rowrank[i] = 0.0;
      M.x_done[i]  = 0;
    }
    
    M.rowptr[M.m] = 0;
    M.colptr[M.m] = 0;
    
    for(i=0;i<M.nnz;++i){
      M.colidx[i] = 0;
      M.rowidx[i] = 0;
      M.A[i]      = 0.0;
    }
  }

  //std::cout<<"!!!Allocation done!!!"<<std::endl;
  
  MPI_Bcast(M.rowptr,M.m+1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(M.colidx,M.nnz,MPI_INT,0,MPI_COMM_WORLD);
  if(argc>2){
    MPI_Bcast(M.colptr,M.m+1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(M.rowidx,M.nnz,MPI_INT,0,MPI_COMM_WORLD);
  }
  MPI_Bcast(M.A,M.nnz,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(M.B,M.m,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(M.x,M.m,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
  
  for(i=0;i<M.m;++i){
    
    M.x_done[i] = 0;
    M.x[i]      = 0.0;
    
  }
  
#if(debug==1)
  
  if(myrank>0){
    //
    for(i=0;i<M.m+1;++i){
      std::cout<<"row "<<i<<" "<<M.rowptr[i]<<"\n";
    }

    for(i=0;i<M.nnz;++i){
      std::cout<<"col "<<i<<" "<<M.colidx[i]<<"\n";
    }

    for(i=0;i<M.nnz;++i){
      std::cout<<"val "<<i<<" "<<M.A[i]<<"\n";
    }   
  }
  
#endif


  std::clock_t c_start = std::clock();
  // store starting rows in each processor
  M.loc_m_size[0] = 0;
  
  for (i=1;i<M.size;++i){
    
    if(M.m%M.size!=0){
      if(i==1){
	M.loc_m_size[i] = M.loc_m_size[i-1]+floor(M.m/M.size)+int(M.m%M.size);
      }else{
	M.loc_m_size[i] = M.loc_m_size[i-1]+floor(M.m/M.size);
      }

    }else{
      
      M.loc_m_size[i] = M.loc_m_size[i-1]+floor(M.m/M.size);

      
    }
    
  }
  
  M.loc_m_size[M.size] = M.m;
  

  for(j=0;j<M.size;++j){
    
    for(i=0;i<M.m;++i){
      
      if(i>=M.loc_m_size[j] && i <M.loc_m_size[j+1]){

	M.rowrank[i] = j;

      }       
    }
  }
 

  depepndency_map_rank(M);
  

  std::clock_t c_end1 = std::clock();
  double time_elapsed_ms1 = double(c_end1-c_start) / CLOCKS_PER_SEC;
  std::clock_t c_start2 = std::clock();
  
  find_x(M);
  
  std::clock_t c_end = std::clock();
  double time_elapsed_ms =  double(c_end-c_start) / CLOCKS_PER_SEC;

  
  double time_elapsed_ms2 =  double(c_end-c_start2) / CLOCKS_PER_SEC;


  time_elapsed_ms1=time_elapsed_ms1/M.size;
  time_elapsed_ms=time_elapsed_ms/M.size;
  time_elapsed_ms2=time_elapsed_ms2/M.size;
  
  MPI_Reduce(&time_elapsed_ms1,&time_elapsed_ms_dep,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&time_elapsed_ms,&time_elapsed_ms_tot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&time_elapsed_ms2,&time_elapsed_ms_for,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  if(M.myrank==0){
    std::cout<<"Dependnecy Map took =" <<time_elapsed_ms_dep<<" s\n";
    std::cout<<"Forward Substitution took =" <<time_elapsed_ms_for<<" s\n";
    std::cout<<"Total time took =" <<time_elapsed_ms_tot<<" s\n";
  }

#if(output==1)
  
  if(M.myrank==0){
    outdata.open(filename);  
    for(i=0;i<M.m;++i){
      if(i>=M.loc_m_size[M.myrank] && i<M.loc_m_size[M.myrank+1]){
	outdata<<M.myrank<<" "<<i+1<<" "<<M.x[i]<<std::endl;
      }
    }

    if(M.size>1){
      MPI_Send(&dummy,1,MPI_INT,M.myrank+1,0,MPI_COMM_WORLD);
    }
    
  }else{
    
    MPI_Recv(&dummy,1,MPI_INT,M.myrank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    outdata.open(filename,std::ios_base::app);
    for(i=0;i<M.m;++i){
      if(i>=M.loc_m_size[M.myrank] && i<M.loc_m_size[M.myrank+1]){
	outdata<<M.myrank<<" "<<i+1<<" "<<M.x[i]<<std::endl;
      }
    }
    
    if(M.myrank<M.size-1){
      
      MPI_Send(&dummy,1,MPI_INT,M.myrank+1,0,MPI_COMM_WORLD);
      
    }
  }
  
#endif
  
  MPI_Finalize();
  return(0); 	  
}



