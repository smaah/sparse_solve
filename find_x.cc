
#include"read_matrix.h"
#include <cmath>
#include "mpi.h"
#include <ctime>
#include <algorithm>
#include <vector>
#include <omp.h>

#define debug 0
#define output 1
#define timing 0

int find_x(matrix &M){

  int i,j,k,x_done_min,dummy,count,errcode,tag;
  int rank,total,v_s,sum_v,v_count;
  double *sum = new double[M.m]; 
  bool    *x_sent = new bool[M.m];
  int    *recv_count = new int[M.size];
  bool *trace_sum = new bool[M.rowptr[M.loc_m_size[M.myrank+1]]-M.rowptr[M.loc_m_size[M.myrank]]];
  int    *recv_count_v = new int[M.size];
  int    *v_size = new int[M.size];
  int    *v_displs = new int[M.size];
  int    *v_glob_displs = new int[M.size];
  int srank,rrank,local_size;
  double dummyD;
  int* v = new int[M.loc_m_size[M.myrank+1] - M.loc_m_size[M.myrank]];
  int* v_glob = new int[M.m];
  
  for(i=0;i<(M.rowptr[M.loc_m_size[M.myrank+1]]-M.rowptr[M.loc_m_size[M.myrank]]);++i){
    trace_sum[i] = false;
  }

  
  for(i=0;i<M.m;++i){
    sum[i] = 0.0;
    x_sent[i] = false;
    if(i<M.loc_m_size[M.myrank+1] - M.loc_m_size[M.myrank]){
      v[i] = 0;
    }
    v_glob[i]=0;
  }


#if(debug==1)
  
  std::cout<<"myrank ="<<M.myrank<<" m ="<<M.m<<std::endl;
  std::cout<<"myrank ="<<M.myrank<<" nnz ="<<M.nnz<<std::endl;
  
#endif

  
  for(i=0;i<M.size;++i){
    recv_count_v[i] = 1;
    v_size[i] = 0;
    v_displs[i] = i;
  }
  
  for(i=0;i<M.size;++i){

    if(M.m%M.size!=0){
	if(i==0){
	  recv_count[i]=floor(M.m/M.size)+int(M.m%M.size);	  
	}else{
	  recv_count[i]=floor(M.m/M.size);	  
	}
      }else{
	recv_count[i]=floor(M.m/M.size);	  
      }     
  }

  int *frac_x_done = new int[recv_count[M.myrank]];
  
  v_count = 0;
  for(i=M.loc_m_size[M.myrank];i<M.loc_m_size[M.myrank+1];++i){
    
    if(M.Ndpdnc[i]==0){ 
      M.x[i]=M.B[i]/M.A[M.rowptr[i]];
      M.x_done[i] = 1;
      v[v_count]=i;
      ++v_count;
    }
    
    frac_x_done[i-M.loc_m_size[M.myrank]] = M.x_done[i];
    
  }
  
  MPI_Allgatherv(&v_count,1,MPI_INT,v_size,recv_count_v,v_displs,MPI_INT,MPI_COMM_WORLD);

  v_glob_displs[0] = 0;
  sum_v = 0;
  for(i=0;i<M.size;++i){
    sum_v = sum_v+v_size[i];
    if(i+1<M.size){
      v_glob_displs[i+1]=v_glob_displs[i]+v_size[i];
    }
  }

  MPI_Allgatherv(v,v_size[M.myrank],MPI_INT,v_glob,v_size,v_glob_displs,MPI_INT,MPI_COMM_WORLD);
  

  for(i=0;i<sum_v;++i){
    M.x_done[v_glob[i]] = 1;
  }


#if(timing==1)
  std::clock_t c_start = std::clock();
#endif

  x_done_min = 1;
  for(i=0;i<M.m;++i){
    if(M.x_done[i]<x_done_min){
      
      x_done_min = M.x_done[i];
          
    }    
    if(M.x_done[i]==1 && x_sent[i] == false && M.SendTo[i]>0){
      for(j=0;j<M.SendTo[i];++j){
	local_size = 1;

	  while(local_size<M.SendTo[i]){
	    
	    if(j<local_size){
	      if((j+local_size)<M.SendTo[i] && M.myrank == M.MapSend[i*M.size+j]){
		MPI_Send(&M.x[i],1,MPI_DOUBLE,M.MapSend[i*M.size+j+local_size],0,MPI_COMM_WORLD);
	      }
	    }
	    
	    if(j >= local_size && j < (local_size*2) && M.myrank == M.MapSend[i*M.size+j]){
	      if((j-local_size)<M.SendTo[i]){		
		MPI_Recv(&M.x[i],1,MPI_DOUBLE,M.MapSend[i*M.size+j-local_size],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      }
	    }
	  
	    local_size = local_size*2;
	  
	  }	  
      }
      x_sent[i] =true;
    }
  }
 
  
#if(timing==1)
  std::clock_t c_end1 = std::clock();
  double time_elapsed_ms1 = double(c_end1-c_start) / CLOCKS_PER_SEC;
  std::cout << "Communication time used: " << time_elapsed_ms1 << " s\n";
#endif
  

#if(timing==1)
  std::clock_t c_start1 = std::clock();
  
#endif

#if(timing==1)
  std::clock_t time_s,time_f;
  double tot_time;
  tot_time = 0.0;
#endif
  

  while(x_done_min==0){

    v_count = 0;    
    for(i=M.loc_m_size[M.myrank];i<M.loc_m_size[M.myrank+1];++i){
      
      if(M.Ndpdnc[i]>0 && M.x_done[i] ==0 ){
	
	for(j=M.rowptr[i];j<=M.rowptr[i+1]-2;++j){
	  if(M.x_done[M.colidx[j]]==1 && trace_sum[j-M.rowptr[M.loc_m_size[M.myrank]]] == false){
	    
	    sum[i] =sum[i]+M.A[j]*M.x[M.colidx[j]];
	    M.Ndpdnc[i] = M.Ndpdnc[i]-1;
	    trace_sum[j-M.rowptr[M.loc_m_size[M.myrank]]] = true;
	    if(M.Ndpdnc[i]==0 && M.x_done[i]==0){
	      
	      M.x[i]=(M.B[i]-sum[i])/M.A[M.rowptr[i+1]-1];
	      M.x_done[i] = 1;

	      v[v_count] = i;
	      ++v_count;
	      frac_x_done[i-M.loc_m_size[M.myrank]] = 1;

	    } 	    
	  }
	}
      }
    }

    
    MPI_Allgatherv(&v_count,1,MPI_INT,v_size,recv_count_v,v_displs,MPI_INT,MPI_COMM_WORLD);
    
    v_glob_displs[0] = 0;
    sum_v = 0;
    for(i=0;i<M.size;++i){
      sum_v = sum_v+v_size[i];
      if(i+1<M.size){
    	v_glob_displs[i+1]=v_glob_displs[i]+v_size[i];
      }
    }
    

    
    MPI_Allgatherv(v,v_size[M.myrank],MPI_INT,v_glob,v_size,v_glob_displs,MPI_INT,MPI_COMM_WORLD);
    
    
    for(i=0;i<sum_v;++i){
      M.x_done[v_glob[i]] = 1;
    }
    
    

#if(timing==1)
    time_s = std::clock();
#endif
    
    

    x_done_min = 1;
    for(i=0;i<M.m;++i){
      if(M.x_done[i]<x_done_min){
	
      x_done_min = M.x_done[i];
      
      
      }    
      if(M.x_done[i]==1 && x_sent[i] == false && M.SendTo[i]>0){
	for(j=0;j<M.SendTo[i];++j){
	  local_size = 1;
	  
	  while(local_size<M.SendTo[i]){
	    
	    if(j<local_size){
	      if((j+local_size)<M.SendTo[i] && M.myrank == M.MapSend[i*M.size+j]){
		MPI_Send(&M.x[i],1,MPI_DOUBLE,M.MapSend[i*M.size+j+local_size],0,MPI_COMM_WORLD);
	      }
	    }
	    
	    if(j >= local_size && j < (local_size*2) && M.myrank == M.MapSend[i*M.size+j]){
	      if((j-local_size)<M.SendTo[i]){
		MPI_Recv(&M.x[i],1,MPI_DOUBLE,M.MapSend[i*M.size+j-local_size],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      }
	    }
	    
	    local_size = local_size*2;
	    
	  }

      }
	x_sent[i] =true;
      }
    }
    

#if(timing==1)
    time_f = std::clock();
    tot_time=tot_time+double(time_f-time_s) / CLOCKS_PER_SEC;
#endif
    
  } // while loop ends
  
#if(timing==1)
  std::cout<<"MPI Communication took = "<<tot_time<<std::endl;
#endif
  
#if(timing==1)
  std::clock_t c_end = std::clock();
  double time_elapsed_ms2 = double(c_end-c_start1) / CLOCKS_PER_SEC;
  std::cout << "While loop time used: " << time_elapsed_ms2 << " s\n";
#endif
  
#if(output==1)
  
  for(i=0;i<M.m;++i){
    if(i>=M.loc_m_size[M.myrank] && i<M.loc_m_size[M.myrank+1]){
      std::cout<<"rank ="<<M.myrank<<" x"<<i+1<<" =  "<<M.x[i]<<std::endl;
    }
  }
  
#endif
  
  delete[] sum;
  delete[] frac_x_done;
  delete[] trace_sum;
  delete[] x_sent;
  delete[] recv_count;
  delete[] recv_count_v;
  delete[] v_size;
  delete[] v_displs;
  delete[] v_glob_displs;
  delete[] v;
  delete[] v_glob;
    
  return 0;
  
}











int dependency_map(matrix &M){

  int i,maxD,n_col,maxDD,j,count,k;
  int *idx=new int[M.m];
  int *coldep = new int[M.m];

    
  //maxD = 0;
  for (i=0;i<M.m;++i){
    
    coldep[i] = 0;
    
    if(M.rowptr[i+1]-M.rowptr[i] == 1){
      M.Ndpdnc[i] = 0;
    }else{      
      M.Ndpdnc[i] = M.rowptr[i+1]-M.rowptr[i]-1;
    }
    
  }
  


#if(timing==1)
  std::clock_t c_startt = std::clock();
#endif
  
  maxDD = 0;

  
  for(i=0;i<M.nnz;++i){
    
    ++coldep[M.colidx[i]];
    
  }
  
  for(i=0;i<M.m;++i){
    if(coldep[i]>maxDD){
      maxDD = coldep[i];
    }
  }

  M.maxDD = maxDD;
  
#if(timing==1)
  std::clock_t c_endd = std::clock();
  double time_elapsed_msd = 1000.0 * (c_endd-c_startt) / CLOCKS_PER_SEC;
  std::cout << "MaxDD time used: " << time_elapsed_msd / 1000.0 << " s\n";
#endif  
  
  // Allocate arrays
#if(timing==1)
  std::clock_t c_start2 = std::clock();
#endif

  
  M.MapSend  = new int[M.m*maxDD];
  M.SendTo   = new int[M.m];
  
  for(i=0;i<M.m;++i){
    M.SendTo[i] = 0;
    idx[i] = 0;
 
  }

#if(timing==1)
  std::clock_t c_end = std::clock();
  double time_elapsed_ms = 1000.0 * (c_end-c_start2) / CLOCKS_PER_SEC;
  std::cout << "Allocation time used: " << time_elapsed_ms / 1000.0 << " s\n";
#endif
 
#if(timing==1)
  std::clock_t c_start4 = std::clock();
#endif
  
  
 
  for(i=0;i<M.m;++i){
    if(M.Ndpdnc[i]>0){
      for(j = M.rowptr[i];j<=M.rowptr[i+1]-2;++j){
 
	M.MapSend[M.colidx[j]*maxDD+idx[M.colidx[j]]] = i;
  	M.SendTo[M.colidx[j]] = M.SendTo[M.colidx[j]]+1; 
 
  	++idx[M.colidx[j]];
      }
    }
  }
  
#if(timing==1)
  std::clock_t c_end3 = std::clock();
  double time_elapsed_ms3 = 1000.0 * (c_end3-c_start4) / CLOCKS_PER_SEC;
  std::cout << "Mapping time used: " << time_elapsed_ms3 / 1000.0 << " s\n";
#endif
  
   
  delete[] idx;
  delete[] coldep;
  
  return 0;
}





int dependency_map_ccs(matrix &M){

  int i,maxDD,j,idx;
  int* MapSendLocal;
  int* MapSendCount = new int[(M.loc_m_size[M.myrank+1]-M.loc_m_size[M.myrank])];
  int* recv_count = new int[M.size];
  int* recv_count2 = new int[M.size];
  int* displs     = new int[M.size];\
  int* displs2     = new int[M.size];

  for (i=0;i<M.m;++i){
      
    if(M.rowptr[i+1]-M.rowptr[i] == 1){
      M.Ndpdnc[i] = 0;
    }else{      
      M.Ndpdnc[i] = M.rowptr[i+1]-M.rowptr[i]-1;
    }    
  }
  
  
  maxDD = 0;
  
  for(i=M.loc_m_size[M.myrank];i<M.loc_m_size[M.myrank+1];++i){
    MapSendCount[i-M.loc_m_size[M.myrank]] = 0;
    maxDD = std::max(maxDD,M.colptr[i+1]-M.colptr[i]-1);         
  }

  
  MPI_Allreduce(&maxDD,&M.maxDD,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  M.MapSend  = new int[M.m*M.maxDD];
  M.SendTo   = new int[M.m];

  displs[0] = 0;
  recv_count[0] = (M.loc_m_size[1]-M.loc_m_size[0])*M.maxDD;

  displs2[0] = 0;
  recv_count2[0] = (M.loc_m_size[1]-M.loc_m_size[0]);

  
  for(i=1;i<M.size;++i){
    
    recv_count[i] = (M.loc_m_size[i+1]-M.loc_m_size[i])*M.maxDD;
    recv_count2[i] = (M.loc_m_size[i+1]-M.loc_m_size[i]);

    displs[i] = displs[i-1]+recv_count[i-1];
    displs2[i] = displs2[i-1]+recv_count2[i-1];
  }


  MapSendLocal = new int[(M.loc_m_size[M.myrank+1]-M.loc_m_size[M.myrank])*M.maxDD];


  
  idx=0;
  for(i=M.loc_m_size[M.myrank];i<M.loc_m_size[M.myrank+1];++i){
    idx =0;
    for(j = M.colptr[i]+1;j<M.colptr[i+1];++j){
      
      MapSendLocal[(i-M.loc_m_size[M.myrank])*M.maxDD+idx]= M.rowidx[j];
      MapSendCount[(i-M.loc_m_size[M.myrank])] = MapSendCount[(i-M.loc_m_size[M.myrank])]+1;

      ++idx;
    }
  }

  MPI_Allgatherv(MapSendLocal,recv_count[M.myrank],MPI_INT,M.MapSend,recv_count,displs,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgatherv(MapSendCount,recv_count2[M.myrank],MPI_INT,M.SendTo,recv_count2,displs2,MPI_INT,MPI_COMM_WORLD);



  delete[] MapSendLocal;
  delete[] MapSendCount;
  delete[] recv_count;
  delete[] recv_count2;
  delete[] displs;
  delete[] displs2;

  return 0;
}


int depepndency_map_rank(matrix &M){

  int i,j;
  int *idx=new int[M.m];

  M.MapSend  = new int[M.m*M.size];
  M.SendTo   = new int[M.m];

  for (i=0;i<M.m;++i){
    idx[i] = 1;
    M.SendTo[i] = 1;
    if(M.rowptr[i+1]-M.rowptr[i] == 1){
      M.Ndpdnc[i] = 0;
    }else{      
      M.Ndpdnc[i] = M.rowptr[i+1]-M.rowptr[i]-1;
    }    
    M.MapSend[i*M.size] = M.rowrank[i];
  }
  

  for(i=0;i<M.m;++i){
    for(j = M.rowptr[i];j<=M.rowptr[i+1]-2;++j){
      if(idx[M.colidx[j]]==0){
	M.MapSend[M.colidx[j]*M.size+idx[M.colidx[j]]] = M.rowrank[i];
	M.SendTo[M.colidx[j]] = M.SendTo[M.colidx[j]]+1;
	++idx[M.colidx[j]];
	
      }else if(M.rowrank[i]>M.MapSend[M.colidx[j]*M.size+idx[M.colidx[j]]-1]){
	M.MapSend[M.colidx[j]*M.size+idx[M.colidx[j]]] = M.rowrank[i];
	M.SendTo[M.colidx[j]] = M.SendTo[M.colidx[j]]+1;
	++idx[M.colidx[j]];
	
      }
    }
  }

  
  delete[] idx;
  return 0;
}
