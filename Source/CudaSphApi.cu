/*
 <DUALSPHYSICS>  Copyright (C) 2012 by Jose M. Dominguez, Dr. Alejandro Crespo, Prof. M. Gomez Gesteira, Anxo Barreiro, Dr. Benedict Rogers.

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file CudaSphApi.cu \brief Implements all the functions used on GPU execution.

#include "CudaSphApi.h"
#include <string>
#include <float.h>
#include <thrust/device_vector.h>  //cu40
#include <thrust/sort.h>           //cu40
#include <cuda.h>  


__constant__ StDeviceCte CTE;

#define CheckErrorCuda(text)  __CheckErrorCuda(text,__FILE__,__LINE__)
void __CheckErrorCuda(const char *test,const char *file,const int line);
dim3 CsGetGridSize(unsigned n,unsigned blocksize);

template <unsigned blockSize> __device__ void WarpReduMaxF(volatile float* sax,unsigned tid);
template <unsigned blockSize> __device__ void WarpReduMaxU(volatile unsigned* sax,unsigned tid);

float ReduMaxF(unsigned ndata,float* data,float* resu);
unsigned ReduMaxU(unsigned ndata,unsigned* data,unsigned* resu);

unsigned CsGetCountRhopOut(StDeviceContext *dc,float rhopoutmin,float rhopoutmax);


#include "CudaSphNL.cu"
#include "CudaSphFC.cu"
#include "CudaSphSU.cu"


//==============================================================================
/// Initialisation of CUDA device. 
//==============================================================================
int CsInitCuda(int gpuid,int *devdef,char *devname,StDeviceContext *dc){
  *devdef=-1;
  devname[0]='\0';
  bool err=false;
  printf("[Initialising of CUDA device]\n");
  int devcount;
  cudaGetDeviceCount(&devcount);
  CheckErrorCuda("Failure getting devices info.");
  for(int dev=0;dev<devcount;dev++){
    cudaDeviceProp devp;
    cudaGetDeviceProperties(&devp,dev);
    printf("\nDevice %d: \"%s\"\n",dev,devp.name);
    printf("  Compute capbility:         %d.%d\n",devp.major,devp.minor);
    int corebymp=(devp.major==1? 8: (devp.major==2? 32: -1));
    printf("  Multiprocessors:           %d (%d cores)\n",devp.multiProcessorCount,devp.multiProcessorCount*corebymp);
    printf("  Global memory :            %d MB\n",int(devp.totalGlobalMem/(1024*1024)));
    printf("  Clock rate:                %.2f GHz\n",devp.clockRate*1e-6f);
    #if CUDART_VERSION >= 2020
    printf("  Run time limit on kernels: %s\n",(devp.kernelExecTimeoutEnabled? "Yes": "No"));
    #endif
  }
  printf("\n");
  if(devcount){
    cudaDeviceProp devp;
    int dev;
    if(gpuid<0){
      cudaGetDevice(&dev);
      cudaGetDeviceProperties(&devp,dev);
      *devdef=dev;
      int len=int(strlen(devp.name)); if(len>=64)len=63;
      for(int c=0;c<len;c++)devname[c]=devp.name[c]; devname[len]='\0';

      dc->mglobal=devp.totalGlobalMem;
      dc->mshared=int(devp.sharedMemPerBlock);
      dc->compute=devp.major*10+devp.minor;
    }
    else{ 
      dev=gpuid;
      cudaSetDevice(dev);
      cudaGetDeviceProperties(&devp,dev);
      int len=int(strlen(devp.name)); if(len>=64)len=63;
      for(int c=0;c<len;c++)devname[c]=devp.name[c]; devname[len]='\0';
      dc->mglobal=devp.totalGlobalMem;
      dc->mshared=int(devp.sharedMemPerBlock);
      dc->compute=devp.major*10+devp.minor;
    }
  } 
  else{ //printf("\nERROR: There is no CUDA devices available.\n");
    err=true;
  }
  return(err);
}

//==============================================================================
/// Initialisation of variables and pointers to NULL.
//==============================================================================
void CsInit(StDeviceContext *dc){
  memset(dc,0,sizeof(StDeviceContext));
  TmgCreation(dc->timers,false);
}

//==============================================================================
/// Releases memory and reinitialises values.
//==============================================================================
void CsReset(StDeviceContext *dc){
  //-Data arrays (device).
  cudaFree(dc->idp);
  cudaFree(dc->pos);
  cudaFree(dc->vel);
  cudaFree(dc->rhop);
  cudaFree(dc->ridpmv);
  cudaFree(dc->ridpft);
  cudaFree(dc->idsortb);    cudaFree(dc->idsortf);
  cudaFree(dc->cellb);      cudaFree(dc->cellf);
  cudaFree(dc->ar);
  cudaFree(dc->ace);
  cudaFree(dc->velxcor);
  cudaFree(dc->viscdt);
  cudaFree(dc->csound);
  cudaFree(dc->rhopm1);
  cudaFree(dc->velm1);   cudaFree(dc->velm1_2);
  cudaFree(dc->velnew);
  cudaFree(dc->rhopnew);
  cudaFree(dc->pospre);  cudaFree(dc->pospre_2);
  cudaFree(dc->velpre);  cudaFree(dc->velpre_2);
  cudaFree(dc->rhoppre);
  cudaFree(dc->tau);     cudaFree(dc->csph);
  cudaFree(dc->matkgc);
  cudaFree(dc->ftdist);
  cudaFree(dc->face);    cudaFree(dc->fomegavel);
  cudaFree(dc->fdrhop);
  cudaFree(dc->reduaux);
  cudaFree(dc->cellbegb); cudaFree(dc->cellbegf);
  cudaFree(dc->cellnvb); cudaFree(dc->cellnvf);
  cudaFree(dc->pospres); cudaFree(dc->velrhop);
  //-Excluded particles Out (host).
  delete[] dc->outidp;
  delete[] dc->outpos;
  delete[] dc->outvel;
  delete[] dc->outrhop;
  //-Destroys used cudaEvent to measure times.
  TmgDestruction(dc->timers);
  //-DeviceContext to NULL.
  memset(dc,0,sizeof(StDeviceContext));
  CheckErrorCuda("Failed releasing GPU memory.");
}

//==============================================================================
/// Memory allocation.
//==============================================================================
unsigned AllocMemory(unsigned count,unsigned ** pt1,unsigned ** pt2=NULL,unsigned ** pt3=NULL,unsigned ** pt4=NULL){
  unsigned m=0;
  if(pt1){ unsigned s=sizeof(unsigned)*count; cudaMalloc(pt1,s); m+=s; }
  if(pt2){ unsigned s=sizeof(unsigned)*count; cudaMalloc(pt2,s); m+=s; }
  if(pt3){ unsigned s=sizeof(unsigned)*count; cudaMalloc(pt3,s); m+=s; }
  if(pt4){ unsigned s=sizeof(unsigned)*count; cudaMalloc(pt4,s); m+=s; }
  return(m);
}
//==============================================================================
unsigned AllocMemory(unsigned count,float3 ** pt1,float3 ** pt2=NULL,float3 ** pt3=NULL,float3 ** pt4=NULL){
  unsigned m=0;
  if(pt1){ unsigned s=sizeof(float3)*count; cudaMalloc(pt1,s); m+=s; }
  if(pt2){ unsigned s=sizeof(float3)*count; cudaMalloc(pt2,s); m+=s; }
  if(pt3){ unsigned s=sizeof(float3)*count; cudaMalloc(pt3,s); m+=s; }
  if(pt4){ unsigned s=sizeof(float3)*count; cudaMalloc(pt4,s); m+=s; }
  return(m);
}
//==============================================================================
unsigned AllocMemory(unsigned count,float4 ** pt1,float4 ** pt2=NULL,float4 ** pt3=NULL,float4 ** pt4=NULL){
  unsigned m=0;
  if(pt1){ unsigned s=sizeof(float4)*count; cudaMalloc(pt1,s); m+=s; }
  if(pt2){ unsigned s=sizeof(float4)*count; cudaMalloc(pt2,s); m+=s; }
  if(pt3){ unsigned s=sizeof(float4)*count; cudaMalloc(pt3,s); m+=s; }
  if(pt4){ unsigned s=sizeof(float4)*count; cudaMalloc(pt4,s); m+=s; }
  return(m);
}
//==============================================================================
unsigned AllocMemory(unsigned count,tsymatrix3f ** pt1,tsymatrix3f ** pt2=NULL){
  unsigned m=0;
  if(pt1){ unsigned s=sizeof(tsymatrix3f)*count; cudaMalloc(pt1,s); m+=s; }
  if(pt2){ unsigned s=sizeof(tsymatrix3f)*count; cudaMalloc(pt2,s); m+=s; }
  return(m);
}
//==============================================================================
unsigned AllocMemory(unsigned count,float ** pt1,float ** pt2=NULL,float ** pt3=NULL,float ** pt4=NULL){
  unsigned m=0;
  if(pt1){ unsigned s=sizeof(float)*count; cudaMalloc(pt1,s); m+=s; }
  if(pt2){ unsigned s=sizeof(float)*count; cudaMalloc(pt2,s); m+=s; }
  if(pt3){ unsigned s=sizeof(float)*count; cudaMalloc(pt3,s); m+=s; }
  if(pt4){ unsigned s=sizeof(float)*count; cudaMalloc(pt4,s); m+=s; }
  return(m);
}
//==============================================================================
unsigned AllocMemory(unsigned count,int2 ** pt1,int2 ** pt2=NULL,int2 ** pt3=NULL,int2 ** pt4=NULL){
  unsigned m=0;
  if(pt1){ unsigned s=sizeof(int2)*count; cudaMalloc(pt1,s); m+=s; }
  if(pt2){ unsigned s=sizeof(int2)*count; cudaMalloc(pt2,s); m+=s; }
  if(pt3){ unsigned s=sizeof(int2)*count; cudaMalloc(pt3,s); m+=s; }
  if(pt4){ unsigned s=sizeof(int2)*count; cudaMalloc(pt4,s); m+=s; }
  return(m);
}
//==============================================================================
unsigned AllocMemory(unsigned count,uint2 ** pt1,uint2 ** pt2=NULL,uint2 ** pt3=NULL,uint2 ** pt4=NULL){
  unsigned m=0;
  if(pt1){ unsigned s=sizeof(uint2)*count; cudaMalloc(pt1,s); m+=s; }
  if(pt2){ unsigned s=sizeof(uint2)*count; cudaMalloc(pt2,s); m+=s; }
  if(pt3){ unsigned s=sizeof(uint2)*count; cudaMalloc(pt3,s); m+=s; }
  if(pt4){ unsigned s=sizeof(uint2)*count; cudaMalloc(pt4,s); m+=s; }
  return(m);
}

//==============================================================================
///  Reserves the needed memory according to the number of particles.
//==============================================================================
void CsAllocMemoryBasic(StDeviceContext *dc,StDeviceCte *cte,bool svtimers){
  unsigned int m=0;
  CsOutResize(dc,2000);
  //-Data arrays (device).   
  m+=AllocMemory(dc->np,&dc->idp);
  m+=AllocMemory(dc->np,&dc->pos,&dc->vel);
  m+=AllocMemory(dc->np,&dc->rhop);
  if(dc->nmoving){
    m+=AllocMemory(dc->nmoving,&dc->ridpmv);
  }
  if(dc->nfloat){
    m+=AllocMemory(dc->nfloat,&dc->ridpft);
  }
  m+=AllocMemory(dc->npb,&dc->idsortb,&dc->cellb);
  m+=AllocMemory(dc->npf,&dc->idsortf,&dc->cellf);
  m+=AllocMemory(dc->np,&dc->ace,&dc->velxcor);
  m+=AllocMemory(dc->np,&dc->ar,&dc->viscdt,&dc->csound);
  if(dc->tstep==STEP_Verlet){
    m+=AllocMemory(dc->np,&dc->rhopm1,&dc->rhopnew);
    m+=AllocMemory(dc->np,&dc->velm1,&dc->velm1_2,&dc->velnew);
  }
  else if(dc->tstep==STEP_Symplectic){
    m+=AllocMemory(dc->np,&dc->pospre,&dc->pospre_2,&dc->velpre,&dc->velpre_2);
    m+=AllocMemory(dc->np,&dc->rhoppre);
  }
  if(dc->tvisco==VISCO_LaminarSPS){
    m+=AllocMemory(dc->npf,&dc->tau,&dc->csph);
  }
  if(dc->kgc){
    m+=AllocMemory(dc->npf,&dc->matkgc);
  }
  if(dc->nfloat){
    m+=AllocMemory(dc->nfloat,&dc->ftdist);
    m+=AllocMemory(1,&dc->face,&dc->fomegavel);
  }
  if(dc->shepardsteps){
    m+=AllocMemory(dc->npf,&dc->fdrhop);
  }
  m+=AllocMemory(dc->np,&dc->reduaux);
  dc->memgpu=m; //-Allocated GPU memory.
  CheckErrorCuda("Failed allocation of basic GPU memory.");
  TmgCreation(dc->timers,svtimers);
  CheckErrorCuda("Failed assignment of basic GPU memory.");
}

//==============================================================================
///  Reserves the needed memory according to the number of cells.
//==============================================================================
bool CsAllocMemoryCellMode(StDeviceContext *dc,StDeviceCte *cte){
  bool err=false;
  unsigned int m=0;
  //-Releases memory of the previous reservation.
  cudaFree(dc->cellbegb); dc->cellbegb=NULL;
  cudaFree(dc->cellbegf); dc->cellbegf=NULL;
  cudaFree(dc->cellnvb);  dc->cellnvb=NULL;
  cudaFree(dc->cellnvf);  dc->cellnvf=NULL;
  cudaFree(dc->pospres);  dc->pospres=NULL;
  cudaFree(dc->velrhop);  dc->velrhop=NULL;
  //-Allocates memory according to Cellmode.
  m+=AllocMemory(dc->nctotmax,&dc->cellbegb,&dc->cellbegf);
  dc->use_cellnv=(dc->cellmode==CELLMODE_Hneigs? 1: 0);
  if(dc->use_cellnv){
    m+=AllocMemory((dc->nctotmax-1)*dc->ncellnv,&dc->cellnvb,&dc->cellnvf);
  }
  if(dc->cellmode==CELLMODE_Hneigs||dc->cellmode==CELLMODE_2H||dc->cellmode==CELLMODE_H){
    m+=AllocMemory(dc->np,&dc->pospres,&dc->velrhop);
  } 
  dc->memgpu2=m;
  cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess){
    err=true;
    if(cuerr!=cudaErrorMemoryAllocation){
      char *cad=new char[2048]; 
      sprintf(cad,"%s (CUDA error: %s -> %s:%i).\n","Failed GPU memory allocation.",cudaGetErrorString(cuerr),__FILE__,__LINE__); 
      throw std::string(cad);
    }
  }
  return(err);
}

//------------------------------------------------------------------------------
/// Computes FtDist[] for a particles of a floating body.
//------------------------------------------------------------------------------
__global__ void KerCalcFtDist(unsigned pini,unsigned n,float3 center,float3 *ftdist,unsigned *ridpft,float3 *pos)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    p+=pini;
    float3 rpos=pos[ridpft[p]];
    ftdist[p]=make_float3(rpos.x-center.x,rpos.y-center.y,rpos.z-center.z); 
  }
}
//==============================================================================
/// Adjusts variables of particles floating body.
//==============================================================================
void CsInitFloating(StDeviceContext *dc){
  CsCalcRidpft(dc);
  for(unsigned cf=0;cf<dc->ftcount;cf++){
    StFloatingData *fobj=dc->ftobjs+cf;
    dim3 sgrid=CsGetGridSize(fobj->count,BLOCKSIZE);
    KerCalcFtDist<<<sgrid,BLOCKSIZE>>>(fobj->begin-dc->npb,fobj->count,Float3(fobj->center),dc->ftdist,dc->ridpft,dc->pos);
  }
}

//==============================================================================
/// Stores constants in GPU.
//==============================================================================
void CsUpCte(StDeviceCte *cte){
  cudaMemcpyToSymbol(CTE,cte,sizeof(StDeviceCte));
  CheckErrorCuda("Failed copying constants in GPU.");
}

//==============================================================================
/// Copies data to GPU to start simulation.
//==============================================================================
void CsUpData(StDeviceContext *dc,StDeviceCte *cte,unsigned *idp,float3 *pos,float3 *vel,float *rhop){
  cudaMemcpy(dc->pos,pos,sizeof(float3)*dc->npok,cudaMemcpyHostToDevice);
  cudaMemcpy(dc->vel,vel,sizeof(float3)*dc->npok,cudaMemcpyHostToDevice);
  cudaMemcpy(dc->rhop,rhop,sizeof(float)*dc->npok,cudaMemcpyHostToDevice);
  cudaMemcpy(dc->idp,idp,sizeof(unsigned)*dc->npok,cudaMemcpyHostToDevice);
  CheckErrorCuda("Failed copying data to GPU.");
  //-Initialises data for simulation.
  if(dc->tstep==STEP_Verlet){
    cudaMemcpy(dc->velm1,dc->vel,sizeof(float3)*dc->npok,cudaMemcpyDeviceToDevice);
    cudaMemcpy(dc->rhopm1,dc->rhop,sizeof(float)*dc->npok,cudaMemcpyDeviceToDevice);
    cudaMemset(dc->velnew,0,sizeof(float3)*dc->npok);
  }
  else if(dc->tstep==STEP_Symplectic){
    cudaMemcpy(dc->pospre,dc->pos,sizeof(float3)*dc->npok,cudaMemcpyDeviceToDevice);
    cudaMemcpy(dc->velpre,dc->vel,sizeof(float3)*dc->npok,cudaMemcpyDeviceToDevice);
    cudaMemcpy(dc->rhoppre,dc->rhop,sizeof(float)*dc->npok,cudaMemcpyDeviceToDevice);
  }
  if(dc->tvisco==VISCO_LaminarSPS)cudaMemset(dc->tau,0,sizeof(tsymatrix3f)*dc->npf);
  if(dc->nfloat)CsInitFloating(dc);
  CheckErrorCuda("Initialising variables for simulation.");
}


//##############################################################################
//# Generic functions.
//##############################################################################

//==============================================================================
/// Checks error and finishes execution.
//==============================================================================
void __CheckErrorCuda(const char *text,const char *file,const int line){
  cudaError_t err=cudaGetLastError();
  if(cudaSuccess!=err){
    char *cad=new char[2048]; 
    sprintf(cad,"%s (CUDA error: %s -> %s:%i).\n",text,cudaGetErrorString(err),file,line); 
    throw std::string(cad);
  }
}

//==============================================================================
/// Returns size of gridsize according to parameters.
//==============================================================================
dim3 CsGetGridSize(unsigned n,unsigned blocksize){
  dim3 sgrid;//=dim3(1,2,3);
  unsigned nb=unsigned(n+blocksize-1)/blocksize;//-Total number of blocks to be launched.
  sgrid.x=(nb<=GRIDSIZEDIM? nb: unsigned(sqrt(float(nb))));
  sgrid.y=(nb<=GRIDSIZEDIM? 1: unsigned((nb+sgrid.x-1)/sgrid.x));
  sgrid.z=1;
  return(sgrid);
}

//==============================================================================
/// Reduction of maximum of values of type float in shared memory for one warp.
//==============================================================================
template <unsigned blockSize> __device__ void WarpReduMaxF(volatile float* sax,unsigned tid) {
  if(blockSize>=64)sax[tid]=max(sax[tid],sax[tid+32]);
  if(blockSize>=32)sax[tid]=max(sax[tid],sax[tid+16]);
  if(blockSize>=16)sax[tid]=max(sax[tid],sax[tid+8]);
  if(blockSize>=8)sax[tid]=max(sax[tid],sax[tid+4]);
  if(blockSize>=4)sax[tid]=max(sax[tid],sax[tid+2]);
  if(blockSize>=2)sax[tid]=max(sax[tid],sax[tid+1]);
}

//==============================================================================
/// Acumulates the sum of n values of type float of dat[] and stores the result in res[0].
//==============================================================================
template <unsigned blockSize> __global__ void KerReduMaxF(unsigned n,float *dat,float *res){
  extern __shared__ float sax[];
  unsigned tid=threadIdx.x;
  unsigned c=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  sax[tid]=(c<n? dat[c]: -FLT_MAX);
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sax[tid]=max(sax[tid],sax[tid+256]);  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sax[tid]=max(sax[tid],sax[tid+128]);  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sax[tid]=max(sax[tid],sax[tid+64]);   __syncthreads(); }
  if(tid<32)WarpReduMaxF<blockSize>(sax,tid);
  if(tid==0)res[blockIdx.y*gridDim.x + blockIdx.x]=sax[0];
}

//==============================================================================
/// Returns the maximum of an array of type float using resu[] as auxiliar array.
//  size of resu[] = a (N/BLOCKSIZE+1)+(N/(BLOCKSIZE*BLOCKSIZE)+BLOCKSIZE)
//==============================================================================
float ReduMaxF(unsigned ndata,float* data,float* resu){
  unsigned n=ndata;
  unsigned smemSize=BLOCKSIZE*sizeof(float);
  dim3 sgrid=CsGetGridSize(n,BLOCKSIZE);
  unsigned n_blocks=sgrid.x*sgrid.y;
  float *dat=data;
  float *resu1=resu,*resu2=resu+n_blocks;
  float *res=resu1;
  while(n>1){
    KerReduMaxF<BLOCKSIZE><<<sgrid,BLOCKSIZE,smemSize>>>(n,dat,res);
    n=n_blocks;
    sgrid=CsGetGridSize(n,BLOCKSIZE);  
    n_blocks=sgrid.x*sgrid.y;
    if(n>1){
      dat=res; res=(dat==resu1? resu2: resu1); 
    }
  }
  float resf;
  if(ndata>1)cudaMemcpy(&resf,res,sizeof(float),cudaMemcpyDeviceToHost);
  else cudaMemcpy(&resf,data,sizeof(float),cudaMemcpyDeviceToHost);
  return(resf);
}

//==============================================================================
/// Reduction of maximum of values of type unsigned in shared memory for one warp.
//==============================================================================
template <unsigned blockSize> __device__ void WarpReduMaxU(volatile unsigned* sax,unsigned tid) {
  if(blockSize>=64)sax[tid]=max(sax[tid],sax[tid+32]);
  if(blockSize>=32)sax[tid]=max(sax[tid],sax[tid+16]);
  if(blockSize>=16)sax[tid]=max(sax[tid],sax[tid+8]);
  if(blockSize>=8)sax[tid]=max(sax[tid],sax[tid+4]);
  if(blockSize>=4)sax[tid]=max(sax[tid],sax[tid+2]);
  if(blockSize>=2)sax[tid]=max(sax[tid],sax[tid+1]);
}

//==============================================================================
/// Acumulates the sum of n values of type unsigned of dat[] and stores the result in res[0].
//==============================================================================
template <unsigned blockSize> __global__ void KerReduMaxU(unsigned n,unsigned *dat,unsigned *res){
  extern __shared__ unsigned sux[];
  unsigned tid=threadIdx.x;
  unsigned c=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  sux[tid]=(c<n? dat[c]: 0);
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sux[tid]=max(sux[tid],sux[tid+256]);  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sux[tid]=max(sux[tid],sux[tid+128]);  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sux[tid]=max(sux[tid],sux[tid+64]);   __syncthreads(); }
  if(tid<32)WarpReduMaxU<blockSize>(sux,tid);
  if(tid==0)res[blockIdx.y*gridDim.x + blockIdx.x]=sux[0];
}

//==============================================================================
/// Returns the maximum of an array of type unsigned using resu[] as auxiliar array.
//  size of resu[] = a (N/BLOCKSIZE+1)+(N/(BLOCKSIZE*BLOCKSIZE)+BLOCKSIZE)
//==============================================================================
unsigned ReduMaxU(unsigned ndata,unsigned* data,unsigned* resu){
  unsigned n=ndata;
  unsigned smemSize=BLOCKSIZE*sizeof(unsigned);
  dim3 sgrid=CsGetGridSize(n,BLOCKSIZE);
  unsigned n_blocks=sgrid.x*sgrid.y;
  unsigned *dat=data;
  unsigned *resu1=resu,*resu2=resu+n_blocks;
  unsigned *res=resu1;
  while(n>1){
    KerReduMaxU<BLOCKSIZE><<<sgrid,BLOCKSIZE,smemSize>>>(n,dat,res);
    n=n_blocks;
    sgrid=CsGetGridSize(n,BLOCKSIZE);  
    n_blocks=sgrid.x*sgrid.y;
    if(n>1){ dat=res; res=(dat==resu1? resu2: resu1); }
  }
  unsigned resf;
  if(ndata>1)cudaMemcpy(&resf,res,sizeof(unsigned),cudaMemcpyDeviceToHost);
  else cudaMemcpy(&resf,data,sizeof(unsigned),cudaMemcpyDeviceToHost);
  return(resf);
}







