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

/// \file CudaSphNL.cu \brief Implements all the functions to compute the neighbour list on GPU.

//==============================================================================
/// Modifies allocated memory for excluded particles Out.
//==============================================================================
void CsOutResize(StDeviceContext *dc,unsigned size){
  size=max(size,dc->nfluid);
  if(size>dc->outsize){
    unsigned *outidp=new unsigned[size];
    if(dc->outcount)memcpy(outidp,dc->outidp,sizeof(unsigned)*dc->outcount);
    delete[] dc->outidp; dc->outidp=outidp;
    float3 *outpos=new float3[size];
    if(dc->outcount)memcpy(outpos,dc->outpos,sizeof(float3)*dc->outcount);
    delete[] dc->outpos; dc->outpos=outpos;
    float3 *outvel=new float3[size];
    if(dc->outcount)memcpy(outvel,dc->outvel,sizeof(float3)*dc->outcount);
    delete[] dc->outvel; dc->outvel=outvel;
    float *outrhop=new float[size];
    if(dc->outcount)memcpy(outrhop,dc->outrhop,sizeof(float)*dc->outcount);
    delete[] dc->outrhop; dc->outrhop=outrhop;
    const unsigned sizepart=sizeof(unsigned)+sizeof(float3)+sizeof(float3)+sizeof(float);
    dc->memcpu-=dc->outsize*sizepart;
    dc->outsize=size;
    dc->memcpu+=dc->outsize*sizepart;
  }
}

//==============================================================================
/// Copies particles Out stored in buffer and empties the buffer.
//==============================================================================
void CsOutGetData(StDeviceContext *dc,unsigned *idp,float3 *pos,float3 *vel,float *rhop){
  memcpy(idp,dc->outidp,sizeof(unsigned)*dc->outcount);
  memcpy(pos,dc->outpos,sizeof(float3)*dc->outcount);
  memcpy(vel,dc->outvel,sizeof(float3)*dc->outcount);
  memcpy(rhop,dc->outrhop,sizeof(float)*dc->outcount);
  dc->outcount=0;
  CheckErrorCuda("Failed copying data to GPU.");
}

//------------------------------------------------------------------------------
/// CUDA Kernel that computes maximum Z of particles within the domain.
//------------------------------------------------------------------------------
template <unsigned int blockSize> __global__ void KerLimitZ(unsigned pini,unsigned n,float3 posmin,float3 difmax,const float3 *pos,float *limitz)
{
  extern __shared__ float spz[];
  unsigned tid=threadIdx.x;
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    float3 rpos=pos[p+pini];
    float dx=rpos.x-posmin.x,dy=rpos.y-posmin.y,dz=rpos.z-posmin.z;
    spz[tid]=((dx>=0 && dy>=0 && dz>=0 && dx<difmax.x && dy<difmax.y && dz<difmax.z)? rpos.z: -FLT_MAX);
  }
  else spz[tid]=-FLT_MAX;
  __syncthreads();
  if(blockSize>=512){ if(tid<256)spz[tid]=max(spz[tid],spz[tid+256]);  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)spz[tid]=max(spz[tid],spz[tid+128]);  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) spz[tid]=max(spz[tid],spz[tid+64]);   __syncthreads(); }
  if(tid<32)WarpReduMaxF<blockSize>(spz,tid);
  if(tid==0)limitz[blockIdx.y*gridDim.x+blockIdx.x]=spz[0];
}

//------------------------------------------------------------------------------
/// CUDA Kernel that loads idsort[] and cell[] to sort particles using radixsort.
//------------------------------------------------------------------------------
__global__ void KerPreSort(unsigned pini,unsigned n,unsigned nct,int ncx,int nsheet,float ovscell,float3 posmin,float3 difmax,const float3 *pos,unsigned *idsort,unsigned *cell)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    idsort[p]=p+pini;
    float3 rpos=pos[p+pini];
    float dx=rpos.x-posmin.x,dy=rpos.y-posmin.y,dz=rpos.z-posmin.z;
    cell[p]=((dx>=0 && dy>=0 && dz>=0 && dx<difmax.x && dy<difmax.y && dz<difmax.z)? int(dx*ovscell)+int(dy*ovscell)*ncx+int(dz*ovscell)*nsheet: nct);
  }
}

//------------------------------------------------------------------------------
/// CUDA Kernel that marks particle as excluded if it is out of the range of valid values for density.
//------------------------------------------------------------------------------
template<bool floating> __global__ void KerRhopOut(unsigned pini,unsigned n,unsigned nct,unsigned *cell,float *rhop,const float rhopmin,const float rhopmax,unsigned *idp)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    float rrhop=rhop[p+pini];
    if(rrhop<rhopmin||rrhop>rhopmax){
      if(floating? idp[p+pini]>=CTE.nbound: true)cell[p]=nct;//-Particles: Fluid.
    }
  }
}

//------------------------------------------------------------------------------
/// CUDA Kernel that reorders particle data according to idsort[].
//------------------------------------------------------------------------------
template<TpStep tstep> __global__ void KerSortData(unsigned pini,unsigned n,unsigned *idsort,unsigned *idp,float3 *pos,float3 *vel,float *rhop,float3 *velm1,float *rhopm1,float3 *pospre,float3 *velpre,float *rhoppre
    ,unsigned *idp2,float3 *pos2,float3 *vel2,float *rhop2,float3 *velm12,float *rhopm12,float3 *pospre2,float3 *velpre2,float *rhoppre2)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    unsigned oldpos=idsort[p];
    p+=pini;
    idp2[p]=idp[oldpos];
    pos2[p]=pos[oldpos];
    vel2[p]=vel[oldpos];
    rhop2[p]=rhop[oldpos];
    if(tstep==STEP_Verlet){
      velm12[p]=velm1[oldpos];
      rhopm12[p]=rhopm1[oldpos];
    }
    else if(tstep==STEP_Symplectic){
      pospre2[p]=pospre[oldpos];
      velpre2[p]=velpre[oldpos];
      rhoppre2[p]=rhoppre[oldpos];
    }
  }
}


//------------------------------------------------------------------------------
/// CUDA Kernel that reorders particle data according to idsort[] with Laminar+SPS viscosity.
//------------------------------------------------------------------------------
__global__ void KerSortDataLaminar(unsigned n,unsigned npb,unsigned *idsort,tsymatrix3f *tau,tsymatrix3f *tau2)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    unsigned oldpos=idsort[p];
    tau2[p]=tau[oldpos-npb];
  }
}

//------------------------------------------------------------------------------
/// CUDA Kernel that computes first particle and (last+1) of each cell.
//------------------------------------------------------------------------------
__global__ void KerCellsBeg(int pini,int n,int2 *cellbeg,unsigned *cell)
{
  extern __shared__ unsigned scell[];   // [blockDim.x+1}
  unsigned int p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  unsigned cel;
  if(p<n){
    cel=cell[p];
    scell[threadIdx.x+1]=cel;
    if(p&&!threadIdx.x)scell[0]=cell[p-1];
  }
  __syncthreads();
  if(p<n){
    unsigned p1=p+pini;
    if(!p||cel!=scell[threadIdx.x]){
      cellbeg[cel].x=p1;
      if(p)cellbeg[scell[threadIdx.x]].y=p1;
    }
    if(p==n-1)cellbeg[cel].y=p1+1;
  }
}

//------------------------------------------------------------------------------
/// CUDA Kernel that computes ranges of neighbouring particles for each cell.
//------------------------------------------------------------------------------
template<unsigned hdiv> __global__ void KerCellBegNv(unsigned nct,unsigned ncells,int ncx,int ncy,int ncz,uint2 *cellnv,int2 *cellbeg)
{
  unsigned cel=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(cel<nct){
    int nsheet=ncx*ncy;
    int cx=cel%ncx;
    int cz=int(cel/nsheet);
    int cy=int((cel%nsheet)/ncx);
    int cxini=cx-min(cx,hdiv);
    int cxfin=cx+min(ncx-cx-1,hdiv)+1;
    int yini=cy-min(cy,hdiv);
    int yfin=cy+min(ncy-cy-1,hdiv);
    int zini=cz-min(cz,hdiv);
    int zfin=cz+min(ncz-cz-1,hdiv);
    unsigned cnv=cel;
    uint2 rg;
    for(int z=zini;z<=zfin;z++){
      int zmod=nsheet*z;
      for(int y=yini;y<=yfin;y++){
        int ymod=zmod+ncx*y;
        rg.y=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=cellbeg[x+ymod];
          if(cbeg.y){
            if(!rg.y)rg.x=cbeg.x;
            rg.y=cbeg.y;
          }
        }
        if(rg.y){ cellnv[cnv]=rg; cnv+=ncells; }
      }
    }
    if(cnv<ncells*(hdiv==1? 9: 25)+cel){
      rg.y=0; cellnv[cnv]=rg;
    }
  }
}


//==============================================================================
/// Computes the Neighbour List of particles.
//==============================================================================
template<TpStep tstep,CsTypeDivide tdiv> void CsDivide(StDeviceContext *dc,StDeviceCte *cte){
  const char met[]="CsDivide";

  //-Obtains configuration of grid for executions of the kernels.
  //----------------------------------------------------------
  dim3 sgridb,sgridf=CsGetGridSize(dc->npf,BLOCKSIZE);
  if(tdiv==DV_All)sgridb=CsGetGridSize(dc->npb,BLOCKSIZE);

  //-Computes ncz and nct according to current position of the particles.
  //------------------------------------------------------------
  TmgStart(dc->timers,TMG_NlLimits);
  unsigned smemsize=BLOCKSIZE*sizeof(float);
  const float modz=(cte->h*2/dc->hdiv)*0.1f;  //-"modz" is used to get a margin due to precision error with real values.
  if(tdiv==DV_All){
    KerLimitZ<BLOCKSIZE> <<<sgridb,BLOCKSIZE,smemsize>>> (0,dc->npb,dc->posmin,dc->difmax,dc->pos,dc->csound);
    float maxz=ReduMaxF(sgridb.x*sgridb.y,dc->csound,dc->reduaux);
    dc->nczbound=unsigned(ceil((maxz-dc->posmin.z+modz)*dc->ovscell));
  }
  KerLimitZ<BLOCKSIZE> <<<sgridf,BLOCKSIZE,smemsize>>> (dc->npb,dc->npf,dc->posmin,dc->difmax,dc->pos,dc->csound);
  float maxz=ReduMaxF(sgridf.x*sgridf.y,dc->csound,dc->reduaux);
  unsigned nczfluid=unsigned(ceil((maxz-dc->posmin.z+modz)*dc->ovscell));
  unsigned ncznew=max(nczfluid,dc->nczbound);  
  if(dc->ncz!=ncznew){
    dc->ncz=ncznew;
    dc->nct=(dc->ncx*dc->ncy)*ncznew;
  }
  TmgStop(dc->timers,TMG_NlLimits);

  //-Prepares idsort[] and cell[] for reordering with radixsort (marks particles that will be excluded).
  //-------------------------------------------------------------------------------------------------------
  TmgStart(dc->timers,TMG_NlPreSort);
  if(tdiv==DV_All)KerPreSort <<<sgridb,BLOCKSIZE>>> (0,dc->npb,dc->nct,dc->ncx,(dc->ncx*dc->ncy),dc->ovscell,dc->posmin,dc->difmax,dc->pos,dc->idsortb,dc->cellb);
  KerPreSort <<<sgridf,BLOCKSIZE>>> (dc->npb,dc->npf,dc->nct,dc->ncx,(dc->ncx*dc->ncy),dc->ovscell,dc->posmin,dc->difmax,dc->pos,dc->idsortf,dc->cellf);

  //-Correction of density RhopOut (marks as Out particles that are out of the valid range of density).
  //-------------------------------------------------------------------------------------------------------
  if(dc->rhopout){
    if(dc->nfloat)KerRhopOut<true> <<<sgridf,BLOCKSIZE>>> (dc->npb,dc->npf,dc->nct,dc->cellf,dc->rhop,dc->rhopoutmin,dc->rhopoutmax,dc->idp);
    else KerRhopOut<false> <<<sgridf,BLOCKSIZE>>> (dc->npb,dc->npf,dc->nct,dc->cellf,dc->rhop,dc->rhopoutmin,dc->rhopoutmax,dc->idp);
  }
  TmgStop(dc->timers,TMG_NlPreSort);

  //-Orders idsort[] according to cell with radixsort.
  //-------------------------------------------------------------------------------------------------------
  TmgStart(dc->timers,TMG_NlRadixSort);
  if(tdiv==DV_All){  //cu40
    thrust::device_ptr<unsigned> dev_keysg(dc->cellb);
    thrust::device_ptr<unsigned> dev_valuesg(dc->idsortb);
    thrust::sort_by_key(dev_keysg,dev_keysg+dc->npb,dev_valuesg);
  }
  thrust::device_ptr<unsigned> dev_keysg(dc->cellf);
  thrust::device_ptr<unsigned> dev_valuesg(dc->idsortf);
  thrust::sort_by_key(dev_keysg,dev_keysg+dc->npf,dev_valuesg); 
  TmgStop(dc->timers,TMG_NlRadixSort);

  {
    //-Prepares auxiliar variables for data reordering.
    //--------------------------------------------------------
    if(sizeof(unsigned)!=sizeof(float))throw std::string("CsDivide> Unsigned and float types do not match in size.");
    unsigned *idp_2=(unsigned*)dc->viscdt;
    float3 *pos_2=dc->ace;
    float3 *vel_2=dc->velxcor;
    float *rhop_2=dc->ar;
    float *rhopm1_2=NULL;
    float *rhoppre_2=NULL;
    if(tstep==STEP_Verlet)rhopm1_2=dc->csound;
    if(tstep==STEP_Symplectic)rhoppre_2=dc->csound;

    //-Reorders data according to idsort[].
    //------------------------------
    TmgStart(dc->timers,TMG_NlSortData);
    if(tdiv==DV_All){
      KerSortData<tstep> <<<sgridb,BLOCKSIZE>>> 
        (0,dc->npb,dc->idsortb,dc->idp,dc->pos,dc->vel,dc->rhop,dc->velm1,dc->rhopm1,dc->pospre,dc->velpre,dc->rhoppre,idp_2,pos_2,vel_2,rhop_2,dc->velm1_2,rhopm1_2,dc->pospre_2,dc->velpre_2,rhoppre_2);
    }
    KerSortData<tstep> <<<sgridf,BLOCKSIZE>>> 
      (dc->npb,dc->npf,dc->idsortf,dc->idp,dc->pos,dc->vel,dc->rhop,dc->velm1,dc->rhopm1,dc->pospre,dc->velpre,dc->rhoppre,idp_2,pos_2,vel_2,rhop_2,dc->velm1_2,rhopm1_2,dc->pospre_2,dc->velpre_2,rhoppre_2);
    if(dc->tvisco==VISCO_LaminarSPS)KerSortDataLaminar <<<sgridf,BLOCKSIZE>>> (dc->npf,dc->npb,dc->idsortf,dc->tau,dc->csph);

    //-Places reordered data in the main variables.
    //----------------------------------------------------------
    if(tdiv!=DV_All){
      int size=sizeof(float)*dc->npb;
      int size3=sizeof(float3)*dc->npb;
      cudaMemcpy(idp_2,dc->idp,sizeof(unsigned)*dc->npb,cudaMemcpyDeviceToDevice);
      cudaMemcpy(pos_2,dc->pos,size3,cudaMemcpyDeviceToDevice);
      cudaMemcpy(vel_2,dc->vel,size3,cudaMemcpyDeviceToDevice);
      cudaMemcpy(rhop_2,dc->rhop,size,cudaMemcpyDeviceToDevice);
      if(tstep==STEP_Verlet)cudaMemcpy(rhopm1_2,dc->rhopm1,size,cudaMemcpyDeviceToDevice);
      if(tstep==STEP_Symplectic)cudaMemcpy(rhoppre_2,dc->rhoppre,size,cudaMemcpyDeviceToDevice);
    }
    unsigned* aux=dc->idp; dc->idp=idp_2; idp_2=aux;               // idp[] <=    idp_2[]
    float3* auxf3=dc->pos; dc->pos=pos_2; pos_2=auxf3;             // pos[] <=    pos_2[]
    auxf3=dc->vel; dc->vel=vel_2; vel_2=auxf3;                     // vel[] <=    vel_2[]
    float* auxf=dc->rhop; dc->rhop=rhop_2; rhop_2=auxf;            // rhop[] <=   rhop_2[]
    if(tstep==STEP_Verlet){
      auxf3=dc->velm1; dc->velm1=dc->velm1_2; dc->velm1_2=auxf3;   // velm1[] <=  velm1_2[]
      auxf=dc->rhopm1; dc->rhopm1=rhopm1_2; rhopm1_2=auxf;         // rhopm1[] <= rhopm1_2[]
    }
    else if(tstep==STEP_Symplectic){
      auxf3=dc->pospre; dc->pospre=dc->pospre_2; dc->pospre_2=auxf3;   // pospre[] <=  pospre_2[]
      auxf3=dc->velpre; dc->velpre=dc->velpre_2; dc->velpre_2=auxf3;   // velpre[] <=  velpre_2[]
      auxf=dc->rhoppre; dc->rhoppre=rhoppre_2; rhoppre_2=auxf;         // rhoppre[] <= rhoppre_2[]
    }
    if(dc->tvisco==VISCO_LaminarSPS){
      tsymatrix3f *auxm=dc->tau; dc->tau=dc->csph; dc->csph=auxm;      // tau[] <= csph[]
    }

    if(dc->lastdivall&&tdiv!=DV_All){
      if(tstep==STEP_Verlet)cudaMemcpy(dc->velm1,dc->velm1_2,sizeof(float3)*dc->npb,cudaMemcpyDeviceToDevice);
      if(tstep==STEP_Symplectic){
        cudaMemcpy(dc->pospre,dc->pospre_2,sizeof(float3)*dc->npb,cudaMemcpyDeviceToDevice);
        cudaMemcpy(dc->velpre,dc->velpre_2,sizeof(float3)*dc->npb,cudaMemcpyDeviceToDevice);
      }
    }
    TmgStop(dc->timers,TMG_NlSortData);

    //-Moves variables used as auxiliar to reorder data.
    //-----------------------------------------------------------------------
    dc->viscdt=(float*)idp_2;
    dc->ace=pos_2;
    dc->velxcor=vel_2;
    dc->ar=rhop_2;
    if(tstep==STEP_Verlet)dc->csound=rhopm1_2;
    if(tstep==STEP_Symplectic)dc->csound=rhoppre_2;
  }

  //-Computes first particle of each cell.
  //------------------------------------------
  TmgStart(dc->timers,TMG_NlCellBegin);
  if(tdiv==DV_All){
    cudaMemset(dc->cellbegb,0,sizeof(int2)*(dc->nctotmax));
    KerCellsBeg <<<sgridb,BLOCKSIZE,sizeof(int2)*(BLOCKSIZE+1)>>> (0,dc->npb,dc->cellbegb,dc->cellb);
  }
  cudaMemset(dc->cellbegf,0,sizeof(int2)*(dc->nct+1));
  KerCellsBeg <<<sgridf,BLOCKSIZE,sizeof(int2)*(BLOCKSIZE+1)>>> (dc->npb,dc->npf,dc->cellbegf,dc->cellf);
  TmgStop(dc->timers,TMG_NlCellBegin);

  //-Manages excluded particles.
  //-------------------------
  TmgStart(dc->timers,TMG_NlOutCheck);
  if(tdiv==DV_All){
    int2 lastcell2;
    cudaMemcpy(&lastcell2,dc->cellbegb+dc->nct,sizeof(int2),cudaMemcpyDeviceToHost);
    unsigned noutb=lastcell2.y-lastcell2.x;
    if(noutb){
      unsigned *vidp=new unsigned[noutb];
      float3 *vpos=new float3[noutb];
      char cad[256]; 
      printf("%s>%u boundary particles are out of the domain. (ndiv:%u)\n",met,noutb,dc->ndiv); 
      cudaMemcpy(vidp,dc->idp+lastcell2.x,noutb*sizeof(unsigned),cudaMemcpyDeviceToHost);
      cudaMemcpy(vpos,dc->pos+lastcell2.x,noutb*sizeof(float3),cudaMemcpyDeviceToHost);
      for(unsigned p=0;p<noutb;p++)printf("%s>Particle %u  id:%u  pos:(%G,%G,%G).\n",met,lastcell2.x+p,vidp[p],vpos[p].x,vpos[p].y,vpos[p].z); 
      fflush(stdout);
      sprintf(cad,"%s>%u boundary particles are out of the domain.",met,noutb); 
      throw std::string(cad);
    }
  }

  int2 lastcell2;
  cudaMemcpy(&lastcell2,dc->cellbegf+dc->nct,sizeof(int2),cudaMemcpyDeviceToHost);
  unsigned noutf=lastcell2.y-lastcell2.x;

  if(noutf){ //-There are new excluded particles.
    //-Copies data of excluded particles.
    if(dc->outcount+noutf>dc->outsize)CsOutResize(dc,max(dc->outcount+noutf+1000,dc->outsize*2));
    const unsigned pini=dc->npok-noutf;
    const unsigned pdes=dc->outcount;
    cudaMemcpy(dc->outidp+pdes,dc->idp+pini,sizeof(unsigned)*noutf,cudaMemcpyDeviceToHost);
    cudaMemcpy(dc->outpos+pdes,dc->pos+pini,sizeof(float3)*noutf,cudaMemcpyDeviceToHost);
    cudaMemcpy(dc->outvel+pdes,dc->vel+pini,sizeof(float3)*noutf,cudaMemcpyDeviceToHost);
    cudaMemcpy(dc->outrhop+pdes,dc->rhop+pini,sizeof(float)*noutf,cudaMemcpyDeviceToHost);
    dc->outcount+=noutf;
    dc->npok-=noutf; dc->npf-=noutf;
    if(dc->rhopout||dc->nfloat){
      const unsigned pfin=pdes+noutf;
      unsigned nerr=0;
      for(unsigned p=pdes;p<pfin;p++){
        if(dc->nfloat&&dc->outidp[p]<dc->nbound){
          nerr++;
          printf("FloatingOut[%u]> id:%u pos:(%f,%f,%f) rhop:%f\n",dc->ndiv,dc->outidp[p],dc->outpos[p].x,dc->outpos[p].y,dc->outpos[p].z,dc->outrhop[p]);
        }
        else if(dc->rhopout&&(dc->outrhop[p]<dc->rhopoutmin||dc->outrhop[p]>dc->rhopoutmax))dc->rhopoutcount++;
        if(nerr){
          char cad[256]; sprintf(cad,"%s>%u floating particles are out of the domain.",met,nerr); 
          throw std::string(cad);
        }
      }
    }
    if(dc->npf<=dc->nfloat){
      printf("\n**** Particles out: %d  (total: %d)\n",noutf,(dc->np-dc->npok));
      char cad[256]; sprintf(cad,"%s>There are no fluid particles in the domain.",met); 
      throw std::string(cad);
    } 
  }
  TmgStop(dc->timers,TMG_NlOutCheck);

  if(dc->use_cellnv){
    TmgStart(dc->timers,TMG_NlCellNv);
    const unsigned ncells=dc->nctotmax-1;
    if(tdiv==DV_All){
      dim3 sgrid=CsGetGridSize(ncells,BLOCKSIZE);
      cudaMemset(dc->cellnvb,0,sizeof(uint2)*ncells*dc->ncellnv);
      if(dc->hdiv==1)KerCellBegNv<1> <<<sgrid,BLOCKSIZE>>> (dc->nct,ncells,dc->ncx,dc->ncy,dc->ncz,dc->cellnvb,dc->cellbegb);
      else KerCellBegNv<2> <<<sgrid,BLOCKSIZE>>> (dc->nct,ncells,dc->ncx,dc->ncy,dc->ncz,dc->cellnvb,dc->cellbegb);
    }
    dim3 sgrid=CsGetGridSize(dc->nct,BLOCKSIZE);
    if(dc->hdiv==1)KerCellBegNv<1> <<<sgrid,BLOCKSIZE>>> (dc->nct,ncells,dc->ncx,dc->ncy,dc->ncz,dc->cellnvf,dc->cellbegf);
    else KerCellBegNv<2> <<<sgrid,BLOCKSIZE>>> (dc->nct,ncells,dc->ncx,dc->ncy,dc->ncz,dc->cellnvf,dc->cellbegf);
    TmgStop(dc->timers,TMG_NlCellNv);
  }

  dc->lastdivall=(tdiv==DV_All? 1: 0);
  if(tdiv==DV_All)dc->bounddivver=dc->bounddatver;
  CheckErrorCuda("Failed while neighbour list.");
}

//==============================================================================
/// Calls CsDivide() to compute the initial neighbour list.
//==============================================================================
void CsCallDivide(CsTypeDivide tdiv,StDeviceContext *dc,StDeviceCte *cte){
  switch(dc->tstep){
    case STEP_Verlet:
      switch(tdiv){
        case DV_All:    CsDivide<STEP_Verlet,DV_All>(dc,cte);    break;
        case DV_Fluid:  CsDivide<STEP_Verlet,DV_Fluid>(dc,cte);  break;
      }
    break;
    case STEP_Symplectic:
      switch(tdiv){
        case DV_All:    CsDivide<STEP_Symplectic,DV_All>(dc,cte);   break;
        case DV_Fluid:  CsDivide<STEP_Symplectic,DV_Fluid>(dc,cte); break;
      }
    break;
  }
}






