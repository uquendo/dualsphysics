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

/// \file CudaSphSU.cu \brief Implements all the functions to update system on GPU.

//------------------------------------------------------------------------------
/// CUDA Kernel that computes position of moving boundary particles according to idp[].
//------------------------------------------------------------------------------
__global__ void KerCalcRidp(unsigned pfirst,unsigned n,unsigned pini,unsigned pfin,unsigned *idp,unsigned *ridp)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    p+=pfirst; 
    unsigned id=idp[p];
    if(id>=pini&&id<pfin)ridp[id-pini]=p;
  }
}
//==============================================================================
/// Computes position of moving boundary particles according to idp[].
//==============================================================================
void CsCalcRidpmv(StDeviceContext *dc){
  dim3 sgrid=CsGetGridSize(dc->npb,BLOCKSIZE);
  KerCalcRidp <<<sgrid,BLOCKSIZE>>> (0,dc->npb,dc->nfixed,dc->nfixed+dc->nmoving,dc->idp,dc->ridpmv);
}
//==============================================================================
/// Computes position of particles of floating objects according to idp[].
//==============================================================================
void CsCalcRidpft(StDeviceContext *dc){
  dim3 sgrid=CsGetGridSize(dc->npf,BLOCKSIZE);
  KerCalcRidp <<<sgrid,BLOCKSIZE>>> (dc->npb,dc->npf,dc->npb,dc->npb+dc->nfloat,dc->idp,dc->ridpft);
}

//------------------------------------------------------------------------------
/// CUDA Kernel that computes values for a floating body.
//------------------------------------------------------------------------------
__global__ void KerCalcFtOmega(unsigned pini,unsigned n,float3 gravity,const float3 *ftdist,const unsigned *ridpft,const float3 *ace,float fmass,float3 *face,float3* fomegavel)
{
  unsigned tid=threadIdx.x;
  if(!tid){
    float rfacex=0,rfacey=0,rfacez=0;
    float rfomegavelx=0,rfomegavely=0,rfomegavelz=0;
    const unsigned pfin=pini+n;
    for(unsigned p=pini;p<pfin;p++){
      float3 race=ace[ridpft[p]];
      race.x-=gravity.x; race.y-=gravity.y; race.z-=gravity.z;
      rfacex+=race.x; rfacey+=race.y; rfacez+=race.z;
      float3 rdist=ftdist[p];
      rfomegavelx+=(race.z*rdist.y - race.y*rdist.z);
      rfomegavely+=(race.x*rdist.z - race.z*rdist.x);
      rfomegavelz+=(race.y*rdist.x - race.x*rdist.y);
    }
    rfacex=(rfacex+fmass*gravity.x)/fmass;
    rfacey=(rfacey+fmass*gravity.y)/fmass;
    rfacez=(rfacez+fmass*gravity.z)/fmass;
    face[0]=make_float3(rfacex,rfacey,rfacez);
    fomegavel[0]=make_float3(rfomegavelx,rfomegavely,rfomegavelz);
  }
}
//------------------------------------------------------------------------------
/// CUDA Kernel that updates particles of a floating body.
//------------------------------------------------------------------------------
template<bool predictor> __global__ void KerFtUpdate(unsigned pini,unsigned n,float dt,float3 center,float3 fvel,float3 fomega,unsigned *ridpft,float3 *pos,float3 *vel,float3 *ftdist)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    p+=pini;
    unsigned cp=ridpft[p];
    float3 rpos=pos[cp];
    float3 rvel=vel[cp];
    rpos.x+=dt*rvel.x;  rpos.y+=dt*rvel.y;  rpos.z+=dt*rvel.z;
    pos[cp]=rpos;
    float3 rdist=make_float3(rpos.x-center.x,rpos.y-center.y,rpos.z-center.z);  
    if(!predictor)ftdist[p]=rdist;
    rvel.x=fvel.x+(fomega.y*rdist.z-fomega.z*rdist.y);
    rvel.y=fvel.y+(fomega.z*rdist.x-fomega.x*rdist.z);
    rvel.z=fvel.z+(fomega.x*rdist.y-fomega.y*rdist.x);
    vel[cp]=rvel;
  }
}

//==============================================================================
/// Processes movement of particles of floating objects.
//==============================================================================
template<bool predictor> void CsRunFloating(StDeviceContext *dc,StDeviceCte *cte,float dt2){
  TmgStart(dc->timers,TMG_SuFloating);
  CsCalcRidpft(dc);
  for(unsigned cf=0;cf<dc->ftcount;cf++){
    StFloatingData *fobj=dc->ftobjs+cf;
    float3 center=Float3(fobj->center);
    float3 fvel=Float3(fobj->fvel);
    float3 fomega=Float3(fobj->fomega);
    float3 face,fomegavel;
    KerCalcFtOmega<<<1,32>>>(fobj->begin-dc->npb,fobj->count,dc->gravity,dc->ftdist,dc->ridpft,dc->ace,fobj->mass,dc->face,dc->fomegavel);
    cudaMemcpy(&face,dc->face,sizeof(float3),cudaMemcpyDeviceToHost);
    cudaMemcpy(&fomegavel,dc->fomegavel,sizeof(float3),cudaMemcpyDeviceToHost);
    fomegavel.x/=fobj->inertia.x;
    fomegavel.y/=fobj->inertia.y;
    fomegavel.z/=fobj->inertia.z;    
    center.x+=dt2*fvel.x;
    center.y+=dt2*fvel.y;
    center.z+=dt2*fvel.z;
    fvel.x+=dt2*face.x;
    fvel.y+=dt2*face.y;
    fvel.z+=dt2*face.z;
    fomega.x+=dt2*fomegavel.x;
    fomega.y+=dt2*fomegavel.y;
    fomega.z+=dt2*fomegavel.z;
    dim3 sgrid=CsGetGridSize(fobj->count,BLOCKSIZE);
    KerFtUpdate<predictor> <<<sgrid,BLOCKSIZE>>>(fobj->begin-dc->npb,fobj->count,dt2,center,fvel,fomega,dc->ridpft,dc->pos,dc->vel,dc->ftdist);
    if(!predictor){//-Stores data.
      fobj->center=TFloat3(center);
      fobj->fvel=TFloat3(fvel);
      fobj->fomega=TFloat3(fomega);
    }
  }
  TmgStop(dc->timers,TMG_SuFloating);
}

//==============================================================================
/// Copies particle data in GPU to the host.
//==============================================================================
void CsDownData(StDeviceContext *dc,unsigned pini,unsigned *idp,float3 *pos,float3 *vel,float *rhop){
  const unsigned n=dc->npok-pini;
  TmgStart(dc->timers,TMG_SuDownData);
  cudaMemcpy((pos+pini),(dc->pos+pini),sizeof(float3)*n,cudaMemcpyDeviceToHost);
  cudaMemcpy((vel+pini),(dc->vel+pini),sizeof(float3)*n,cudaMemcpyDeviceToHost);
  cudaMemcpy(rhop,dc->rhop,sizeof(float)*dc->npok,cudaMemcpyDeviceToHost);
  cudaMemcpy((idp+pini),(dc->idp+pini),sizeof(int)*n,cudaMemcpyDeviceToHost);
  TmgStop(dc->timers,TMG_SuDownData);
  CheckErrorCuda("Failed copying data to GPU.");
}

//==============================================================================
/// Recovers information about Rhop.
//==============================================================================
void CsDownDataRhop(StDeviceContext *dc,float *rhop){
  cudaMemcpy(rhop,dc->rhop,sizeof(float)*dc->npok,cudaMemcpyDeviceToHost);
}

//------------------------------------------------------------------------------
/// CUDA Kernel that computes Fa2=ace[].x*ace[].x+ace[].y*ace[].y+ace[].z*ace[].z
//------------------------------------------------------------------------------
__global__ void KerCalcFa2(unsigned n,float3 *ace,float *fa2)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    float3 race=ace[p];
    fa2[p]=race.x*race.x+race.y*race.y+race.z*race.z;
  }
}

//==============================================================================
/// Calculates average of values of DT.
//==============================================================================
void InfoDtMean(StInfoDtMean &sdt){
  if(sdt.count){
    float v=0; 
    for(int c=0;c<sdt.count;c++)v+=sdt.vdt[c];
    sdt.mean=((sdt.mean*sdt.meancount)+v)/(sdt.meancount+sdt.count);
    sdt.meancount+=sdt.count;
    sdt.count=0;
  }
}

//==============================================================================
/// Adds a new value of DT to StInfoDtMean.
//==============================================================================
void InfoDtAdd(float dt,StInfoDtMean &sdt){
  if(sdt.count>=INFODT_MAX)InfoDtMean(sdt);
  sdt.vdt[sdt.count]=dt;
  sdt.count++;
}

//==============================================================================
/// Adds new values of dt1 and dt2 to InfoDt.
//==============================================================================
void InfoDtsAdd(StInfoDt *idt,float dt1,float dt2){
  InfoDtAdd(dt1,idt->dt1);
  InfoDtAdd(dt2,idt->dt2);
}

//==============================================================================
/// Adds new value of dtstep to InfoDt.
//==============================================================================
void InfoDtStepAdd(StInfoDt *idt,float dtstep){
  if(idt->dt.count||idt->dt.meancount){
    if(idt->dtmin>dtstep)idt->dtmin=dtstep;
    if(idt->dtmax<dtstep)idt->dtmax=dtstep;
  }
  else{
    idt->dtmin=dtstep;
    idt->dtmax=dtstep;
  }
  InfoDtAdd(dtstep,idt->dt);
}

//==============================================================================
/// Computes a variable DT.
//==============================================================================
float CsDtVariable(StDeviceContext *dc,StDeviceCte *cte){
  TmgStart(dc->timers,TMG_SuCalcDt);
  const bool ALL=false;
  unsigned np=(ALL? dc->npok: dc->npf);
  unsigned pini=(ALL? 0: dc->npb);
  float maxcsound=ReduMaxF(np,dc->csound+pini,dc->reduaux);
  dim3 sgrid=CsGetGridSize(np,BLOCKSIZE);
  KerCalcFa2 <<<sgrid,BLOCKSIZE>>> (np,dc->ace+pini,dc->csound);
  float maxfa2=ReduMaxF(np,dc->csound,dc->reduaux);
  float maxviscdt=ReduMaxF(np,dc->viscdt+pini,dc->reduaux);
  float dt1=sqrt(cte->h)/sqrt(sqrt(maxfa2));
  float dt2=cte->h/(maxcsound+cte->h*maxviscdt);
  float dt=dc->cflnumber*min(dt1,dt2);
  if(dc->infodt)InfoDtsAdd(dc->infodt,dt1,dt2);
  if(dt<dc->dtmin){
    dt=dc->dtmin;
    dc->dtmodif++;
  }
  TmgStop(dc->timers,TMG_SuCalcDt);
  CheckErrorCuda("Failed computing the new value of dt.");
  return(dt);
}

//==============================================================================
/// Manages changes of velocity in Symplectic for boundary particles.
//==============================================================================
void CsSymplecticMotionPre(StDeviceContext *dc,bool all){
  cudaMemset(dc->vel,0,sizeof(float3)*dc->npb);
  if(all)cudaMemset(dc->velpre,0,sizeof(float3)*dc->npb);
}
//==============================================================================
void CsSymplecticMotionPost(StDeviceContext *dc){
  cudaMemcpy(dc->pospre,dc->pos,sizeof(float3)*dc->npb,cudaMemcpyDeviceToDevice);
  cudaMemcpy(dc->velpre,dc->vel,sizeof(float3)*dc->npb,cudaMemcpyDeviceToDevice);
}

//==============================================================================
/// Computes the new Neighbour List.
//==============================================================================
void CsCallComputeStepDivide(StDeviceContext *dc,StDeviceCte *cte){
  switch(dc->tstep){
    case STEP_Verlet:     
      if(dc->bounddivver>=dc->bounddatver)CsDivide<STEP_Verlet,DV_Fluid>(dc,cte);
      else CsDivide<STEP_Verlet,DV_All>(dc,cte);
    break;
    case STEP_Symplectic:
      if(dc->bounddivver>=dc->bounddatver)CsDivide<STEP_Symplectic,DV_Fluid>(dc,cte);
      else CsDivide<STEP_Symplectic,DV_All>(dc,cte);
    break;
  }
}

//------------------------------------------------------------------------------
/// CUDA Kernel that updates pos, velnew and rhopnew for CsComputeStep<STEP_Verlet> with floating.
//------------------------------------------------------------------------------
template<bool rhopbound,bool floating> __global__ void KerComputeStepVerlet(unsigned n,unsigned npb,float *ar,float *rhop,float *rhopnew,unsigned *idp,float3 *pos,float3 *vel,float3 *velxcor,float3 *ace,float3* velv, float3* velnew,float dt,float dt205,float dtdt)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    if(p<npb){//-Particles: Fixed & Moving
      if(rhopbound){
        float rrhop=rhop[p]+ar[p]*dtdt;
        rhopnew[p]=(rrhop<RHOPZERO? RHOPZERO: rrhop); //-To avoid unphysical behaviour of boundaries.
      }
      else rhopnew[p]=rhop[p];
    }
    else{ //-Particles: Floating & Fluid
      float rrhop=rhop[p]+ar[p]*dtdt;
      float3 rvel=vel[p];
      if(floating? idp[p]>=CTE.nbound: true){//-Particles: Fluid
        float3 rpos=pos[p],rvelxcor=velxcor[p],race=ace[p];
        rpos.x+=(rvel.x+rvelxcor.x*CTE.eps)*dt + race.x*dt205;
        rpos.y+=(rvel.y+rvelxcor.y*CTE.eps)*dt + race.y*dt205;
        rpos.z+=(rvel.z+rvelxcor.z*CTE.eps)*dt + race.z*dt205;
        pos[p]=rpos;
        float3 rvelv=velv[p];
        rvel.x=rvelv.x+race.x*dtdt;
        rvel.y=rvelv.y+race.y*dtdt;
        rvel.z=rvelv.z+race.z*dtdt;
        velnew[p]=rvel;
      }
      else{//-Particles: Floating
        rrhop=(rrhop<RHOPZERO? RHOPZERO: rrhop); //-To avoid unphysical behaviour of boundaries.
      }
      rhopnew[p]=rrhop;
      velnew[p]=rvel;
    }
  }
}
//------------------------------------------------------------------------------
/// CUDA Kernel that updates pos, vel and rhop in predictor of CsComputeStep<STEP_Symplectic>.
//------------------------------------------------------------------------------
template<bool rhopbound,bool floating> __global__ void KerComputeStepSymplectic(unsigned n,unsigned npb,unsigned *idp,float3 *pos,float3 *vel,float *rhop,float3 *pospre,float3 *velpre,float *rhoppre,float *ar,float3 *velxcor,float3 *ace,float dt2)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    if(p<npb){//-Particles: Fixed & Moving
      if(rhopbound){
        float rrhop=rhoppre[p]+ar[p]*dt2;
        rhop[p]=(rrhop<RHOPZERO? RHOPZERO: rrhop); //-To avoid unphysical behaviour of boundaries.
      }
      else rhop[p]=rhoppre[p];
    }
    else{ //-Particles: Floating & Fluid
      float rrhop=rhoppre[p]+ar[p]*dt2;
      float3 rpos=pospre[p],rvel=velpre[p];
      if(floating? idp[p]>=CTE.nbound: true){//-Particles: Fluid
        float3 rvelxcor=velxcor[p],race=ace[p];
        rpos.x+=(rvel.x+rvelxcor.x*CTE.eps)*dt2;
        rpos.y+=(rvel.y+rvelxcor.y*CTE.eps)*dt2;
        rpos.z+=(rvel.z+rvelxcor.z*CTE.eps)*dt2;
        rvel.x+=race.x*dt2;
        rvel.y+=race.y*dt2;
        rvel.z+=race.z*dt2;
      }
      else{//-Particles: Floating
        rrhop=(rrhop<RHOPZERO? RHOPZERO: rrhop); //-To avoid unphysical behaviour of boundaries.
      }
      rhop[p]=rrhop;
      pos[p]=rpos; 
      vel[p]=rvel;
    }
  }
}
//------------------------------------------------------------------------------
/// CUDA Kernel that updates pos, vel and rhop in corrector of CsComputeStep<STEP_Symplectic>.
//------------------------------------------------------------------------------
template<bool rhopbound,bool floating> __global__ void KerComputeStepSymplectic2(unsigned n,unsigned npb,unsigned *idp,float3 *pos,float3 *vel,float *rhop,float3 *pospre,float3 *velpre,float *rhoppre,float *ar,float3 *ace,float dt2,float dtpre)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    if(p<npb){//-Particles: Fixed & Moving
      if(rhopbound){
        float epsilon_rdot=(-ar[p]/rhop[p])*dtpre;
        float rrhop=rhoppre[p] * (2.f-epsilon_rdot)/(2.f+epsilon_rdot);
        rhop[p]=(rrhop<RHOPZERO? RHOPZERO: rrhop); //-To avoid unphysical behaviour of boundaries.
      }
      else rhop[p]=rhoppre[p];
    }
    else{ //-Particles: Floating & Fluid
      float epsilon_rdot=(-ar[p]/rhop[p])*dtpre;
      float rrhop=rhoppre[p] * (2.f-epsilon_rdot)/(2.f+epsilon_rdot);
      float3 rpos=pospre[p],rvel=velpre[p];
      if(floating? idp[p]>=CTE.nbound: true){//-Particles: Fluid
        float3 rvelp=rvel,race=ace[p];
        rvel.x+=race.x*dtpre;
        rvel.y+=race.y*dtpre;
        rvel.z+=race.z*dtpre;
        rpos.x+=(rvelp.x+rvel.x)*dt2;
        rpos.y+=(rvelp.y+rvel.y)*dt2;
        rpos.z+=(rvelp.z+rvel.z)*dt2;
      }
      else{//-Particles: Floating
        rrhop=(rrhop<RHOPZERO? RHOPZERO: rrhop); //-To avoid unphysical behaviour of boundaries.
      }
      rhop[p]=rrhop;
      pos[p]=rpos; 
      vel[p]=rvel;
    }
  }
}

//==============================================================================
/// Computes particle interaction and updates system using Verlet algorithm.
//==============================================================================
template<bool rhopbound,bool floating> float CsComputeStep_Ver(StDeviceContext *dc,StDeviceCte *cte){
  CsCallInteraction_Forces<true,floating>(dc,cte);  //-Particle interaction.
  float dt=CsDtVariable(dc,cte);                    //-Computes new dt.

  //-Updates position, velocity and density.
  TmgStart(dc->timers,TMG_SuComputeStep);
  float dtsq_05=0.5f*dt*dt;
  dim3 sgridf=CsGetGridSize(dc->npf,BLOCKSIZE);
  dim3 sgrid=CsGetGridSize(dc->npok,BLOCKSIZE);
  dc->verletstep++;
  if(dc->verletstep<dc->verletsteps){
    KerComputeStepVerlet<rhopbound,floating> <<<sgrid,BLOCKSIZE>>> (dc->npok,dc->npb,dc->ar,dc->rhopm1,dc->rhopnew,dc->idp,dc->pos,dc->vel,dc->velxcor,dc->ace,dc->velm1,dc->velnew,dt,dtsq_05,dt+dt);
  }
  else{
    KerComputeStepVerlet<rhopbound,floating> <<<sgrid,BLOCKSIZE>>> (dc->npok,dc->npb,dc->ar,dc->rhop,dc->rhopnew,dc->idp,dc->pos,dc->vel,dc->velxcor,dc->ace,dc->vel,dc->velnew,dt,dtsq_05,dt);
    dc->verletstep=0;
  }
  float3* auxf3=dc->velm1; dc->velm1=dc->vel; dc->vel=auxf3;   // VelM1[] <= Vel[]
  float* auxf=dc->rhopm1; dc->rhopm1=dc->rhop; dc->rhop=auxf;  // RhopM1[] <= Rhop[]
  auxf3=dc->vel; dc->vel=dc->velnew; dc->velnew=auxf3;         // Vel[] <= VelNew[]
  auxf=dc->rhop; dc->rhop=dc->rhopnew; dc->rhopnew=auxf;       // Rhop[] <= RhopNew[]
  if(dc->verletresetvel){
    cudaMemset(dc->velnew,0,sizeof(float3)*dc->npb);  //velnew[]=0
    dc->verletresetvel--;
  }
  TmgStop(dc->timers,TMG_SuComputeStep);
  
  if(floating)CsRunFloating<false>(dc,cte,dt);      //-Processes movement of floating objects.
  if(dc->infodt)InfoDtStepAdd(dc->infodt,dt);       //-Stores info of dt.
  CheckErrorCuda("Failed while step computing.");   
  return(dt);
}

//==============================================================================
/// Computes particle interaction and updates system using Symplectic algorithm.
//==============================================================================
template<bool rhopbound,bool floating> float CsComputeStep_Sym(StDeviceContext *dc,StDeviceCte *cte){
  const float dt=dc->dtpre;

  //-Predictor
  //-----------
  CsCallInteraction_Forces<true,floating>(dc,cte);    //-Particle interaction.
  const float dt_p=CsDtVariable(dc,cte);              //-Computes dt of predictor step. 
  //-Changes data to variables Pre to compute new data.
  TmgStart(dc->timers,TMG_SuComputeStep);
  float3* auxf3=dc->pospre; dc->pospre=dc->pos; dc->pos=auxf3;   // pospre[]  <= pos[]
  auxf3=dc->velpre; dc->velpre=dc->vel; dc->vel=auxf3;           // velpre[]  <= vel[]
  float* auxf=dc->rhoppre; dc->rhoppre=dc->rhop; dc->rhop=auxf;  // rhoppre[] <= rhop[]
  //-Updates position, velocity and density.
  const float dt2=dt*.5f;
  dim3 sgrid=CsGetGridSize(dc->npok,BLOCKSIZE);
  KerComputeStepSymplectic<rhopbound,floating> <<<sgrid,BLOCKSIZE>>> (dc->npok,dc->npb,dc->idp,dc->pos,dc->vel,dc->rhop,dc->pospre,dc->velpre,dc->rhoppre,dc->ar,dc->velxcor,dc->ace,dt2);
  TmgStop(dc->timers,TMG_SuComputeStep);
  if(floating)CsRunFloating<true>(dc,cte,dt2);        //-Processes movement of floating objects.

  //-Corrector
  //-----------
  CsDivide<STEP_Symplectic,DV_Fluid>(dc,cte);
  CsCallInteraction_Forces<false,floating>(dc,cte);   //-Particle interaction.

  //-Updates position, velocity and density.
  TmgStart(dc->timers,TMG_SuComputeStep);
  sgrid=CsGetGridSize(dc->npok,BLOCKSIZE);
  KerComputeStepSymplectic2<rhopbound,floating> <<<sgrid,BLOCKSIZE>>> (dc->npok,dc->npb,dc->idp,dc->pos,dc->vel,dc->rhop,dc->pospre,dc->velpre,dc->rhoppre,dc->ar,dc->ace,dt2,dc->dtpre);
  TmgStop(dc->timers,TMG_SuComputeStep);
  if(floating)CsRunFloating<false>(dc,cte,dt2);       //-Processes movement of floating objects.
  if(dc->nmoving)dc->bounddatver++;

  //-Computes dt for the next ComputeStep
  const float dt_c=CsDtVariable(dc,cte);
  dc->dtpre=min(dt_p,dt_c);
  if(dc->infodt)InfoDtStepAdd(dc->infodt,dt);
  CheckErrorCuda("Failed while step computing.");
  return(dt);
}

//==============================================================================
/// Calls CsComputeStep() depending on the time algorithm.
//==============================================================================
template<bool rhopbound,bool floating>float CsCallComputeStep2(StDeviceContext *dc,StDeviceCte *cte){
  if(dc->tstep==STEP_Verlet)return(CsComputeStep_Ver<rhopbound,floating>(dc,cte));
  else if(dc->tstep==STEP_Symplectic)return(CsComputeStep_Sym<rhopbound,floating>(dc,cte));
  return(0);
}
//==============================================================================
template<bool rhopbound>float CsCallComputeStep1(StDeviceContext *dc,StDeviceCte *cte){
  if(dc->nfloat)return(CsCallComputeStep2<rhopbound,true>(dc,cte));
  else return(CsCallComputeStep2<rhopbound,false>(dc,cte));
}
//==============================================================================
float CsCallComputeStep(StDeviceContext *dc,StDeviceCte *cte,bool rhopbound){
  if(rhopbound)return(CsCallComputeStep1<true>(dc,cte));
  else return(CsCallComputeStep1<false>(dc,cte));
}

//------------------------------------------------------------------------------
/// CUDA Kernel that modifies position and velocity of boundary particles due to simple movement.
//------------------------------------------------------------------------------
__global__ void KerMoveLinBound(unsigned n,unsigned pini,float3 mvsimple,float3 mvvel,unsigned *ridpmv,float3 *pos,float3 *vel)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    int pid=ridpmv[pini+p];
    float3 rpos=pos[pid];
    rpos.x+=mvsimple.x; rpos.y+=mvsimple.y; rpos.z+=mvsimple.z;
    pos[pid]=rpos;
    vel[pid]=mvvel;
  }
}
//==============================================================================
/// Modifies position and velocity of boundary particles due to simple movement.
//==============================================================================
void CsMoveLinBound(StDeviceContext *dc,unsigned pini,unsigned npar,float3 mvsimple,float3 mvvel){
  dim3 sgrid=CsGetGridSize(npar,BLOCKSIZE);
  if(dc->simulate2d){ mvsimple.y=0; mvvel.y=0; }
  KerMoveLinBound <<<sgrid,BLOCKSIZE>>> (npar,pini-dc->nfixed,mvsimple,mvvel,dc->ridpmv,dc->pos,dc->vel);
}

//------------------------------------------------------------------------------
/// CUDA Kernel that modifies position and velocity of boundary particles due to movement with matrix.
//------------------------------------------------------------------------------
template<bool simulate2d> __global__ void KerMoveMatBound(unsigned n,unsigned pini,tmatrix4f m,float dt,unsigned *ridpmv,float3 *pos,float3 *vel)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    unsigned pid=ridpmv[pini+p];
    float3 rpos=pos[pid];
    float3 rps;
    rps.x= m.a11*rpos.x + m.a12*rpos.y + m.a13*rpos.z + m.a14;
    rps.y= m.a21*rpos.x + m.a22*rpos.y + m.a23*rpos.z + m.a24;
    rps.z= m.a31*rpos.x + m.a32*rpos.y + m.a33*rpos.z + m.a34;
    pos[pid]=rps;
    rps.x=(rps.x-rpos.x)/dt;
    rps.y=(rps.y-rpos.y)/dt;
    rps.z=(rps.z-rpos.z)/dt;
    vel[pid]=rps;
    if(simulate2d){
      pos[pid].y=rpos.y;
      vel[pid].y=0;
    }
  }
}
//==============================================================================
/// Modifies position and velocity of boundary particles due to movement with matrix.
//==============================================================================
void CsMoveMatBound(StDeviceContext *dc,unsigned pini,unsigned npar,tmatrix4f mvmatrix,float dt){
  dim3 sgrid=CsGetGridSize(npar,BLOCKSIZE);
  if(dc->simulate2d)KerMoveMatBound<true> <<<sgrid,BLOCKSIZE>>> (npar,pini-dc->nfixed,mvmatrix,dt,dc->ridpmv,dc->pos,dc->vel);
  else KerMoveMatBound<false> <<<sgrid,BLOCKSIZE>>> (npar,pini-dc->nfixed,mvmatrix,dt,dc->ridpmv,dc->pos,dc->vel);
}






