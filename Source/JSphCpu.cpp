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

/// \file JSphCpu.cpp \brief Implements the class \ref JSphCpu

#include "JSphCpu.h"
#include "Functions.h"
#include "JSphMotion.h"
#include "JFormatFiles2.h"

//==============================================================================
/// Constructor.
//==============================================================================
JSphCpu::JSphCpu(){
  ClassName="JSphCpu";
  Cpu=true;
  Ace=NULL; Ar=NULL; VelXcor=NULL; 
  Csound=NULL; PrRhop=NULL; Tensil=NULL;
  VelM1=NULL; RhopM1=NULL; VelNew=NULL; RhopNew=NULL;  //-Verlet.
  PosPre=NULL; VelPre=NULL; RhopPre=NULL;              //-Symplectic.
  FdWab=NULL; FdRhop=NULL; FdVol=NULL;
  Tau=NULL; Csph=NULL;   //-Laminar+SPS.
  MatKgc=NULL;           //-KGC.
  Div=NULL;
  FtDist=NULL;           //-Floating objects.
  TmcCreation(Timers,false);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphCpu::~JSphCpu(){
  Reset();
  TmcDestruction(Timers);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphCpu::Reset(){
  JSph::Reset();
  delete[] Ace;     Ace=NULL;
  delete[] Ar;      Ar=NULL;
  delete[] VelXcor; VelXcor=NULL;
  delete[] Csound;  Csound=NULL;
  delete[] PrRhop;  PrRhop=NULL;
  delete[] Tensil;  Tensil=NULL;
  delete[] VelM1;   VelM1=NULL;   //-Verlet
  delete[] RhopM1;  RhopM1=NULL;  //-Verlet
  delete[] VelNew;  VelNew=NULL;  //-Verlet
  delete[] RhopNew; RhopNew=NULL; //-Verlet
  delete[] PosPre;  PosPre=NULL;  //-Symplectic
  delete[] VelPre;  VelPre=NULL;  //-Symplectic
  delete[] RhopPre; RhopPre=NULL; //-Symplectic
  delete[] FdWab;   FdWab=NULL;
  delete[] FdRhop;  FdRhop=NULL;
  delete[] FdVol;   FdVol=NULL;
  delete[] Tau;     Tau=NULL;     //-Laminar+SPS
  delete[] Csph;    Csph=NULL;    //-Laminar+SPS
  delete[] MatKgc;  MatKgc=NULL;  //-KGC
  delete Div; Div=NULL;
  delete[] FtDist; FtDist=NULL;
  BoundDivVer=0;  
  RunMode="";
}

//==============================================================================
/// Allocates memory of main data.
//==============================================================================
void JSphCpu::AllocMemory(){
  try{
    unsigned npite=Np*(OmpMode!=OMPM_Single? OmpThreads: 1);
    Ace=new tfloat3[npite];       MemCpuStatic+=sizeof(tfloat3)*npite;
    Ar=new float[npite];          MemCpuStatic+=sizeof(float)*npite;
    VelXcor=new tfloat3[npite];   MemCpuStatic+=sizeof(tfloat3)*npite;
    Csound=new float[Np];         MemCpuStatic+=sizeof(float)*Np;
    PrRhop=new float[Np];         MemCpuStatic+=sizeof(float)*Np;
    if(TKernel==KERNEL_Cubic){
      Tensil=new float[Np];    MemCpuStatic+=sizeof(float)*Np;
    } 
    if(TStep==STEP_Verlet){
      VelM1=new tfloat3[Np];   MemCpuStatic+=sizeof(tfloat3)*Np;
      RhopM1=new float[Np];    MemCpuStatic+=sizeof(float)*Np; 
      VelNew=new tfloat3[Np];  MemCpuStatic+=sizeof(tfloat3)*Np;
      RhopNew=new float[Np];   MemCpuStatic+=sizeof(float)*Np;
    } 
    else if(TStep==STEP_Symplectic){
      PosPre=new tfloat3[Np];  MemCpuStatic+=sizeof(tfloat3)*Np;
      VelPre=new tfloat3[Np];  MemCpuStatic+=sizeof(tfloat3)*Np;
      RhopPre=new float[Np];   MemCpuStatic+=sizeof(float)*Np;
    }
    if(ShepardSteps){   
      unsigned npf=Np-Npb;
      FdVol=new float[npf];         MemCpuStatic+=sizeof(float)*npf;
      unsigned nshepard=npf*(OmpMode!=OMPM_Single? OmpThreads: 1);
      FdWab=new float[nshepard];    MemCpuStatic+=sizeof(float)*nshepard;  
      FdRhop=new float[nshepard];   MemCpuStatic+=sizeof(float)*nshepard;
    }
    if(Nfloat){
      FtDist=new tfloat3[Nfloat];  MemCpuStatic+=sizeof(tfloat3)*Np;
    }
    if(TVisco==VISCO_LaminarSPS){ 
      const unsigned npf=Np-Npb;
      const unsigned npfomp=npf*(OmpMode!=OMPM_Single? OmpThreads: 1);
      Tau=new tsymatrix3f[npf];      MemCpuStatic+=sizeof(tsymatrix3f)*npf;  
      Csph=new tsymatrix3f[npfomp];  MemCpuStatic+=sizeof(tsymatrix3f)*npfomp;  
    }
    if(Kgc){     
      const unsigned npf=Np-Npb;
      const unsigned npfomp=npf*(OmpMode!=OMPM_Single? OmpThreads: 1);
      MatKgc=new tsymatrix3f[npfomp];   MemCpuStatic+=sizeof(tsymatrix3f)*npfomp;
    }
  }
  catch(const std::bad_alloc){
    RunException("AllocMemory","Could not allocate the requested memory.");
  }
  TmcCreation(Timers,SvTimers);
}

//==============================================================================
/// Computes the Neighbour List of boundary particles.
//==============================================================================
void JSphCpu::RunDivideBoundary(){
  TmcStart(Timers,TMC_NlDivideBound);
  Div->DivBoundary(Pos);
  TmcStop(Timers,TMC_NlDivideBound);
  TmcStart(Timers,TMC_NlSortDataBound);
  Div->SortBoundary(Idp);
  Div->SortBoundary(Pos);
  Div->SortBoundary(Vel);
  Div->SortBoundary(Rhop);
  if(TStep==STEP_Verlet){
    Div->SortBoundary(VelM1);
    Div->SortBoundary(RhopM1);
  }
  else if(TStep==STEP_Symplectic){
    Div->SortBoundary(PosPre);
    Div->SortBoundary(VelPre);
    Div->SortBoundary(RhopPre);
  }
  ClearRidpBound();
  BoundDivVer=BoundDatVer;
  TmcStop(Timers,TMC_NlSortDataBound);
}

//==============================================================================
/// Computes the Neighbour List of fluid particles.
//==============================================================================
void JSphCpu::RunDivideFluid(){
  TmcStart(Timers,TMC_NlDivideFluid);
  bool datmodif=Div->DivFluid(Idp,Pos,Vel,Rhop,TimeStep);
  TmcStop(Timers,TMC_NlDivideFluid);
  TmcStart(Timers,TMC_NlSortDataFluid);
  Div->SortFluid(Idp);
  Div->SortFluid(Pos);
  Div->SortFluid(Vel);
  Div->SortFluid(Rhop);
  if(TStep==STEP_Verlet){
    Div->SortFluid(VelM1);
    Div->SortFluid(RhopM1);
  }
  else if(TStep==STEP_Symplectic){
    Div->SortFluid(PosPre);
    Div->SortFluid(VelPre);
    Div->SortFluid(RhopPre);
  }
  if(TVisco==VISCO_LaminarSPS){ 
    Div->SortOnlyFluid(Tau);
  }
  if(datmodif){
    NpOk=Np-Div->GetNfluidOut();
    RhopOutCount=Div->GetRhopOutCount();
    Ncells=Div->GetNcells();
    Nsheet=Ncells.x*Ncells.y; Nct=Nsheet*Ncells.z;
  }
  ClearRidpFluid();
  TmcStop(Timers,TMC_NlSortDataFluid);
}

//==============================================================================
/// Computes the new Neighbour List.
//==============================================================================
void JSphCpu::ComputeStepDivide(){
  if(BoundDivVer<BoundDatVer)RunDivideBoundary();
  RunDivideFluid(); 
}

//==============================================================================
/// Cells interaction.
//==============================================================================
template<TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void JSphCpu::InteractionCells(){
#ifdef USE_OPENMP
  if(OmpMode==OMPM_Dynamic)InteractionCells_Dynamic<true,tinter,tker,tvis,floating>();
  if(OmpMode==OMPM_Static)InteractionCells_Static<true,tinter,tker,tvis,floating>();
#endif
  if(OmpMode==OMPM_Single)InteractionCells_Single<USE_SYMMETRY,tinter,tker,tvis,floating>();
}

//==============================================================================
/// Cells interaction with Single-Core.
//==============================================================================
template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void JSphCpu::InteractionCells_Single(){
  //-Shepard only for Fluid-Fluid
  for(byte kind=(tinter==INTER_Shepard? 2: 1);kind<=2;kind++){ //-1:Boundary, 2:Fluid
    byte kind2ini=(tinter==INTER_Shepard? 2: 3-kind); //(1->2, 2->1 y 2->2)
    if(tsym){//-Applies symmetry in force computation.
      for(int cz=0;cz<Ncells.z;cz++)for(int cy=0;cy<Ncells.y;cy++){
        int box_zy=cz*Nsheet+cy*Ncells.x;
        for(int cx=0;cx<Ncells.x;cx++){
          int box=box_zy+cx;
          if(Div->CellNoEmpty(box,kind)){
            int box2;
            int ibegin=Div->CellBegin(box,kind);
            int iend=Div->CellBegin(box+1,kind)-1;        
            int xneg=-min(cx,int(Hdiv));    
            int xpos=min(Ncells.x-cx-1,int(Hdiv));
            int yneg=-min(cy,int(Hdiv));
            int ypos=min(Ncells.y-cy-1,int(Hdiv));
            int zneg=-min(cz,int(Hdiv));
            int zpos=min(Ncells.z-cz-1,int(Hdiv));
            if(cx+1<Ncells.x)InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box+1,box+xpos,kind2ini); //East
            if(cy+1<Ncells.y){ //North
              for(int y=1;y<=ypos;y++){
                box2=box+Ncells.x*y;
                InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box2+xneg,box2+xpos,kind2ini);
              }
            }
            if(cz+1<Ncells.z){ //Up
              for(int z=1;z<=zpos;z++)for(int y=yneg;y<=ypos;y++){
                box2=box+Nsheet*z+Ncells.x*y;
                InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box2+xneg,box2+xpos,kind2ini);
              }
            }
          }
        }
      } 
      for(int box=0;box<Nct;box++)if(Div->CellNoEmpty(box,kind))InteractSelf<tsym,tinter,tker,tvis,floating>(box,kind,kind2ini);
    }
    else{//-No symmetry in force computation.
      int box=0;
      for(int cz=0;cz<Ncells.z;cz++)for(int cy=0;cy<Ncells.y;cy++)for(int cx=0;cx<Ncells.x;cx++){
        if(Div->CellNoEmpty(box,kind)){
          int box2;
          int ibegin=Div->CellBegin(box,kind);
          int iend=Div->CellBegin(box+1,kind)-1;
          int xneg=-min(cx,int(Hdiv));
          int xpos=min(Ncells.x-cx-1,int(Hdiv));
          int yneg=-min(cy,int(Hdiv));
          int ypos=min(Ncells.y-cy-1,int(Hdiv));
          int zneg=-min(cz,int(Hdiv));
          int zpos=min(Ncells.z-cz-1,int(Hdiv));
          for(int z=zneg;z<=zpos;z++)for(int y=yneg;y<=ypos;y++){
            box2=box+Nsheet*z+Ncells.x*y;
            InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box2+xneg,box2+xpos,kind2ini);
          }
        }
        box++;
      }
    }
  }
}

#ifdef USE_OPENMP
//==============================================================================
/// Cells interaction with OpenMp-Dynamic (with symmetry).
//==============================================================================
template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void JSphCpu::InteractionCells_Dynamic(){
  omp_set_num_threads(OmpThreads);
  //-Shepard only for Fluid-Fluid
  for(byte kind=(tinter==INTER_Shepard? 2: 1);kind<=2;kind++){ //-1:Boundary, 2:Fluid
    byte kind2ini=(tinter==INTER_Shepard? 2: 3-kind); //(1->2, 2->1 y 2->2)
    #pragma omp parallel for schedule (dynamic,10)
    for(int box=0;box<Nct;box++){
      const int cz=int(box/Nsheet);
      int bx=box-(cz*Nsheet);
      const int cy=int(bx/Ncells.x);
      const int cx=bx-(cy*Ncells.x);
      if(Div->CellNoEmpty(box,kind)){
        int box2;
        int ibegin=Div->CellBegin(box,kind);
        int iend=Div->CellBegin(box+1,kind)-1;        
        int xneg=-min(cx,int(Hdiv));
        int xpos=min(Ncells.x-cx-1,int(Hdiv));
        int yneg=-min(cy,int(Hdiv));
        int ypos=min(Ncells.y-cy-1,int(Hdiv));
        int zneg=-min(cz,int(Hdiv));
        int zpos=min(Ncells.z-cz-1,int(Hdiv));        
        InteractSelf<tsym,tinter,tker,tvis,floating>(box,kind,kind2ini);//Self
        if(cx+1<Ncells.x)InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box+1,box+xpos,kind2ini); //East
        if(cy+1<Ncells.y){ //North
          for(int y=1;y<=ypos;y++){
            box2=box+Ncells.x*y;
            InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box2+xneg,box2+xpos,kind2ini);
          }
        }
        if(cz+1<Ncells.z){ //Up
          for(int z=1;z<=zpos;z++)for(int y=yneg;y<=ypos;y++){
            box2=box+Nsheet*z+Ncells.x*y;
            InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box2+xneg,box2+xpos,kind2ini);
          }
        }
      }
    }
  }
}
//==============================================================================
/// Cells interaction with OpenMp-Static (with symmetry).
//==============================================================================
template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void JSphCpu::InteractionCells_Static(){
  omp_set_num_threads(OmpThreads);
  //-Shepard only for Fluid-Fluid
  for(byte kind=(tinter==INTER_Shepard? 2: 1);kind<=2;kind++){ //-1:Boundary, 2:Fluid
    byte kind2ini=(tinter==INTER_Shepard? 2: 3-kind); //(1->2, 2->1 y 2->2)
    #pragma omp parallel for schedule (static,10)
    for(int box=0;box<Nct;box++){
      const int cz=int(box/Nsheet);
      int bx=box-(cz*Nsheet);
      const int cy=int(bx/Ncells.x);
      const int cx=bx-(cy*Ncells.x);
      if(Div->CellNoEmpty(box,kind)){
        int box2;
        int ibegin=Div->CellBegin(box,kind);
        int iend=Div->CellBegin(box+1,kind)-1;        
        int xneg=-min(cx,int(Hdiv));
        int xpos=min(Ncells.x-cx-1,int(Hdiv));
        int yneg=-min(cy,int(Hdiv));
        int ypos=min(Ncells.y-cy-1,int(Hdiv));
        int zneg=-min(cz,int(Hdiv));
        int zpos=min(Ncells.z-cz-1,int(Hdiv));        
        InteractSelf<tsym,tinter,tker,tvis,floating>(box,kind,kind2ini);//Self
        if(cx+1<Ncells.x)InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box+1,box+xpos,kind2ini); //East
        if(cy+1<Ncells.y){ //North
          for(int y=1;y<=ypos;y++){
            box2=box+Ncells.x*y;
            InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box2+xneg,box2+xpos,kind2ini);
          }
        }
        if(cz+1<Ncells.z){ //Up
          for(int z=1;z<=zpos;z++)for(int y=yneg;y<=ypos;y++){
            box2=box+Nsheet*z+Ncells.x*y;
            InteractCelij<tsym,tinter,tker,tvis,floating>(ibegin,iend,box2+xneg,box2+xpos,kind2ini);
          }
        }
      }
    }
  }
}
#endif

//==============================================================================
/// Computes particle interactions of the same cell.
//==============================================================================
template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void JSphCpu::InteractSelf(int box,byte kind,byte kind2ini){
  float viscdt=0; 
#ifdef USE_OPENMP
  const int cth=omp_get_thread_num(),offset=(OmpMode!=OMPM_Single? Np*cth: 0);
  const int offsetf=(OmpMode!=OMPM_Single? (Np-Npb)*cth: 0);
#else
  const int cth=0,offset=0,offsetf=0;
#endif
  int ibegin=Div->CellBegin(box,kind);
  int iend=Div->CellBegin(box+1,kind)-1;
  for(byte kind2=max(kind,kind2ini);kind2<=2;kind2++){
    int jbegin,jend;
    if(kind!=kind2){
      jbegin=Div->CellBegin(box,kind2);
      jend=Div->CellBegin(box+1,kind2)-1;
    }
    else jend=iend;
    for(int i=ibegin;i<=iend;i++){
      TpParticle tpi=GetTpParticlePos(i);
      if(kind==kind2)jbegin=(tsym? i+1: ibegin);
      if(tinter==INTER_Shepard){//-Shepard only for Fluid-Fluid (no floating).
        if(tpi==PART_Fluid)for(int j=jbegin;j<=jend;j++){
          float drx=Pos[i].x-Pos[j].x, dry=Pos[i].y-Pos[j].y, drz=Pos[i].z-Pos[j].z;
          float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=Fourh2&&rr2>=1e-18f&&GetTpParticlePos(j)==PART_Fluid)ComputeForcesShepard<tsym,tker>(i,j,drx,dry,drz,rr2,offsetf);
        }
      }
      else{
        for(int j=jbegin;j<=jend;j++){
          float drx=Pos[i].x-Pos[j].x, dry=Pos[i].y-Pos[j].y, drz=Pos[i].z-Pos[j].z;
          float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=Fourh2&&rr2>=1e-18f){
            if(tinter==INTER_MatrixKgc)ComputeMatrixKgc<tsym,tker>(tpi,GetTpParticlePos(j),i,j,drx,dry,drz,rr2,offsetf);       
            else{
              float viscdtcf=ComputeForces<tsym,tinter,tker,tvis,floating>(tpi,GetTpParticlePos(j),i,j,drx,dry,drz,rr2,offset,offsetf);
              viscdt=max(viscdt,viscdtcf);
            }
          }
        }
      }
    }
  }
  ViscDtThread[cth*STRIDE_OMP]=max(ViscDtThread[cth*STRIDE_OMP],viscdt);
}

//==============================================================================
/// Computes particle interactions of a cell with other cells.
//==============================================================================
template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void JSphCpu::InteractCelij(int ibegin,int iend,int box1,int box2,byte kind2ini){
  float viscdt=0;
#ifdef USE_OPENMP
  const int cth=omp_get_thread_num(),offset=(OmpMode!=OMPM_Single? Np*cth: 0);
  const int offsetf=(OmpMode!=OMPM_Single? (Np-Npb)*cth: 0);
#else
  const int cth=0,offset=0,offsetf=0;
#endif
  for(byte kind=kind2ini;kind<=2;kind++){
    int jbegin=Div->CellBegin(box1,kind);
    int jend=Div->CellBegin(box2+1,kind)-1;
    if(jbegin<=jend){
      for(int i=ibegin;i<=iend;i++){
        TpParticle tpi=GetTpParticlePos(i);
        if(tinter==INTER_Shepard){//-Shepard only for Fluid-Fluid (no floating).
          if(tpi==PART_Fluid)for(int j=jbegin;j<=jend;j++){
            float drx=Pos[i].x-Pos[j].x, dry=Pos[i].y-Pos[j].y, drz=Pos[i].z-Pos[j].z;
            float rr2=drx*drx+dry*dry+drz*drz;
            if(rr2<=Fourh2&&rr2>=1e-18f&&GetTpParticlePos(j)==PART_Fluid)ComputeForcesShepard<tsym,tker>(i,j,drx,dry,drz,rr2,offsetf);
          }
        }
        else{
          for(int j=jbegin;j<=jend;j++){
            float drx=Pos[i].x-Pos[j].x, dry=Pos[i].y-Pos[j].y, drz=Pos[i].z-Pos[j].z;
            float rr2=drx*drx+dry*dry+drz*drz;
            if(rr2<=Fourh2&&rr2>=1e-18f){
              if(tinter==INTER_MatrixKgc)ComputeMatrixKgc<tsym,tker>(tpi,GetTpParticlePos(j),i,j,drx,dry,drz,rr2,offsetf);           
              else{
                float viscdtcf=ComputeForces<tsym,tinter,tker,tvis,floating>(tpi,GetTpParticlePos(j),i,j,drx,dry,drz,rr2,offset,offsetf);
                viscdt=max(viscdt,viscdtcf);
              }
            }
          }
        }
      }
    }
  }
  ViscDtThread[cth*STRIDE_OMP]=max(ViscDtThread[cth*STRIDE_OMP],viscdt);
}

//==============================================================================
/// Computes kernel and density summation for Shepard filter.
//==============================================================================
template<bool tsym,TpKernel tker> void JSphCpu::ComputeForcesShepard(int i,int j,float drx,float dry,float drz,float rr2,int offsetsh){
  float *fdwab=FdWab+offsetsh;
  float *fdrhop=FdRhop+offsetsh;
  const float rad=sqrt(rr2);
  const float qq=rad/H;
  float wab;   
  if(tker==KERNEL_Cubic){//-Cubic kernel.
    if(rad>H){
      float wqq1=2.0f-qq;
      wab=CubicCte.a24*(wqq1*wqq1*wqq1);
    }
    else wab=CubicCte.a2*(1.0f+(0.75f*qq-1.5f)*qq*qq);  //float wqq2=qq*qq;  wab=CubicCte.a2*(1.0f-1.5f*wqq2+0.75f*(wqq2*qq));   
  }
  if(tker==KERNEL_Wendland){//-Wendland kernel.
    float wqq=2*qq+1;
    float wqq1=1.f-0.5f*qq;
    float wqq2=wqq1*wqq1;
    wab=WendlandCte.awen*wqq*wqq2*wqq2;
  }
  const int i2=i-Npb,j2=j-Npb;
  fdwab[i2]+=wab*FdVol[j2];   
  fdrhop[i2]+=wab;//-MassFluid is multiplied in the end FdRhop[i]+=wab*MassFluid.
  if(tsym){
    fdwab[j2]+=wab*FdVol[i2];
    fdrhop[j2]+=wab;//-MassFluid is multiplied in the end  FdRhop[j]+=wab*MassFluid.
  }
}

//==============================================================================
/// Computes elements of the matrix used to apply the  kernel gradient correction (KGC).
//==============================================================================
template<bool tsym,TpKernel tker> void JSphCpu::ComputeMatrixKgc(TpParticle tpi,TpParticle tpj,int i,int j,float drx,float dry,float drz,float rr2,int offsetf){       
  tsymatrix3f *matkgc=MatKgc+offsetf;  
  //===== Kernel =====
  const float rad=sqrt(rr2);
  const float qq=rad/H;
  float fac;
  if(tker==KERNEL_Cubic){//-Cubic kernel.
    const float wqq1=2.0f-qq;
    fac=(rad>H? CubicCte.c2*(wqq1*wqq1): (CubicCte.c1+CubicCte.d1*qq)*qq) /rad;
  }
  if(tker==KERNEL_Wendland){//-Wendland kernel.
    const float wqq1=1.f-0.5f*qq;
    fac=WendlandCte.bwen*qq*wqq1*wqq1*wqq1/rad;
  }
  float frx=fac*drx,fry=fac*dry,frz=fac*drz;
  //===== KGC Matrix  =====
  if(tpi==PART_Fluid){//-Particle i is Fluid.
    const float massj=(tpj==PART_Fluid? MassFluid: MassBound);
    const float volj=-massj/Rhop[j];
    const unsigned i2=i-Npb;
    float fr=frx*volj; matkgc[i2].xx+=fr*drx; matkgc[i2].xy+=fr*dry; matkgc[i2].xz+=fr*drz; 
          fr=fry*volj; matkgc[i2].xy+=fr*drx; matkgc[i2].yy+=fr*dry; matkgc[i2].yz+=fr*drz; 
          fr=frz*volj; matkgc[i2].xz+=fr*drx; matkgc[i2].yz+=fr*dry; matkgc[i2].zz+=fr*drz; 
    // since KGC matrix must be symmetric to be inverted 
    // we already assume that matkgc[i2].yx=matkgc[i2].xy , matkgc[i2].zx=matkgc[i2].xz , matkgc[i2].zy=matkgc[i2].yz
    // so only 6 elements are needed instead of 3x3.
  }
  if(tsym&&tpj==PART_Fluid){//-If symmetry and particle j is Fluid.
    const float massi=(tpi==PART_Fluid? MassFluid: MassBound);
    const float voli=-massi/Rhop[i];
    const unsigned j2=j-Npb;
    float fr=frx*voli; matkgc[j2].xx+=fr*drx; matkgc[j2].xy+=fr*dry; matkgc[j2].xz+=fr*drz;
          fr=fry*voli; matkgc[j2].xy+=fr*drx; matkgc[j2].yy+=fr*dry; matkgc[j2].yz+=fr*drz;
          fr=frz*voli; matkgc[j2].xz+=fr*drx; matkgc[j2].yz+=fr*dry; matkgc[j2].zz+=fr*drz;
  }
}

//==============================================================================
/// Computes acceleration and density derivative during particle interactions.
/// [\ref INTER_Forces,\ref INTER_ForcesCorr,\ref INTER_ForcesKgc,\ref INTER_ForcesCorrKgc]
//==============================================================================
template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> float JSphCpu::ComputeForces(TpParticle tpi,TpParticle tpj,int i,int j,float drx,float dry,float drz,float rr2,int offset,int offsetf){
  float viscdt=0;
  float *ar=Ar+offset;  
  tfloat3 *ace=Ace+offset;     
  tfloat3 *velxcor=VelXcor+offset;
  const tsymatrix3f *matkgc=MatKgc; 
  const bool computei=(tpi&PART_BoundFt_Fluid)!=0;          //-Particle i is Fluid or BoundFt.   
  const bool computej=tsym&&((tpj&PART_BoundFt_Fluid)!=0);  //-Particle j is Fluid or BoundFt.   
  float prs=PrRhop[i]+PrRhop[j];       
  float wab;
  float frxi,fryi,frzi; //-new force values that will be corrected or not.
  float frxj,fryj,frzj;      
  {//===== Kernel =====
    float rad=sqrt(rr2);
    float qq=rad/H;
    float fac;
    if(tker==KERNEL_Cubic){//-Cubic kernel.
      if(rad>H){
        float wqq1=2.0f-qq;
        float wqq2=wqq1*wqq1;
        wab=CubicCte.a24*(wqq2*wqq1);
        fac=CubicCte.c2*wqq2/rad;
      }
      else{
        float wqq2=qq*qq;
        float wqq3=wqq2*qq;
        wab=CubicCte.a2*(1.0f-1.5f*wqq2+0.75f*wqq3);
        fac=(CubicCte.c1*qq+CubicCte.d1*wqq2)/rad;
      }
      //-Tensile correction.
      float fab=wab*CubicCte.od_wdeltap;
      fab*=fab; fab*=fab; //fab=fab^4
      prs+=fab*(Tensil[i]+Tensil[j]);
    }
    if(tker==KERNEL_Wendland){//-Wendland kernel.
      float wqq=2.f*qq+1.f;
      float wqq1=1.f-0.5f*qq;
      float wqq2=wqq1*wqq1;
      wab=WendlandCte.awen*wqq*wqq2*wqq2;
      fac=WendlandCte.bwen*qq*wqq2*wqq1/rad;
    }
    if(tinter==INTER_ForcesKgc||tinter==INTER_ForcesCorrKgc){
      float frx=fac*drx,fry=fac*dry,frz=fac*drz;
      if(tpi==PART_Fluid){
        const unsigned i2=i-Npb;
        frxi=matkgc[i2].xx*frx + matkgc[i2].xy*fry + matkgc[i2].xz*frz;      
        fryi=matkgc[i2].xy*frx + matkgc[i2].yy*fry + matkgc[i2].yz*frz;
        frzi=matkgc[i2].xz*frx + matkgc[i2].yz*fry + matkgc[i2].zz*frz;
      }
      else{ frxi=frx; fryi=fry; frzi=frz; }
      if(tpj==PART_Fluid){
        const unsigned j2=j-Npb;
        frxj=matkgc[j2].xx*frx + matkgc[j2].xy*fry + matkgc[j2].xz*frz;
        fryj=matkgc[j2].xy*frx + matkgc[j2].yy*fry + matkgc[j2].yz*frz;
        frzj=matkgc[j2].xz*frx + matkgc[j2].yz*fry + matkgc[j2].zz*frz;
      }
      else{ frxj=frx; fryj=fry; frzj=frz; }
    }
    else{
      frxi=frxj=fac*drx; fryi=fryj=fac*dry; frzi=frzj=fac*drz;
    }
  }    
  const float massi=(tpi==PART_Fluid? MassFluid: MassBound);
  const float massj=(tpj==PART_Fluid? MassFluid: MassBound);
  {//===== Aceleration ===== 
    const float p_vpmi=-prs*massj;
    if(tpi==PART_Fluid){                         //-Particle i is Fluid.
      ace[i].x+=p_vpmi*frxi; ace[i].y+=p_vpmi*fryi; ace[i].z+=p_vpmi*frzi;
    }
    else if(tpi==PART_BoundFt){                  //-Particle i is BoundFt. 
      ace[i].x+=p_vpmi*frxi*massi; ace[i].y+=p_vpmi*fryi*massi; ace[i].z+=p_vpmi*frzi*massi;          
    }
    if(tsym){
      const float p_vpmj=prs*massi;
      if(tpj==PART_Fluid){
        ace[j].x+=p_vpmj*frxj; ace[j].y+=p_vpmj*fryj; ace[j].z+=p_vpmj*frzj;    
      }
      else if(tpj==PART_BoundFt){
        ace[j].x+=p_vpmj*frxj*massj; ace[j].y+=p_vpmj*fryj*massj; ace[j].z+=p_vpmj*frzj*massj; 
      }
    }
  }
  const float dvx=Vel[i].x-Vel[j].x, dvy=Vel[i].y-Vel[j].y, dvz=Vel[i].z-Vel[j].z;
  const float robar=(Rhop[i]+Rhop[j])*0.5f;    
  {//===== Viscosity ===== 
    const float dot=drx*dvx + dry*dvy + drz*dvz;    
    //-Artificial viscosity. 
    if(tvis==VISCO_Artificial && dot<0){
      float cbar=(Csound[i]+Csound[j])*0.5f;
      float amubar=H*dot/(rr2+Eta2);
      float pi_visc=-Visco*cbar*amubar/robar;
      if(computei){
        const float v=-massj*pi_visc; ace[i].x+=v*frxi; ace[i].y+=v*fryi; ace[i].z+=v*frzi;  
      }
      if(computej){ 
        const float v=massi*pi_visc;  ace[j].x+=v*frxj; ace[j].y+=v*fryj; ace[j].z+=v*frzj;  
      }
    }
    //-Laminar+SPS viscosity.
    if(tvis==VISCO_LaminarSPS){
      tsymatrix3f *csph=Csph+offsetf;
      tsymatrix3f *tau=Tau;

      const float temp=2.0f*Visco/((rr2+Eta2)*robar);
      if(computei){
        const float vtemp=massj*temp*(drx*frxi+dry*fryi+drz*frzi);  
        ace[i].x+=vtemp*dvx; ace[i].y+=vtemp*dvy; ace[i].z+=vtemp*dvz;
      }
      if(computej){ 
        const float vtemp=-massi*temp*(drx*frxj+dry*fryj+drz*frzj); 
        ace[j].x+=vtemp*dvx; ace[j].y+=vtemp*dvy; ace[j].z+=vtemp*dvz; 
      }
      //-SPS turbulence model.  
      float tau_xx=0,tau_xy=0,tau_xz=0,tau_yy=0,tau_yz=0,tau_zz=0;
      const int i2=i-Npb;
      const int j2=j-Npb;
      if(computei){ 
        tau_xx+=tau[i2].xx; tau_xy+=tau[i2].xy; tau_xz+=tau[i2].xz;
        tau_yy+=tau[i2].yy; tau_yz+=tau[i2].yz; tau_zz+=tau[i2].zz;
      }
      if(computej){ 
        tau_xx+=tau[j2].xx; tau_xy+=tau[j2].xy; tau_xz+=tau[j2].xz;
        tau_yy+=tau[j2].yy; tau_yz+=tau[j2].yz; tau_zz+=tau[j2].zz;
      }
      if(computei){  
        ace[i].x+=massj*(tau_xx*frxi+tau_xy*fryi+tau_xz*frzi);
        ace[i].y+=massj*(tau_xy*frxi+tau_yy*fryi+tau_yz*frzi);
        ace[i].z+=massj*(tau_xz*frxi+tau_yz*fryi+tau_zz*frzi);
      }
      if(computej){ 
        ace[j].x+=-massi*(tau_xx*frxj+tau_xy*fryj+tau_xz*frzj);
        ace[j].y+=-massi*(tau_xy*frxj+tau_yy*fryj+tau_yz*frzj);
        ace[j].z+=-massi*(tau_xz*frxj+tau_yz*fryj+tau_zz*frzj);
      }
      //-Velocity gradients.       
      if(computei){
        const float volj=-massj/Rhop[j];
        float dv=dvx*volj; csph[i2].xx+=dv*frxi; csph[i2].xy+=dv*fryi; csph[i2].xz+=dv*frzi;
              dv=dvy*volj; csph[i2].xy+=dv*frxi; csph[i2].yy+=dv*fryi; csph[i2].yz+=dv*frzi;
              dv=dvz*volj; csph[i2].xz+=dv*frxi; csph[i2].yz+=dv*fryi; csph[i2].zz+=dv*frzi;
        // to compute tau terms we assume that csph.xy=csph.dudy+csph.dvdx, csph.xz=csph.dudz+csph.dwdx, csph.yz=csph.dvdz+csph.dwdy
        // so only 6 elements are needed instead of 3x3.
      }
      if(computej){ 
        const float voli=-massi/Rhop[i];
        float dv=dvx*voli; csph[j2].xx+=dv*frxj; csph[j2].xy+=dv*fryj; csph[j2].xz+=dv*frzj;
              dv=dvy*voli; csph[j2].xy+=dv*frxj; csph[j2].yy+=dv*fryj; csph[j2].yz+=dv*frzj;
              dv=dvz*voli; csph[j2].xz+=dv*frxj; csph[j2].yz+=dv*fryj; csph[j2].zz+=dv*frzj;
      }
    }
    viscdt=dot/(rr2+Eta2);
  }//===== Viscosity ===== 
  //-Density derivative.
  const float dot2i=(dvx*frxi+dvy*fryi+dvz*frzi);
  ar[i]+=massj*dot2i;
  if(tsym){
    const float dot2j=(dvx*frxj+dvy*fryj+dvz*frzj);
    ar[j]+=massi*dot2j;
  }
  //-XSPH correction.
  if(tinter==INTER_Forces||tinter==INTER_ForcesKgc){ //-XSPH correction only for Verlet or Symplectic-Predictor.
    const float wab_rhobar=wab/robar;
    if(computei){
      const float wab_rhobar_mass=-massj*wab_rhobar;
      velxcor[i].x+=wab_rhobar_mass * dvx;
      velxcor[i].y+=wab_rhobar_mass * dvy;
      velxcor[i].z+=wab_rhobar_mass * dvz;
    }
    if(computej){
      const float wab_rhobar_mass=massi*wab_rhobar;
      velxcor[j].x+=wab_rhobar_mass * dvx;
      velxcor[j].y+=wab_rhobar_mass * dvy;
      velxcor[j].z+=wab_rhobar_mass * dvz;
    }
  }
  return(viscdt);
}

#ifdef USE_OPENMP
//==============================================================================
/// Reduction of consecutive arrays of type float to the first one.
//==============================================================================
void JSphCpu::OmpMergeDataSum(int ini,int fin,float *data,int stride,int nthreads){
   #pragma omp parallel for schedule (static)
  for(int c=ini;c<fin;c++){
    for(int t=1,c2=c+stride;t<nthreads;t++,c2+=stride)data[c]+=data[c2];
  }
}
//==============================================================================
/// Reduction of consecutive arrays of type tfloat3 to the first one.
//==============================================================================
void JSphCpu::OmpMergeDataSum(int ini,int fin,tfloat3 *data,int stride,int nthreads){
  #pragma omp parallel for schedule (static)
  for(int c=ini;c<fin;c++){
    for(int t=1,c2=c+stride;t<nthreads;t++,c2+=stride)data[c]=data[c]+data[c2];
  }
}
//==============================================================================
/// Reduction of consecutive arrays of type tsymatrix3f to the first one.
//==============================================================================
void JSphCpu::OmpMergeDataSum(int ini,int fin,tsymatrix3f *data,int stride,int nthreads){
  #pragma omp parallel for schedule (static)
  for(int c=ini;c<fin;c++){
    tsymatrix3f dat=data[c];
    for(int t=1,c2=c+stride;t<nthreads;t++,c2+=stride){
      dat.xx+=data[c2].xx;
      dat.xy+=data[c2].xy;
      dat.xz+=data[c2].xz;
      dat.yy+=data[c2].yy;
      dat.yz+=data[c2].yz;
      dat.zz+=data[c2].zz;
    }
    data[c]=dat;
  }
}

#endif //-#ifdef USE_OPENMP

//==============================================================================
/// Computes a variable DT.
//==============================================================================
float JSphCpu::DtVariable(bool savedts){
  TmcStart(Timers,TMC_SuCalcDt);
  float fa_max=-FLT_MAX;
  float cs_max=-FLT_MAX;
  for(unsigned p=Npb;p<NpOk;p++){     
    float fa=(Ace[p].x*Ace[p].x+Ace[p].y*Ace[p].y+Ace[p].z*Ace[p].z);
    if(fa_max<=fa)fa_max=fa;
    if(cs_max<Csound[p])cs_max=Csound[p];
  } 
  float fa_sqrt=sqrt(sqrt(fa_max));
  float dt1=sqrt(H)/fa_sqrt;         //-dt1 depends on force per unit mass.
  float dt2=H/(cs_max+H*ViscDt);     //-dt2 combines the Courant and the viscous time-step controls.
  float dt=CFLnumber*min(dt1,dt2);   //-new value of time step.
  if(SvDt&&savedts)InfoDtsAdd(dt1,dt2);
  if(dt<DtMin){
    dt=DtMin;
    DtModif++;
  }
  TmcStop(Timers,TMC_SuCalcDt);
  return(dt);
}

//==============================================================================
/// Computes new values of density.
//==============================================================================
void JSphCpu::ComputeRhop(float* rhopnew,const float* rhopold,float armul,bool rhopbound){
  for(unsigned p=Npb;p<NpOk;p++)rhopnew[p]=rhopold[p]+Ar[p]*armul; 
  if(rhopbound)for(unsigned p=0;p<Npb;p++){
    rhopnew[p]=rhopold[p]+Ar[p]*armul;
    if(rhopnew[p]<RHOPZERO)rhopnew[p]=RHOPZERO; //-To avoid unphysical behaviour of boundaries.
  }
  else memcpy(rhopnew,rhopold,sizeof(float)*Npb);
  if(Nfloat){
    UpdateRidpFluid();
    for(unsigned c=Npb;c<Nbound;c++){ 
      const unsigned p=Ridp[c]; 
      if(rhopnew[p]<RHOPZERO)rhopnew[p]=RHOPZERO; //-To avoid unphysical behaviour of particles of floating objects.
    }
  }
}

//==============================================================================
/// Computes new value of density using epsilon (for corrector step of Symplectic).
//==============================================================================
void JSphCpu::ComputeRhopEpsilon(float* rhopnew,const float* rhopold,float armul,bool rhopbound){
  for(unsigned p=Npb;p<NpOk;p++){
    float epsilon_rdot=(-Ar[p]/rhopnew[p])*armul;
    rhopnew[p]=rhopold[p] * (2.f-epsilon_rdot)/(2.f+epsilon_rdot);
  }
  if(rhopbound)for(unsigned p=0;p<Npb;p++){
    float epsilon_rdot=(-Ar[p]/rhopnew[p])*armul;
    rhopnew[p]=rhopold[p] * (2.f-epsilon_rdot)/(2.f+epsilon_rdot);
    if(rhopnew[p]<RHOPZERO)rhopnew[p]=RHOPZERO;//-To avoid unphysical behaviour of boundaries.
  }
  else memcpy(rhopnew,rhopold,sizeof(float)*Npb);
  if(Nfloat){
    UpdateRidpFluid();
    for(unsigned c=Npb;c<Nbound;c++){ 
      const unsigned p=Ridp[c]; 
      if(rhopnew[p]<RHOPZERO)rhopnew[p]=RHOPZERO; //-To avoid unphysical behaviour of particles of floating objects.
    }
  }
}

//==============================================================================
/// Computes sub-particle stress tensor (Tau) for SPS turbulence model.
//==============================================================================
void JSphCpu::SPSCalcTau(){         
  const int npfok=(NpOk-Npb);
#ifdef USE_OPENMP
  #pragma omp parallel for schedule (static)
#endif
  for(int p=0;p<npfok;p++){
    const tsymatrix3f csph=Csph[p];
    const float pow1=csph.xx*csph.xx + csph.yy*csph.yy + csph.zz*csph.zz;
    const float prr=pow1+pow1 + csph.xy*csph.xy + csph.xz*csph.xz + csph.yz*csph.yz;
    const float visc_SPS=Smag*sqrt(prr);
    const float div_u=csph.xx+csph.yy+csph.zz;
    const float sps_k=(2.0f/3.0f)*visc_SPS*div_u;
    const float sps_Blin=Blin*prr;
    const float sumsps=-(sps_k+sps_Blin);
    const float twovisc_SPS=(visc_SPS+visc_SPS);
    const float one_rho2 = 1.0f/Rhop[p+Npb];   
    Tau[p].xx=one_rho2*(twovisc_SPS*csph.xx +sumsps);
    Tau[p].xy=one_rho2*(visc_SPS*csph.xy);
    Tau[p].xz=one_rho2*(visc_SPS*csph.xz);
    Tau[p].yy=one_rho2*(twovisc_SPS*csph.yy +sumsps);
    Tau[p].yz=one_rho2*(visc_SPS*csph.yz);
    Tau[p].zz=one_rho2*(twovisc_SPS*csph.zz +sumsps);
  }
}

//==============================================================================
/// Computes matrix of kernel gradient correction (KGC).
//==============================================================================
void JSphCpu::Interaction_MatrixKgc(){
  TmcStart(Timers,TMC_CfMatrixKgc);
  const unsigned npfomp=(Np-Npb)*(OmpMode!=OMPM_Single? OmpThreads: 1);
  memset(MatKgc,0,npfomp*sizeof(tsymatrix3f));    //MatKgc[]=0
  //-Interaction to compute elements of the matrix for KGC.
  CallInteractionCells(INTER_MatrixKgc);
  //-Reduction of MatKgc in only one matrix..
#ifdef USE_OPENMP
  if(OmpMode!=OMPM_Single)OmpMergeDataSum(0,NpOk-Npb,MatKgc,Np-Npb,OmpThreads);
#endif
  //-Matrix in case of 2D simulation.
  const int npfok=int(NpOk-Npb);
  if(Simulate2D){
    #ifdef USE_OPENMP
      #pragma omp parallel for schedule (static)
    #endif
    for(int p=0;p<npfok;p++)MatKgc[p].yy=1.0f;
  }   
  const tsymatrix3f imat={1,0,0 ,1,0 ,1}; //-Matrix when no correction will be applied.
  #ifdef USE_OPENMP
    #pragma omp parallel for schedule (static)
  #endif
  //-Inversion of a symmetric matrix.
  for(int p=0;p<npfok;p++){
    tsymatrix3f mat=MatKgc[p];
    // To ensure that matrix is symmetric
    // we assumed that matkgc.yx=matkgc.xy , matkgc.zx=matkgc.xz , matkgc.zy=matkgc.yz
    // and now we need the elements matkgc.xy=0.5*(matkgc.xy+matkgc.yx)...
    mat.xy*=0.5f;
    mat.xz*=0.5f;
    mat.yz*=0.5f;
    float det=mat.xx*mat.yy*mat.zz + 2.f*mat.xy*mat.yz*mat.xz - mat.xz*mat.yy*mat.xz - mat.xx*mat.yz*mat.yz - mat.xy*mat.xy*mat.zz; //-Determinant of the matrix.
    //-Matrix inversion.
    if(abs(det)>0.01 && abs(mat.xx)>MATRIXKGC_CONTROL && abs(mat.yy)>MATRIXKGC_CONTROL && abs(mat.zz)>MATRIXKGC_CONTROL){
      tsymatrix3f invmat;
      invmat.xx=(mat.yy*mat.zz-mat.yz*mat.yz)/det;
      invmat.xy=(mat.xz*mat.yz-mat.xy*mat.zz)/det;
      invmat.xz=(mat.xy*mat.yz-mat.yy*mat.xz)/det;
      invmat.yy=(mat.xx*mat.zz-mat.xz*mat.xz)/det;
      invmat.yz=(mat.xy*mat.xz-mat.xx*mat.yz)/det;
      invmat.zz=(mat.xx*mat.yy-mat.xy*mat.xy)/det;
      MatKgc[p]=invmat;
    }
    else MatKgc[p]=imat; //-No correction will be applied.
  }
  TmcStop(Timers,TMC_CfMatrixKgc);
}

//==============================================================================
/// Prepares variables for particle interaction "INTER_Forces" or "INTER_ForcesCorr".
//==============================================================================
void JSphCpu::PreInteraction_Forces(TpInter tinter){
  switch(TKernel){
    case KERNEL_Cubic:     PreInteraction_Forces_<KERNEL_Cubic>(tinter);     break;
    case KERNEL_Wendland:  PreInteraction_Forces_<KERNEL_Wendland>(tinter);  break;
  }
}
//==============================================================================
template<TpKernel tker> void JSphCpu::PreInteraction_Forces_(TpInter tinter){
  TmcStart(Timers,TMC_CfPreForces);
  if(OmpMode==OMPM_Dynamic||OmpMode==OMPM_Static){
    const int npomp=Np*OmpThreads;
    memset(Ar,0,sizeof(float)*npomp); //Ar[]=0
    if(tinter==INTER_Forces||tinter==INTER_ForcesKgc)memset(VelXcor,0,sizeof(tfloat3)*npomp); //VelXcor[]=(0,0,0)
    if(TVisco==VISCO_LaminarSPS)memset(Csph,0,sizeof(tsymatrix3f)*((Np-Npb)*OmpThreads));  //Csph[]={0,..,0}
    memset(Ace,0,sizeof(tfloat3)*npomp); //Ace[]=0 (only interval npb:np is needed)
  }
  else{
    memset(Ar,0,sizeof(float)*NpOk); //Ar[]=0   
    if(tinter==INTER_Forces||tinter==INTER_ForcesKgc)memset(VelXcor+Npb,0,sizeof(tfloat3)*(NpOk-Npb)); //VelXcor[]=(0,0,0)
    if(TVisco==VISCO_LaminarSPS)memset(Csph,0,sizeof(tsymatrix3f)*(NpOk-Npb)); //Csph[]={0,..,0}
  }
  for(unsigned p=Npb;p<NpOk;p++)Ace[p]=Gravity;  //Ace[]=Gravity
  ViscDt=0;
  for(int t=0;t<OmpThreads;t++)ViscDtThread[t*STRIDE_OMP]=0;

  //-Computes values of Csound[], PrRhop[p], Tensil[].
  const int npok=int(NpOk);
  #ifdef USE_OPENMP
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=0;p<npok;p++){
    float rhop=Rhop[p],rhop_r0=rhop*OVERRHOPZERO;
    Csound[p]=Cs0*(rhop_r0*rhop_r0*rhop_r0);
    float press=CteB*(pow(rhop_r0,Gamma)-1.0f);
    PrRhop[p]=press/(rhop*rhop);
    if(tker==KERNEL_Cubic)Tensil[p]=PrRhop[p]*(press>0? 0.01f: -0.2f);   
  }
  TmcStop(Timers,TMC_CfPreForces);
}

//==============================================================================
/// Particle interaction.
//==============================================================================
void JSphCpu::Interaction_Forces(bool forcescorr){
  TpInter tinter;
  if(Kgc){ //-Kernel Gradient Correction.
    Interaction_MatrixKgc();
    tinter=(forcescorr? INTER_ForcesCorrKgc: INTER_ForcesKgc); 
  }
  else tinter=(forcescorr? INTER_ForcesCorr: INTER_Forces); //-Kernel Gradient is not corrected. 
  PreInteraction_Forces(tinter); //-Prepares variables for force computation.
  TmcStart(Timers,TMC_CfForces); 
  CallInteractionCells(tinter); //-Force computation.
#ifdef USE_OPENMP
  if(OmpMode==OMPM_Dynamic||OmpMode==OMPM_Static){
    OmpMergeDataSum(0,NpOk,Ar,Np,OmpThreads); //-Reduction of Ar in only one array..
    OmpMergeDataSum(Npb,NpOk,Ace,Np,OmpThreads); //-Reduction of Ace in only one array..
    if(tinter==INTER_Forces||tinter==INTER_ForcesKgc)OmpMergeDataSum(Npb,NpOk,VelXcor,Np,OmpThreads); //-Reduction of VelXcor in only one array..
    if(TVisco==VISCO_LaminarSPS)OmpMergeDataSum(0,int(NpOk-Npb),Csph,int(Np-Npb),OmpThreads); //-Reduction of Csph in only one array..
  }
#endif
  if(TVisco==VISCO_LaminarSPS)SPSCalcTau(); //-Computes sub-particle stress tensor (Tau) for SPS turbulence model.
  for(int t=0;t<OmpThreads;t++)ViscDt=max(ViscDt,ViscDtThread[t*STRIDE_OMP]); //-Reduction of ViscDt in only one value..
  if(Simulate2D)for(unsigned p=Npb;p<Np;p++)Ace[p].y=0; //-Forces are cancelled in Y direction in case of 2D simulation.
  TmcStop(Timers,TMC_CfForces);
}

//==============================================================================
/// Computes new values of position and velocity using Verlet algorithm.
//==============================================================================
void JSphCpu::ComputeVerletVars(const tfloat3 *vel1,const tfloat3 *vel2,float dt,float dt2,tfloat3 *velnew){
  const float dtsq_05=0.5f*dt*dt;
  const int np=int(Np);
  #ifdef USE_OPENMP
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=int(Npb);p<np;p++){ //-Fluid particle or floating object.
    if(!Nfloat||Idp[p]>=Nbound){ //-PART_Fluid.
      Pos[p].x+=(vel1[p].x+VelXcor[p].x*Eps) * dt + Ace[p].x*dtsq_05;
      Pos[p].y+=(vel1[p].y+VelXcor[p].y*Eps) * dt + Ace[p].y*dtsq_05;
      Pos[p].z+=(vel1[p].z+VelXcor[p].z*Eps) * dt + Ace[p].z*dtsq_05;
      velnew[p].x=vel2[p].x+Ace[p].x*dt2;
      velnew[p].y=vel2[p].y+Ace[p].y*dt2;
      velnew[p].z=vel2[p].z+Ace[p].z*dt2;
    }
    else velnew[p]=vel1[p]; //-PART_BoundFt.
  }
}

//==============================================================================
/// Computes particle interaction and updates system using Verlet algorithm.
//==============================================================================
float JSphCpu::ComputeStep_Ver(bool rhopbound){
  Interaction_Forces(false);   //-Particle interaction.
  float dt=DtVariable(true);   //-Computes new dt.
  
  TmcStart(Timers,TMC_SuComputeStep);
  VerletStep++;
  if(VerletStep<VerletSteps){
    float twodt=dt+dt;
    ComputeVerletVars(Vel,VelM1,dt,twodt,VelNew); //-Updates position and velocity.
    ComputeRhop(RhopNew,RhopM1,twodt,rhopbound);  //-Updates density.
  }
  else{
    ComputeVerletVars(Vel,Vel,dt,dt,VelNew); //-Updates position and velocity.
    ComputeRhop(RhopNew,Rhop,dt,rhopbound);  //-Updates density.
    VerletStep=0;
  }
  swap(VelM1,Vel);    //VelM1[] <= Vel[]
  swap(RhopM1,Rhop);  //RhopM1[] <= Rhop[]
  swap(Vel,VelNew);   //Vel[] <= VelNew[]
  if(VerletResetVel){
    memset(VelNew,0,sizeof(tfloat3)*Npb);  //VelNew[]=0
    VerletResetVel--;
  }
  swap(Rhop,RhopNew); //Rhop[] <= RhopNew[]
  TmcStop(Timers,TMC_SuComputeStep);
  
  if(Nfloat)RunFloating(dt,false);  //-Processes movement of floating objects.
  if(SvDt)InfoDtStepAdd(dt);        //-Stores info of dt.
  return(dt);
}

//==============================================================================
/// Computes particle interaction and updates system using Symplectic algorithm.
//==============================================================================
float JSphCpu::ComputeStep_Sym(bool rhopbound){
  const float dt=DtPre;

  //-Predictor
  //-----------
  Interaction_Forces(false);           //-Particle interaction.
  const float dt_p=DtVariable(false);  //-Computes dt of predictor step.
  //-Changes data to variables Pre to compute new data.
  TmcStart(Timers,TMC_SuComputeStep);
  swap(PosPre,Pos);    //PosPre[] <= Pos[]
  swap(VelPre,Vel);    //VelPre[] <= Vel[]
  swap(RhopPre,Rhop);  //RhopPre[]<= Rhop[]
  const float dt2=dt*.5f;
  //-Updates position and velocity in predictor step.
  const int npok=int(NpOk);            //-Number of particles without the excluded ones.
  #ifdef USE_OPENMP
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=int(Npb);p<npok;p++){
    if(!Nfloat||Idp[p]>=Nbound){ //-PART_Fluid.
      Vel[p].x=VelPre[p].x + Ace[p].x * dt2; 
      Vel[p].y=VelPre[p].y + Ace[p].y * dt2; 
      Vel[p].z=VelPre[p].z + Ace[p].z * dt2; 
      Pos[p].x=PosPre[p].x + (VelPre[p].x+VelXcor[p].x*Eps) * dt2; 
      Pos[p].y=PosPre[p].y + (VelPre[p].y+VelXcor[p].y*Eps) * dt2; 
      Pos[p].z=PosPre[p].z + (VelPre[p].z+VelXcor[p].z*Eps) * dt2; 
    }
    else{ //-PART_BoundFt.
      Vel[p]=VelPre[p];
      Pos[p]=PosPre[p];
    }
  }
  //-Updates density in predictor step.
  ComputeRhop(Rhop,RhopPre,dt2,rhopbound); 
  TmcStop(Timers,TMC_SuComputeStep);
  if(Nfloat)RunFloating(dt2,true);     //-Processes movement of floating objects.

  //-Corrector
  //-----------
  RunDivideFluid(); 
  Interaction_Forces(true);            //-Particle interaction.
  const float dt_c=DtVariable(true);   //-Computes dt of corrector step.
  TmcStart(Timers,TMC_SuComputeStep);
  //-Updates position and velocity in corrector step.
  const int npok2=int(NpOk);
  #ifdef USE_OPENMP
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=Npb;p<npok2;p++){
    if(!Nfloat||Idp[p]>=Nbound){//-PART_Fluid.
      Vel[p].x=VelPre[p].x + Ace[p].x * DtPre; 
      Vel[p].y=VelPre[p].y + Ace[p].y * DtPre; 
      Vel[p].z=VelPre[p].z + Ace[p].z * DtPre; 
      Pos[p].x=PosPre[p].x + (VelPre[p].x+Vel[p].x) * dt2; 
      Pos[p].y=PosPre[p].y + (VelPre[p].y+Vel[p].y) * dt2; 
      Pos[p].z=PosPre[p].z + (VelPre[p].z+Vel[p].z) * dt2; 
    }
    else{//-PART_BoundFt.
      Vel[p]=VelPre[p];
      Pos[p]=PosPre[p];
    }
  }
  //-Updates density in corrector step.
  ComputeRhopEpsilon(Rhop,RhopPre,DtPre,rhopbound); 
  TmcStop(Timers,TMC_SuComputeStep);
  if(Nfloat)RunFloating(dt2,false);    //-Processes movement of floating objects.
  //-Computes dt for the next ComputeStep.
  const float stepdt=DtPre;
  DtPre=min(dt_p,dt_c);
  if(SvDt)InfoDtStepAdd(stepdt);
  return(stepdt);
}


//==============================================================================
/// Computes and updates density and velocity of probe particles.
//==============================================================================
void JSphCpu::ComputeProbes(){
        
}


//==============================================================================
/// Applies Shepard density filter.
//==============================================================================
void JSphCpu::RunShepard(){
  TmcStart(Timers,TMC_CfShepard);
  //-Prepares variables for Shepard interaction.
  const int npfok=int(NpOk-Npb),npf=int(Np-Npb);
  unsigned size=sizeof(float)*npf*(OmpMode!=OMPM_Single? OmpThreads: 1);
  memset(FdWab,0,size);  //FdWab[]=0
  memset(FdRhop,0,size); //FdRhop[]=0
  #ifdef USE_OPENMP
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=0;p<npfok;p++)FdVol[p]=MassFluid/Rhop[p+Npb]; //-Volume of the fluid particles. 
  //-Computes interactions for Shepard filter
  CallInteractionCells(INTER_Shepard);
  //-Reduction of FdWab and FdRhop in only one array. 
#ifdef USE_OPENMP
  if(OmpMode==OMPM_Dynamic||OmpMode==OMPM_Static){
    OmpMergeDataSum(0,npfok,FdWab,npf,OmpThreads);
    OmpMergeDataSum(0,npfok,FdRhop,npf,OmpThreads);
  }
#endif
  //-Reinitiliases density of fluid particles.
  if(!Nfloat){
    #pragma omp parallel for schedule (static)
    for(int p=0;p<npfok;p++){
      Rhop[p+Npb]=((FdRhop[p]+CteShepard)*MassFluid)/(FdWab[p]+(CteShepard*FdVol[p]));   
    }
  }
  else{
    #pragma omp parallel for schedule (static)
    for(int p=0;p<npfok;p++){ 
      if(Idp[p+Npb]>=Nbound){ //-Only fluid particles excluding floating objects.
        Rhop[p+Npb]=((FdRhop[p]+CteShepard)*MassFluid)/(FdWab[p]+(CteShepard*FdVol[p]));
      }
    }
  }
  TmcStop(Timers,TMC_CfShepard);
}

//==============================================================================
/// Processes movement of moving boundary particles.
//==============================================================================
void JSphCpu::RunMotion(float stepdt){
  TmcStart(Timers,TMC_SuMotion);
  if(Motion->ProcesTime(TimeStep+MotionTimeMod,stepdt)){
    unsigned nmove=Motion->GetMovCount();
    if(TStep==STEP_Symplectic){
      memset(Vel,0,sizeof(tfloat3)*Npb);    //Vel[]=0 when there is movement or when it finishes.
      if(!nmove)memset(VelPre,0,sizeof(tfloat3)*Npb); //VelPre[]=0
    }
    if(nmove){
      UpdateRidpBound();
      //-Movevement of boundary particles.
      for(unsigned c=0;c<nmove;c++){
        unsigned ref;
        tfloat3 mvsimple;
        tmatrix4f mvmatrix;
        if(Motion->GetMov(c,ref,mvsimple,mvmatrix)){//-Simple movement. 
          int plast=MotionObjBegin[ref+1];
          if(Simulate2D)for(int id=MotionObjBegin[ref];id<plast;id++){
            int pid=Ridp[id];
            tfloat3 *ps=(Pos+pid);
            ps->x+=mvsimple.x; ps->z+=mvsimple.z;
            Vel[pid]=TFloat3(mvsimple.x/stepdt,0,mvsimple.z/stepdt);
          }
          else for(int id=MotionObjBegin[ref];id<plast;id++){
            int pid=Ridp[id];
            tfloat3 *ps=(Pos+pid);
            ps->x+=mvsimple.x; ps->y+=mvsimple.y; ps->z+=mvsimple.z;
            Vel[pid]=TFloat3(mvsimple.x/stepdt,mvsimple.y/stepdt,mvsimple.z/stepdt);
          }
        }
        else{//-Movement with matrix.
          int plast=MotionObjBegin[ref+1];
          for(int id=MotionObjBegin[ref];id<plast;id++){
            int pid=Ridp[id];
            tfloat3 *ps=(Pos+pid);
            tfloat3 ps2=MatrixMulPoint(mvmatrix,*ps);
            Vel[pid]=TFloat3((ps2.x-ps->x)/stepdt,(ps2.y-ps->y)/stepdt,(ps2.z-ps->z)/stepdt);
            if(Simulate2D){ ps2.y=ps->y; Vel[pid].y=0; }
            *ps=ps2;
          }
        }
      }
      BoundDatVer++;//-Only neighbour list of boundaris must be updated. 
      if(TStep==STEP_Verlet)VerletResetVel=2;
      if(TStep==STEP_Symplectic){
        memcpy(PosPre,Pos,sizeof(tfloat3)*Npb);  //-Swap(PosPre,Pos)
        memcpy(VelPre,Vel,sizeof(tfloat3)*Npb);  //-Swap(VelPre,Vel)
      }
    }
  }
  TmcStop(Timers,TMC_SuMotion);
}

//==============================================================================
/// Adjusts variables of particles of floating bodies.
//==============================================================================
void JSphCpu::InitFloating(){
  UpdateRidp();
  for(unsigned c=0;c<FtCount;c++){
    const StFloatingData* fobj=FtObjs+c;
    unsigned idlast=fobj->begin+fobj->count-1;
    for(unsigned id=fobj->begin;id<=idlast;id++){
      int p=Ridp[id];
      FtDist[id-Npb]=TFloat3(Pos[p].x-fobj->center.x,Pos[p].y-fobj->center.y,Pos[p].z-fobj->center.z);
    }
  }
}

//==============================================================================
/// Processes movement of particles of floating objects.
//==============================================================================
void JSphCpu::RunFloating(float dt2,bool predictor){
  TmcStart(Timers,TMC_SuFloating);
  UpdateRidpFluid();
  for(unsigned cf=0;cf<FtCount;cf++){
    StFloatingData *fobj=FtObjs+cf;
    tfloat3 center=fobj->center;
    tfloat3 fvel=fobj->fvel;
    tfloat3 fomega=fobj->fomega;
    tfloat3 fomegavel=TFloat3(0);
    tfloat3 face=TFloat3(0);
    const unsigned pend=fobj->begin+fobj->count-1;
    for(unsigned id=fobj->begin;id<=pend;id++){
      int p=Ridp[id];      
      float acex=Ace[p].x-Gravity.x,acey=Ace[p].y-Gravity.y,acez=Ace[p].z-Gravity.z;   //-Ace was previously initialised with value of gravity.
      face.x+=acex; face.y+=acey; face.z+=acez;
      tfloat3 *dist=FtDist+(id-Npb);
      fomegavel.x+= acez*dist->y - acey*dist->z;
      fomegavel.y+= acex*dist->z - acez*dist->x;
      fomegavel.z+= acey*dist->x - acex*dist->y;
    }
    face.x=(face.x+fobj->mass*Gravity.x)/fobj->mass;
    face.y=(face.y+fobj->mass*Gravity.y)/fobj->mass;
    face.z=(face.z+fobj->mass*Gravity.z)/fobj->mass;
    fomegavel.x/=fobj->inertia.x;
    fomegavel.y/=fobj->inertia.y;
    fomegavel.z/=fobj->inertia.z;
    if(Simulate2D){ face.y=0; fomegavel.x=0; fomegavel.z=0; fvel.y=0; }
    center.x+=dt2*fvel.x;
    center.y+=dt2*fvel.y;
    center.z+=dt2*fvel.z;
    fvel.x+=dt2*face.x;
    fvel.y+=dt2*face.y;
    fvel.z+=dt2*face.z;
    fomega.x+=dt2*fomegavel.x;
    fomega.y+=dt2*fomegavel.y;
    fomega.z+=dt2*fomegavel.z;
    if(Simulate2D)for(unsigned id=fobj->begin;id<=pend;id++){
      int p=Ridp[id];
      tfloat3 *pos=Pos+p,*vel=Vel+p;
      pos->x+=dt2*vel->x;  pos->z+=dt2*vel->z;
      tfloat3 distaux;
      tfloat3 *dist=(predictor? &distaux: FtDist+(id-Npb)); 
      *dist=TFloat3(pos->x-center.x,0,pos->z-center.z); 
      vel->x=fvel.x+(fomega.y*dist->z-fomega.z*dist->y);
      vel->y=0;
      vel->z=fvel.z+(fomega.x*dist->y-fomega.y*dist->x);
    }
    else for(unsigned id=fobj->begin;id<=pend;id++){
      int p=Ridp[id];
      tfloat3 *pos=Pos+p,*vel=Vel+p;
      pos->x+=dt2*vel->x;  pos->y+=dt2*vel->y;  pos->z+=dt2*vel->z;
      tfloat3 distaux;
      tfloat3 *dist=(predictor? &distaux: FtDist+(id-Npb)); 
      *dist=TFloat3(pos->x-center.x,pos->y-center.y,pos->z-center.z);   
      vel->x=fvel.x+(fomega.y*dist->z-fomega.z*dist->y);
      vel->y=fvel.y+(fomega.z*dist->x-fomega.x*dist->z);
      vel->z=fvel.z+(fomega.x*dist->y-fomega.y*dist->x);
    }
    if(!predictor){//-Stores computed values.
      fobj->center=center;
      fobj->fvel=fvel;
      fobj->fomega=fomega;
    }
  }
  TmcStop(Timers,TMC_SuFloating);
}

//==============================================================================
/// Configuration for Neighbour List.
//==============================================================================
void JSphCpu::DivideConfig(){
  Log->Print("");
  Div=new JDivideCpu(Log,DirOut);
  Div->Config(Dosh/float(Hdiv),IncZ,Np,Nfixed,Nmoving,Nfloat,Nfluid,Pos,Dp);
  Div->ConfigRhopOut(RhopOut,RhopOutMin,RhopOutMax);
  Div->ConfigMemory();
  Ncells=Div->GetNcells();
  Nsheet=Ncells.x*Ncells.y; Nct=Nsheet*Ncells.z;
  JFormatFiles2::SaveVtkCells(DirOut+"MapCells.vtk",Div->GetPosMin(),ToTUint3(Ncells),Dosh/float(Hdiv)); //-MapCells.vtk shows the cells division of the domain.
}

//==============================================================================
/// Initialisation of arrays and variables for the execution.
//==============================================================================
void JSphCpu::InitVars(){
  memset(Ace,0,sizeof(tfloat3)*Np);     //Ace[]=(0,0,0)
  memset(VelXcor,0,sizeof(tfloat3)*Np); //VelXcor[]=(0,0,0)
  if(TStep==STEP_Verlet){
    memcpy(VelM1,Vel,sizeof(tfloat3)*Np);
    memcpy(RhopM1,Rhop,sizeof(float)*Np);
    memset(VelNew,0,sizeof(tfloat3)*Np);   
    VerletStep=0;
    VerletResetVel=0;
  }
  else if(TStep==STEP_Symplectic){
    DtPre=DtIni;
    memcpy(PosPre,Pos,sizeof(tfloat3)*Np); //swap(PosPre,Pos)
    memcpy(VelPre,Vel,sizeof(tfloat3)*Np); //swap(VelPre,Vel)
    memcpy(RhopPre,Rhop,sizeof(float)*Np); //swap(RhopPre,Rhop)
  }
  if(TVisco==VISCO_LaminarSPS){ 
    memset(Tau,0,sizeof(tsymatrix3f)*(Np-Npb));
  }
  if(Nfloat)InitFloating();
  Part=PartIni; Ndiv=0; PartNdiv=0; PartOut=0;
  TimeStep=TimeStepIni; TimeStepM1=TimeStep;
}

//==============================================================================
/// Starts simulation on CPU.
//==============================================================================
void JSphCpu::Run(JCfgRun *cfg,JLog *log){
  const char* met="Run";
  char cad[256];
  if(!cfg||!log)return;
  Hardware="Cpu";

  //-Loads parameters and input data.
  //---------------------------------------------------------------------------
  JSph::Run(cfg,log);
  OmpMode=cfg->OmpMode;
  if(OmpMode==OMPM_Single)RunMode=string("Single core, Symmetry:")+(USE_SYMMETRY? "True": "False");
  else{
    RunMode=string("OpenMP threads:")+fun::IntStr(OmpThreads);
    RunMode=RunMode+", "+GetNameOmpMode(OmpMode);
  }
  RunMode=RunMode+", SSE:False";
  Log->Print(fun::VarStr("RunMode",RunMode));

  if(!CellModeIsCPU(CellMode))CellMode=(Simulate2D? CELLMODE_2H: CELLMODE_H); //-Valid CellMode.
  Hdiv=(CellMode==CELLMODE_2H? 1: 2);
  Log->Print(fun::VarStr("CellMode",string(GetNameCellMode(CellMode))));
  if(Hdiv!=1)Log->Print(fun::VarStr("Hdiv",Hdiv));


  //-Configures device and allocates memory.
  //-------------------------------------------------------------------------
  AllocMemory();
  DivideConfig();
  LoadPartBegin();
  InitVars();
 
  //-First Neighbour List computation and PART_0000.
  //---------------------------------------------------------------------------
  PrintMemoryAlloc();
  RunDivideBoundary();
  RunDivideFluid();
  SaveData(); PartNdiv=-1;
  Part++;

  //-MAIN LOOP.
  //---------------------------------------------------------------------------
  bool partoutstop=false;
  TimerSim.Start();
  TimerPart.Start();   //CkCalc=clock();  CkStartPart=CkCalc;
  sprintf(cad,"\n[Initialising simulation  %s]",fun::GetDateTime().c_str()); Log->Print(cad);
  PrintHeadPart();
  while(TimeStep<TimeMax){
    float stepdt=ComputeStep(Ndiv%DBCSteps==0);
    if(Nmoving)RunMotion(stepdt);
    ComputeStepDivide();
    if(ShepardSteps&&Ndiv&&( !((Ndiv+1)%ShepardSteps) || (TStep==STEP_Verlet && !(Ndiv%ShepardSteps)) ))RunShepard(); //-Shepard density filter.
    //-Using Verlet+Shepard, the density filter is applied at the given step and the following one to smooth Rhop and RhopM1.
    TimeStep+=stepdt;
    partoutstop=unsigned(Np-NpOk)>=PartOutMax;
    if((TimeStep-TimeStepIni)-TimePart*((Part-PartIni)-1)>=TimePart||SvSteps||partoutstop){
      if(partoutstop){
        Log->Print("\n**** Particles OUT limit reached...\n");
        TimeMax=TimeStep;
      }
      SaveData();
      Part++;
      PartNdiv=Ndiv;
      TimeStepM1=TimeStep;
      TimerPart.Start();
    }
    Ndiv++;
  }
  TimerSim.Stop(); TimerTot.Stop();

  //-End of simulation.
  //---------------------------------------------------------------------------
  sprintf(cad,"\n[Simulation %s  %s]",(partoutstop? "INTERRUPTED": "finished"),fun::GetDateTime().c_str());  Log->Print(cad);
  sprintf(cad,"DTs adjusted to DtMin............: %d",DtModif);  Log->Print(cad);
  sprintf(cad,"Excluded particles...............: %d",Np-NpOk);  Log->Print(cad);
  if(RhopOut){ sprintf(cad,"Excluded particles due to RhopOut: %u",RhopOutCount);  Log->Print(cad); }
  float tsim=TimerSim.GetElapsedTimeF()/1000.f;
  float ttot=TimerTot.GetElapsedTimeF()/1000.f;
  float tseg=tsim/TimeStep;
  float ndivseg=float(Ndiv)/tsim;
  sprintf(cad,"Total Runtime....................: %f sec.",ttot);    Log->Print(cad);
  sprintf(cad,"Simulation Runtime...............: %f sec.",tsim);    Log->Print(cad);
  sprintf(cad,"Time per second of simulation....: %f sec.",tseg);    Log->Print(cad);
  sprintf(cad,"Steps per second.................: %f",ndivseg);      Log->Print(cad);
  sprintf(cad,"Steps of simulation..............: %d",Ndiv);         Log->Print(cad);
  sprintf(cad,"PART files.......................: %d",Part-PartIni); Log->Print(cad);
  if(SvTimers)ShowTimers();
  if(SvRes)SaveRes(tsim,tseg,ttot,";RunMode",string(";")+RunMode);
  Log->Print(" ");
}

//==============================================================================
/// Shows active timers.
//==============================================================================
void JSphCpu::ShowTimers(){
  Log->Print("\n[CPU Timers]");
  for(unsigned c=0;c<TimerGetCount();c++){
    string tx=TimerToText(c); if(!tx.empty())Log->Print(tx);
  }
}

//==============================================================================
/// Creates files with output data.
//==============================================================================
void JSphCpu::SaveData(){
  const char met[]="SaveData";
  TmcStart(Timers,TMC_SuSavePart);
  //-Updates Pdata.
  unsigned nout=Div->GetOutCount();
  if(nout)Div->GetDataOut(Idp,Pos,Vel,Rhop,true);
  if(Nprobe)ComputeProbes();
  if(Pdata.SetDataUnsorted(Part,TimeStep,false,NpOk,nout,Idp,Pos,Vel,Rhop,ProbeVel,ProbeRhop))RunException(met,"Some excluded particles appear again in the simulation.");
  JSph::SaveData();
  TmcStop(Timers,TMC_SuSavePart);
}









