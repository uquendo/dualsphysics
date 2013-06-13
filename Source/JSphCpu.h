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

/// \file JSphCpu.h \brief Declares the class \ref JSphCpu

#ifndef _JSphCpu_
#define _JSphCpu_

#include "Types.h"
#include "JSphTimersCpu.h"
#include "JSph.h"
#include "JDivideCpu.h"
#include <string>
using namespace std;

//##############################################################################
//# JSphCpu
//##############################################################################
/// \brief Defines the attributes and functions to be used only in CPU simulations.

class JSphCpu : public JSph
{
protected:
  std::string RunMode;       ///<Stores execution mode (symmetry,openmp,balancing,...).

  class JDivideCpu *Div;     ///<Declares an object of the class \ref JDivideCpu
  unsigned BoundDivVer;      ///<Version of neighbour list (NL) to be compared with \ref BoundDatVer to know when the NL for boundaries must be computed.
  
  float ViscDt;              ///<Value to compute new time step according to viscosity terms dt2=H/(cs_max+H*ViscDt).
  float ViscDtThread[MAXTHREADS_OMP*STRIDE_OMP]; ///<Variable to record the ViscDt of each thread (even without OpenMP).

  TimersCpu Timers;          ///<Declares an array with timers for CPU (type structure \ref StSphTimerCpu).
  
  //-Resulting variables of force computation.                            
  tfloat3 *Ace;              ///<Acceleration of the particles (X,Y,Z).         [\ref INTER_Forces,\ref INTER_ForcesCorr,\ref INTER_ForcesKgc,\ref INTER_ForcesCorrKgc]
  float *Ar;                 ///<Density derivative.                            [\ref INTER_Forces,\ref INTER_ForcesCorr,\ref INTER_ForcesKgc,\ref INTER_ForcesCorrKgc]
  tfloat3 *VelXcor;          ///<XSPH correction of velocity (without Eps).     [\ref INTER_Forces,\ref INTER_ForcesKgc]
  float *FdWab;              ///<Kernel summation in Shepard Filter.            [\ref INTER_Shepard]
  float *FdRhop;             ///<Density summation in Shepard Filter.           [\ref INTER_Shepard]
  
  //-Variables used in the force computation.
  float *Csound;             ///<Speed of sound: Cs0*pow((Rhop[]*OVERRHOPZERO),3).
  float *PrRhop;             ///<Pressure term: Pres[]/(Rhop[]*Rhop[]).
  float *Tensil;             ///<Term for tensile correction with \ref KERNEL_Cubic.
  
  float *FdVol;              ///<Volume of the particle: Pm[]/Rhop[].           [\ref INTER_Shepard] 

  //-Variables used to update system using Verlet algorithm.
  tfloat3 *VelM1;            ///<Verlet: array to store velocity values of the previous time step.
  tfloat3 *VelNew;           ///<Verlet: array to store new velocity values.
  float *RhopM1;             ///<Verlet: array to store density values of the previous time step.
  float *RhopNew;            ///<Verlet: array to store new density values.
  int VerletStep;            ///<Current step of the Verlet algorithm after having applied Eulerian equations.
  byte VerletResetVel;       ///<Allows reinitialise velocity when change was performed the last two steps.

  //-Variables used to update system using Symplectic algorithm.
  tfloat3 *PosPre;           ///<Sympletic: array to store position values in predictor step.
  tfloat3 *VelPre;           ///<Sympletic: array to store velocity values in predictor step.
  float *RhopPre;            ///<Sympletic: array to store density values in predictor step.
  float DtPre;               ///<Sympletic: array to store time step value in predictor step.

  //-Variables for Laminar+SPS viscosity.  
  tsymatrix3f *Tau;          ///<SPS sub-particle stress tensor.
  tsymatrix3f *Csph;         ///<Velocity gradients.
 
  //-Variables for KGC.         
  tsymatrix3f *MatKgc;       ///<Matrix for kernel gradient correction. 

  //-Variables for floating bodies
  tfloat3* FtDist;           ///<Distance of the particles to the centre of the floating object.

  void AllocMemory();
  void Reset();
  void SaveData();
  void RunDivideBoundary();
  void RunDivideFluid();
  void ComputeStepDivide();

  void CallInteractionCells(TpInter tinter){
    switch(tinter){
      case INTER_Forces:        CallInteractionCells1<INTER_Forces>();         break;          //-Interaction to compute forces without correction in predictor step.
      case INTER_ForcesCorr:    CallInteractionCells1<INTER_ForcesCorr>();     break;          //-Interaction to compute forces without correction in corrector step (no XSPH is used).
      case INTER_ForcesKgc:     CallInteractionCells1<INTER_ForcesKgc>();      break;          //-Interaction to compute corrected forces in predictor step.
      case INTER_ForcesCorrKgc: CallInteractionCells1<INTER_ForcesCorrKgc>();  break;          //-Interaction to compute corrected forces in corrector step (no XSPH is used).
      case INTER_Shepard:       CallInteractionCells1<INTER_Shepard>();        break;          //-Interaction to compute density reinitialisation with Shepard filter.
      case INTER_MatrixKgc:     CallInteractionCells1<INTER_MatrixKgc>();      break;          //-Interaction to compute elements of the matrix for kernel gradient correction.
    }
  }
  template<TpInter tinter> void CallInteractionCells1(){
    switch(TKernel){
      case KERNEL_Cubic:     CallInteractionCells2<tinter,KERNEL_Cubic>();     break;          //-Cubic Spline kernel.
      case KERNEL_Wendland:  CallInteractionCells2<tinter,KERNEL_Wendland>();  break;          //-Wendland quintic kernel.
    }
  }
  template<TpInter tinter,TpKernel tker> void CallInteractionCells2(){
    switch(TVisco){
      case VISCO_Artificial: CallInteractionCells3<tinter,tker,VISCO_Artificial>();  break;    //-Artificial viscosity treatment.
      case VISCO_LaminarSPS: CallInteractionCells3<tinter,tker,VISCO_LaminarSPS>();  break;    //-Laminar viscosity and Sub-Partice Scale Turbulence.
    }
  }
  template<TpInter tinter,TpKernel tker,TpVisco tvis> void CallInteractionCells3(){
    if(Nfloat)InteractionCells<tinter,tker,tvis,true>();
    else InteractionCells<tinter,tker,tvis,false>();
  }

  template<TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void InteractionCells();
  template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void InteractionCells_Single();
#ifdef USE_OPENMP
  template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void InteractionCells_Dynamic();
  template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void InteractionCells_Static();
#endif

  template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void InteractSelf(int box,byte kind,byte kind2ini);
  template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> void InteractCelij(int ibegin,int iend,int box1,int box2,byte kind2ini);

  template<bool tsym,TpKernel tker> void ComputeForcesShepard(int i,int j,float drx,float dry,float drz,float rr2,int offsetsh);
  template<bool tsym,TpKernel tker> void ComputeMatrixKgc(TpParticle tpi,TpParticle tpj,int i,int j,float drx,float dry,float drz,float rr2,int offsetf);  
  template<bool tsym,TpInter tinter,TpKernel tker,TpVisco tvis,bool floating> float ComputeForces(TpParticle tpi,TpParticle tpj,int i,int j,float drx,float dry,float drz,float rr2,int offset,int offsetf); 

#ifdef USE_OPENMP
  void OmpMergeDataSum(int ini,int fin,float *data,int stride,int nthreads);
  void OmpMergeDataSum(int ini,int fin,tfloat3 *data,int stride,int nthreads);
  void OmpMergeDataSum(int ini,int fin,tsymatrix3f *data,int stride,int nthreads);
#endif
  
  float DtVariable(bool savedts);
  void ComputeRhop(float* rhopnew,const float* rhopold,float armul,bool rhopbound);
  void ComputeRhopEpsilon(float* rhopnew,const float* rhopold,float armul,bool rhopbound);

  void SPSCalcTau();
  void Interaction_MatrixKgc();
  void PreInteraction_Forces(TpInter tinter);
  template<TpKernel tker> void PreInteraction_Forces_(TpInter tinter);
  void Interaction_Forces(bool forcescorr);

  float ComputeStep(bool rhopbound){ return(TStep==STEP_Verlet? ComputeStep_Ver(rhopbound): ComputeStep_Sym(rhopbound)); }
  void ComputeVerletVars(const tfloat3 *vel1,const tfloat3 *vel2,float dt,float dt2,tfloat3 *velnew);
  float ComputeStep_Ver(bool rhopbound);
  float ComputeStep_Sym(bool rhopbound);

  void RunShepard();
  void RunMotion(float stepdt);
  void InitFloating();
  void RunFloating(float dt2,bool predictor);
  void DivideConfig();
  void InitVars();

public:
  JSphCpu();
  ~JSphCpu();
  void Run(JCfgRun *cfg,JLog *log);
  JDivideCpu* GetDivide(){ return(Div); };    
  tfloat3* GetPos(){ return(Pos); };

  void ShowTimers();
  unsigned GetMemoryCpu()const{ return(JSph::GetMemoryCpu()+Div->GetMemoryAlloc()); }

  unsigned TimerGetCount()const{ return(TmcGetCount()); }
  bool TimerIsActive(unsigned ct)const{ return(TmcIsActive(Timers,(CsTypeTimerCPU)ct)); }
  float TimerGetValue(unsigned ct)const{ return(TmcGetValue(Timers,(CsTypeTimerCPU)ct)); }
  string TimerGetName(unsigned ct)const{ return(TmcGetName((CsTypeTimerCPU)ct)); }
};

#endif







