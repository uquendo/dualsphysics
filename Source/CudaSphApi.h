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

/// \file CudaSphApi.h \brief Declares all the functions and variables used on GPU execution.

#ifndef _CudaSphApi_
#define _CudaSphApi_

#include "JSphTimersGpu.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "Types.h"
#include <ctime>
#include <cstdio>

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

inline float3 Float3(const tfloat3& v){ float3 p={v.x,v.y,v.z}; return(p); }
inline float3 Float3(float x,float y,float z){ float3 p={x,y,z}; return(p); }
inline tfloat3 TFloat3(const float3& v){ return(TFloat3(v.x,v.y,v.z)); }

#define BLOCKSIZE 256 //128 256 320
#define GRIDSIZEDIM 65535 

typedef enum{DV_Fluid,DV_All}CsTypeDivide;
typedef enum{IT_FluidFluid=10,IT_FluidBound=11,IT_BoundFluid=12}CsTypeInterac;

///Structure with variables to be stored in the constant memory of GPU.
typedef struct{
  float h;                   ///<Smoothing length (=coef*sqrt(dx*dx+dy*dy+dz*dz))
  float fourh2;              ///< \ref h * \ref h * 4 
  float cubic_a2,cubic_a24,cubic_c1,cubic_d1,cubic_c2,cubic_odwdeltap; ///<Ctes of Cubic Spline kernel.
  float wendland_awen,wendland_bwen;                                   ///<Ctes of Wendland kernel.
  float cteshepard;          ///<Constant used in Shepard filter to evaluate the self contribution of the particle itself.
  float cs0;                 ///<Speed of sound of reference.
  float visco;               ///<Viscosity value.
  float eta2;                ///<eta*eta being eta=0.1*\ref h
  float eps;                 ///<Epsilon constant for XSPH variant.
  float massb;               ///<Mass of a boundary particle.
  float massf;               ///<Mass of a fluid particle.
  unsigned nbound;           ///<Number of boundary particles.
}StDeviceCte; 

///Structure with variables to be used in the CUDA files.
typedef struct{
  //-Constants on GPU.
  int simulate2d;            ///<Activates or deactivates 2D simulation (Y-forces are cancelled).
  float cflnumber;           ///<Constant for the Courant condition (0.1 - 0.5)
  unsigned hdiv;             ///<Value to divide 2h.
  unsigned ncellnv;          ///<Maximum number of ranges of neighbours of a cell.((\ref hdiv*2)+1)^2. i.e. 2h:9,h:25,2h/3:49
  float ovscell;             ///<One over cells size. 
  float3 posmin;             ///<Lower limit of the domain.
  float3 difmax;             ///<Upper limit of the domain  (with \ref posmin subtracted)
  int ncx;                   ///<Number of cells in X direction.
  int ncy;                   ///<Number of cells in Y direction.
  float gamma;               ///<Politropic constant. (1-7).
  float b;                   ///<Constant that sets a limit for the maximum change in density.
  float3 gravity;            ///<Gravity acceleration.

  unsigned np;               ///<Number of total particles.  
  unsigned nfixed;           ///<Number of fixed boundary particles.
  unsigned nmoving;          ///<Number of moving boundary particles.
  unsigned nfloat;           ///<Number of floating boundary particles. 
  unsigned nbound;           ///<Number of boundary particles ( \ref nfixed + \ref nmoving + \ref nfloat ).
  unsigned nfluid;           ///<Number of fluid particles (including the excluded ones).   

  //-Variables on GPU.
  unsigned npok;             ///<Number of total particles activated at each time step (\ref np - excluded particles). 
  unsigned npb;              ///<Number of particles of the boundary block ( \ref nbound - \ref nfloat ) or ( \ref nfixed + \ref nmoving).
  unsigned npf;              ///<Number of particles treated as fluid no excluded (\ref nfloat + \ref nfluid - \ref nfluidout).

  unsigned memcpu;           ///<Number of bytes of allocated memory in Cpu.
  unsigned memgpu;           ///<Number of bytes of allocated memory in Gpu (basic). 
  unsigned memgpu2;          ///<Number of bytes of allocated memory in Gpu (depends on CellMode).
  
  unsigned ncz;              ///<Number of cells in Z direction.
  unsigned nct;              ///<Number of total cells.
  unsigned nctotmax;         ///<Number of maximum total cells
  unsigned nczbound;         ///<Number of cells in Z direction for boundaries.
  unsigned ndiv;             ///<Current time step.

  float massbound;           ///<Mass of a boundary particle.
  float massfluid;           ///<Mass of a fluid particle.

  //-Arrays to manage excluded particles (out) allocated in CPU.
  unsigned outsize;          ///<Number of particles that can be allocated.
  unsigned outcount;         ///<Current stored number of particles excluded.
  unsigned *outidp;          ///<Idp of current particles out.
  float3 *outpos;            ///<Position of current particles out.
  float3 *outvel;            ///<Velocity of current particles out.
  float *outrhop;            ///<Density of current particles out.

  //-Features of the GPU.
  size_t mglobal;            ///<Size of global memory in bytes.
  unsigned mshared;          ///<Size of shared memory per block in bytes.
  unsigned compute;          ///<Compute capbility: 10,11,12,20

  //-Variables of control of the execution.
  TpStep tstep;              ///<Type of step algorithm.    
  TpKernel tkernel;          ///<Type of kernel function. 
  TpVisco tvisco;            ///<Type of viscosity treatment.
  bool kgc;                  ///<Activates Kernel Gradient Correction. 

  unsigned verletstep;       ///<Current step of the Verlet algorithm after having applied Eulerian equations.
  unsigned verletsteps;      ///<Verlet: number of steps to apply Eulerian equations (def=40).
  unsigned verletresetvel;   ///<Allows reinitialise velocity when change was performed the last two steps.

  unsigned shepardsteps;     ///<Number of steps to apply Shepard density filter.
  unsigned dbcsteps;         ///<Number of steps to update the density of the boundaries.
  unsigned bounddatver;      ///<Version of data to register when movement of boundaries must be computed.
  unsigned bounddivver;      ///<Version of neighbour list (NL) to be compared with \ref BoundDatVer to know when the NL for boundaries must be computed.
   
  unsigned lastdivall;       ///<Indicates the type of the last NL (1==DV_All, 0==DV_Fluid).

  //-Particle data [np].
  unsigned *idp;             ///<Particle identifier according to its position in data.
  float3 *pos;               ///<Position of the particles (X,Y,Z).
  float3 *vel;               ///<Velocity of the particles (X,Y,Z).
  float *rhop;               ///<Density of the particles.
           
  //-Position in data of particle according to idp.
  unsigned *ridpmv;          ///<Position in data of moving boundary particles [\ref nmoving] and when \ref nmoving!=0. 
  unsigned *ridpft;          ///<Position in data of floating particles [\ref nfloat] and when \ref nfloat!=0.

  //-Variables for force computation.
  float3 *ace;               ///<Acceleration of the particles (X,Y,Z). [\ref np]
  float *ar;                 ///<Density derivative.
  float3 *velxcor;           ///<XSPH correction of velocity (without Eps).
  float *viscdt;             ///<Value to compute new time step according to viscosity terms dt2=h/(cs_max+h*viscdt).
  float *csound;             ///<Speed of sound: cs0*pow((rhop[]*OVERRHOPZERO),3).

  //- VERLET
  float3 *velm1,*velm1_2;    ///<Verlet: array to store velocity values of the previous time step.
  float3 *velnew;            ///<Verlet: array to store new velocity values.
  float *rhopm1;             ///<Verlet: array to store density values of the previous time step.
  float *rhopnew;            ///<Verlet: array to store new density values.
  
  //-SYMPLECTIC
  float3 *pospre,*pospre_2;  ///<Sympletic: array to store position values in predictor step.
  float3 *velpre,*velpre_2;  ///<Sympletic: array to store velocity values in predictor step.
  float *rhoppre;            ///<Sympletic: array to store density values in predictor step.
  float dtpre;               ///<Sympletic: array to store time step value in predictor step.

  float4 *pospres;           ///<New array that combines position and pressure (pospres[p]=make_float4(r.x,r.y,r.z,press)).
  float4 *velrhop;           ///<New array that combines velocity and density (velrhop[p]=make_float4(r.x,r.y,r.z,rrhop)).

  //-Variables for Laminar+SPS viscosity.  
  float smag;                ///<Smagorinsky constant used in SPS turbulence model.
  float blin;                ///<Blin constant used in the SPS turbulence model.
  tsymatrix3f *tau;          ///<SPS sub-particle stress tensor.
  tsymatrix3f *csph;         ///<Velocity gradients.

  //-Variables for KGC.  
  tsymatrix3f *matkgc;       ///<Matrix for kernel gradient correction. 
  
  //-Variables for floating bodies
  StFloatingData* ftobjs;    ///<Data of floating object.
  unsigned ftcount;          ///<Number of floating objects.
  float3 *ftdist;            ///<Distance of the particles to the centre of a floating object.
  float3 *face;              ///<Aceleration of the centre of a floating object.
  float3 *fomegavel;         ///<Angular velocity of the centre of a floating object.

  //-Vars. para density filter: SHEPARD
  float *fdrhop;             ///<Density summation in Shepard Filter.  

  float *reduaux;            ///<Array that helps in the reduction of \ref np values.

  //-Arrays for reordering in GPU.
  unsigned *idsortb;         ///<Position in data before reordering [\ref nb].
  unsigned *idsortf;         ///<Position in data before reordering [\ref nf].
  unsigned *cellb;           ///<Cell of the boundary particles after reordering [\ref nb].
  unsigned *cellf;           ///<Cell of the fluid particles after reordering [\ref nf].
  int2 *cellbegb;            ///<First and last boundary particle of each cell [\ref nct+2].
  int2 *cellbegf;            ///<First and last fluid particle of each cell [\ref nct+2].

  unsigned use_cellnv;       ///<Indicates if \ref cellnvb and \ref cellnvf is used or not.
  uint2 *cellnvb;            ///<First and last boundary particle of each neighbouring cell [(\ref nctotmax+1)*9].
  uint2 *cellnvf;            ///<First and last fluid particle of each neighbouring cell [(\ref nctotmax+1)*9].

  float dtmin;               ///<Minimum time step.
  int dtmodif;               ///<Number of modifications performed when the new value of dt is too low.
  StInfoDt *infodt;          ///<Structure to monitor the used dt values.

  TimersGpu timers;          ///<Declares an array with timers for GPU (type structure \ref StSphTimerGpu).

  unsigned rhopout;          ///<Indicates whether the density correction rhopOut is activated or not.
  float rhopoutmin;          ///<Minimum value for the density correction rhopOut.
  float rhopoutmax;          ///<Minimum value for the density correction rhopOut..
  unsigned rhopoutcount;     ///<Number of excluded particles with the density correction RhopOut.

  TpCellMode cellmode;       ///<Modes of cells division.

  unsigned bsforcesfluid;    ///<Block size used in the interaction kernels that compute forces of fluid particles.
  unsigned bsforcesfluidcorr;///<Block size used in the interaction kernels that compute forces of fluid particles in corrector step.
  unsigned bsforcesbound;    ///<Block size used in the interaction kernels that compute forces of boundary particles.
  unsigned bsshepard;        ///<Block size used in the interaction kernels that apply Shepard filter.
  unsigned bsmatrixkgc;      ///<Block size used in the interaction kernels that compute elements of matrix KGC.

}StDeviceContext; 

//-Set of functions called from the JSphGpu.cpp file and implemented in the CUDA files (.cu).
int CsInitCuda(int gpuid,int *devdef,char *devname,StDeviceContext *dc);
void CsInit(StDeviceContext *dc);
void CsReset(StDeviceContext *dc);
void CsAllocMemoryBasic(StDeviceContext *dc,StDeviceCte *cte,bool svtimers);
bool CsAllocMemoryCellMode(StDeviceContext *dc,StDeviceCte *cte);
void CsUpCte(StDeviceCte *cte);
void CsUpData(StDeviceContext *dc,StDeviceCte *cte,unsigned *idp,float3 *pos,float3 *vel,float *rhop);
void CsDownData(StDeviceContext *dc,unsigned pini,unsigned *idp,float3 *pos,float3 *vel,float *rhop);
void CsDownDataRhop(StDeviceContext *dc,float *rhop);
void CsCallDivide(CsTypeDivide tdiv,StDeviceContext *dc,StDeviceCte *cte);
void CsOutGetData(StDeviceContext *dc,unsigned *idp,float3 *pos,float3 *vel,float *rhop);
void CsCallComputeStepDivide(StDeviceContext *dc,StDeviceCte *cte);
float CsCallComputeStep(StDeviceContext *dc,StDeviceCte *cte,bool rhopbound);
void CsSymplecticMotionPre(StDeviceContext *dc,bool all);
void CsSymplecticMotionPost(StDeviceContext *dc);
void CsCalcRidpmv(StDeviceContext *dc);
void CsMoveLinBound(StDeviceContext *dc,unsigned pini,unsigned npar,float3 mvsimple,float3 mvvel);
void CsMoveMatBound(StDeviceContext *dc,unsigned pini,unsigned npar,tmatrix4f mvmatrix,float dt);
void CsCallRunShepard(StDeviceContext *dc,StDeviceCte *cte);

#endif






