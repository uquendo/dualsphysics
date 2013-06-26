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

/// \file JSph.h \brief Declares the class \ref JSph

#ifndef _JSph_
#define _JSph_

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

#include "Types.h"
#include "JObject.h"
#include "JCfgRun.h"
#include "JLog.h"
#include "JPartData.h"
#include "JTimer.h"

#include <float.h>
#include <string>
#include <cmath>
#include <ctime>
#include <sstream>
#include <iostream>
#include <fstream>
using namespace std;

#ifdef USE_OPENMP
  #include <omp.h>  
#endif

class JSphMotion;

//##############################################################################
//# JSph
//##############################################################################
/// \brief Defines all the attributes and functions that CPU and GPU simulations share.

class JSph : protected JObject
{
public:
  static const unsigned int VersionMajor=80;
  static const unsigned int VersionMinor=0;
  static std::string GetVersionStr();

/// Structure with constants for the Cubic Spline Kernel.
  typedef struct {
    float a1,a2,aa,a24,c1,d1,c2;
    float od_wdeltap;        ///<Parameter for tensile instability correction.  
  }StCubicCte;

/// Structure with constants for the Wendland Kernel.
  typedef struct {
    float awen,bwen;
  }StWendlandCte;

private:
  bool RidpbOld;             ///<Indicates whether Ridp for Npb particles (boundaries)  is updated or not.
  bool RidpfOld;             ///<Indicates whether Ridp for Npf particles (fluid ones) is updated or not.

protected:
  bool Simulate2D;           ///<Activates or deactivates 2D simulation (Y-forces are cancelled).
  string Hardware;
  bool Cpu;                  ///<Activates or deactivates only CPU simulation.
  int OmpThreads;            ///<Maximum number of threads for OpenMP CPU executions.
  TpOmpMode OmpMode;         ///<Type of execution with or without OpenMP.
  TpCellMode CellMode;       ///<Modes of cells division.
  unsigned Hdiv;             ///<Value to divide 2h.

  TpStep TStep;              ///<Type of step algorithm. 
  TpKernel TKernel;          ///<Type of kernel function. 
  StCubicCte CubicCte;       ///<Constants for Cubic Spline Kernel.
  StWendlandCte WendlandCte; ///<Constants for Wendland Kernel.
  TpVisco TVisco;            ///<Type of viscosity treatment.
  JLog *Log;                 ///<Declares an object of the class \ref JLog.
  int VerletSteps;           ///<Verlet: number of steps to apply Eulerian equations (def=40).
  bool Kgc;                  ///<Activates Kernel Gradient Correction.  
  int ShepardSteps;          ///<Number of steps to apply Shepard density filter.
  int DBCSteps;              ///<Number of steps to update the density of the boundaries.
  string CaseName;           ///<Name of the case.
  string DirCase;            ///<Name of the input directory.
  string DirOut;             ///<Name of the output directory.
  string RunName;            ///<Name of the execution case.

  string PartBeginDir;       ///<Directory to look for the PART file to start.
  unsigned PartBegin;        ///<Indicates the PART file to start.
  unsigned PartBeginFirst;   ///<Indicates the number of the firest PART file to be created.

  bool CaseLoaded;           ///<Activated once the case has been successfully loaded.

  unsigned Np;               ///<Number of total particles.  
  unsigned Nfixed;           ///<Number of fixed boundary particles. 
  unsigned Nmoving;          ///<Number of moving boundary particles. 
  unsigned Nfloat;           ///<Number of floating boundary particles. 
  unsigned Nbound;           ///<Number of boundary particles ( \ref Nfixed + \ref Nmoving + \ref Nfloat ).
  unsigned Npb;              ///<Number of particles of the boundary block ( \ref Nbound - \ref Nfloat ) or ( \ref Nfixed + \ref Nmoving).
  unsigned Nfluid;           ///<Number of fluid particles (including the excluded ones). 
  unsigned NpOk;             ///<Number of total particles activated at each time step (\ref Np - excluded particles). 
  unsigned Nprobe;      ///<Number of probe particles.
 
  float H;                   ///<Smoothing length (=coef*sqrt(dx*dx+dy*dy+dz*dz))
  float CteB;                ///<Constant that sets a limit for the maximum change in density.
  float Gamma;               ///<Politropic constant. (1-7).
  float Eps;                 ///<Epsilon constant for XSPH variant.
  float Visco;               ///<Viscosity value.
  float CFLnumber;           ///<Constant for the Courant condition (0.1 - 0.5).
  float Cs0;                 ///<Speed of sound of reference.
  float Dp;                  ///<Inter-particle distance.
  float MassBound;           ///<Mass of a boundary particle.
  float MassFluid;           ///<Mass of a fluid particle.
  float TimeMax;             ///<Time of simulation.
  float TimePart;            ///<Time between output files.
  
  tfloat3 Gravity;           ///<Gravity acceleration.
  float DtIni;               ///<Initial time step.
  float DtMin;               ///<Minimum time step.
  float Dosh;                ///<2*\ref H.
  float H2;                  ///<\ref H*\ref H.
  float Fourh2;              ///<\ref H2*\ref H2.
  float Eta2;                ///<eta*eta being eta=0.1*\ref H
  float CteShepard;          ///<Constant used in Shepard filter to evaluate the self contribution of the particle itself.
  float Smag;                ///<Smagorinsky constant used in SPS turbulence model.
  float Blin;                ///<Blin constant used in the SPS turbulence model.
  float IncZ;                ///<Allowed increase in Z direction distZ=(distZ*(1+IncZ)).
  unsigned PartOutMax;       ///<Allowed percentage of fluid particles out the domain.

  StInfoDt InfoDt;           ///<Structure to monitor the used dt values.
  int DtModif;               ///<Number of modifications performed when the new value of dt is too low.

  int SvData;               ///<Indicates the format of the output files.               

  ofstream *DtPf;            ///<Pointer for the file with info of DT.
  bool SvDt;                 ///<Stores a file with info of DT.
  bool SvRes;                ///<Stores a file with info of the execution.  
  bool SvTimers;             ///<Computes the time for each process.
  
  bool SvSteps;              ///<Stores the output data of all time steps. 
  JPartData Pdata;           ///<Declares an object of the class \ref JPartData.

  unsigned *Idp;             ///<Particle identifier according to its position in data.
  unsigned *Ridp;            ///<Position in data according to the particle identifier.
  tfloat3 *Pos;              ///<Position of the particles (X,Y,Z).
  tfloat3 *Vel;              ///<Velocity of the particles (X,Y,Z).
  float *Rhop;               ///<Density of the particles.

  tfloat3 *ProbePos;         ///<Position of the probe particles (X,Y,Z).
  tfloat3 *ProbeVel;         ///<Velocity of the probe particles (X,Y,Z).
  float *ProbeRhop;          ///<Density of the probe particles.

  tint3 Ncells;              ///<Number of cells for each direction.
  int Nct;                   ///<Number of total cells.
  int Nsheet;                ///<Number of cells in X direction * Y direction.
  int PartIni;               ///<First PART file.
  int Part;                  ///<Next PART to be stored.
  int Ndiv;                  ///<Current time step.
  int PartNdiv;              ///<Number of steps when the last PART was stored.           
  int PartOut;               ///<Number of particulas out (excluded) when the last PART was stored.

  float TimeStepIni;         ///<Initial instant of the simulation.
  float TimeStep;            ///<Current instant of the simulation.
  float TimeStepM1;          ///<Instant of the simulation when the last PART was stored.
  JTimer TimerTot;           ///<Total runtime of the simulation.
  JTimer TimerSim;           ///<Total runtime starting when the computation of the main loop starts.
  JTimer TimerPart;          ///<Runtime since the last PART to the next one.

  unsigned MemCpuStatic;     ///<Number of bytes of all the allocated memory.

  unsigned BoundDatVer;      ///<Version of data to register when movement of boundaries must be computed.
  JSphMotion* Motion;        ///<Declares an object of the class \ref JSphMotion.
  float MotionTimeMod;       ///<Modifier of \ref TimeStep for \ref Motion.
  unsigned MotionObjCount;   ///<Number of objects with motion.
  unsigned MotionObjBegin[256];

  StFloatingData* FtObjs;    ///<Data of floating object.
  unsigned FtCount;          ///<Number of floating objects.

  bool RhopOut;              ///<Indicates whether the density correction RhopOut is activated or not.
  float RhopOutMin;          ///<Minimum value for the density correction RhopOut.
  float RhopOutMax;          ///<Minimum value for the density correction RhopOut..
  unsigned RhopOutCount;     ///<Number of excluded particles with the density correction RhopOut.
  
  void Reset();
  void AllocMemory(int np);
  void AllocMemoryFloating(unsigned ftcount);
  void AllocMemoryProbes(unsigned np);
  void LoadCase();
  void LoadPartBegin();

  void SaveData();
  void PrintMemoryAlloc();
  void PrintHeadPart();

  void InfoDtReset();
  void InfoDtMean(StInfoDtMean &sdt);
  void InfoDtAdd(float dt,StInfoDtMean &sdt);
  void InfoDtsAdd(float dt1,float dt2);  
  void InfoDtStepAdd(float dtstep);
  void SaveDt();

  void SaveRes(float tsim,float tseg,float ttot,const string &headplus="",const string &detplus="");

  void ClearRidpBound(){ RidpbOld=true; }
  void ClearRidpFluid(){ RidpfOld=true; }
  void ClearRidp(){ ClearRidpBound(); ClearRidpFluid(); }  
  
  void UpdateRidpBound(){ if(RidpbOld){ for(unsigned p=0;p<Npb;p++)Ridp[Idp[p]]=p; RidpbOld=false; } }
  void UpdateRidpFluid(){ if(RidpfOld){ for(unsigned p=Npb;p<NpOk;p++)Ridp[Idp[p]]=p; RidpfOld=false; } }
  void UpdateRidp(){ UpdateRidpBound(); UpdateRidpFluid(); }   

public:

  JSph();
  ~JSph();
  void Run(JCfgRun *cfg,JLog *log);
  string GetDirOut();
  int GetPart();
  int GetNdiv();
  float GetH(){ return(H); }
  float GetDp(){ return(Dp); }
  int GetNp(){ return(Np); }
  int GetNbound(){ return(Nbound); }
  TpParticle GetTpParticleId(unsigned id)const{ return(id<Nbound? (id<Nfixed? PART_BoundFx: (id<Npb? PART_BoundMv: PART_BoundFt) ): PART_Fluid); }
  TpParticle GetTpParticlePos(unsigned pos)const{ return(Nfixed==Nbound? (pos<Nbound? PART_BoundFx: PART_Fluid): GetTpParticleId(Idp[pos])); }
  void VisuConfig();
  static string GetStepName(TpStep tstep);
  static string GetKernelName(TpKernel tkernel);
  static string GetViscoName(TpVisco tvisco);
  virtual unsigned GetMemoryCpu()const{ return(MemCpuStatic+Pdata.GetMemoryAlloc()); }
  virtual unsigned GetMemoryGpu()const{ return(0); }
  string TimerToText(const string &name,float value)const;
  string TimerToText(unsigned ct)const;
  virtual unsigned TimerGetCount()const{ return(0); }
  virtual bool TimerIsActive(unsigned ct)const{ return(false); }
  virtual float TimerGetValue(unsigned ct)const{ return(0); }
  virtual string TimerGetName(unsigned ct)const{ return("???"); }
};


#endif







