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

/// \file Types.h \brief Defines specific types for the SPH application.

#ifndef _Types_
#define _Types_

#include "TypesDef.h"

#define USE_OPENMP          ///<Enables/Disables OpenMP.   
#define MAXTHREADS_OMP 32
#define STRIDE_OMP 200

//#define USE_SSE             ///<Enables/Disables the use of intruccions SSE (SIMD).

#define USE_SYMMETRY true   ///<Enables/Disables the symmetry in the force computation.

#define OVERPI 0.318309886  ///<Value of 1/PI.

#define RHOPZERO 1000.f
#define OVERRHOPZERO 0.001f

#define MATRIXKGC_CONTROL 0.25f  ///<To prevent forces will increase more than 4 times in KGC matrix. 

#define MK_RANGE 256        ///<Maximum range to label particles.

/// Structure with the information of the floating object.
typedef struct{
  unsigned begin;
  unsigned count;
  float mass;
  tfloat3 inertia;
  tfloat3 center;
  tfloat3 fvel;
  tfloat3 fomega;
}StFloatingData;

#define INFODT_MAX 200      ///<Maximum number of values used in the structure \ref StInfoDt

/// Structure with the information of the average value of the time step (dt).
typedef struct {
  float vdt[INFODT_MAX];    ///<Multiple values of dt can be stored.
  int count;                ///<Number of  values in dt.
  float mean;               ///<Average value of dt.
  int meancount;            ///<Number of values in dtmean.
}StInfoDtMean; 

/// Structure with the information of the time step (dt).
typedef struct {
  StInfoDtMean dt;          ///<Average value of dtstep.
  StInfoDtMean dt1;         ///<Average value of dt1.
  StInfoDtMean dt2;         ///<Average value of dt2.
  float dtmin;              ///<Minimum value of dt.
  float dtmax;              ///<Maximum value of dt.
  unsigned outcount;        ///<Number of excluded particles.
  unsigned rhopoutcount;    ///<Number of excluded particles by RhopOut.
}StInfoDt;

///Controls the output of information on the screen and/or log.
typedef enum{ 
    MOUT_ScrFile=3,
    MOUT_File=2,
    MOUT_Screen=1,
    MOUT_None=0
}TpModeOut;                 

///Options of the format of output files.
typedef enum{ 
    SDAT_Binx2=1,           ///<Bynary format .bi2
    SDAT_Vtk=2,             ///<Vtk format .vtk
    SDAT_Csv=4,             ///<Csv format .csv
    SDAT_Sphysics=8,        ///<Ascii format (SPHysics code).
    SDAT_Flw=16,            ///<Flw format (Poligon code).
    SDAT_None=0            
}TpSaveDat;                 

///Types of step algorithm.
typedef enum{ 
    STEP_Symplectic=2,      ///<Symplectic algorithm.
    STEP_Verlet=1,          ///<Verlet algorithm.
    STEP_None=0 
}TpStep;                    

///Types of kernel function.
typedef enum{ 
    KERNEL_Wendland=2,      ///<Wendland kernel.
    KERNEL_Cubic=1,         ///<Cubic Spline kernel.
    KERNEL_None=0 
}TpKernel;                  

///Types of viscosity treatment.
typedef enum{ 
    VISCO_LaminarSPS=2,     ///<Laminar viscosity and Sub-Partice Scale Turbulence.
    VISCO_Artificial=1,     ///<Artificial viscosity.
    VISCO_None=0 
}TpVisco;                   

///Types of interaction.
typedef enum{ 
    INTER_ForcesCorrKgc=6,  ///<Interaction to compute forces with the corrected kernel gradient using the corrector step of Symplectic algorithm where XSPH variant is not applied. 
    INTER_ForcesKgc=5,      ///<Interaction to compute forces with the corrected kernel gradient using the Verlet algorithm and the predictor step of Symplectic algorithm. 
    INTER_MatrixKgc=4,      ///<Interaction to compute the components of the matrix when using kernel gradient correction.
    INTER_Shepard=3,        ///<Interaction to compute the new density values when using Shepard density filter.
    INTER_ForcesCorr=2,     ///<Interaction to compute forces using the corrector step of Symplectic algorithm where XSPH variant is not applied.
    INTER_Forces=1          ///<Interaction to compute forces using the Verlet algorithm and the predictor step of Symplectic algorithm. 
}TpInter;                   

///Types of particles.
typedef enum{ 
    PART_BoundFx=1,         ///<Fixed boundary particles.
    PART_BoundMv=2,         ///<Moving boundary particles.
    PART_BoundFx_BoundMv=3, ///<Both fixed and moving boundary particles.
    PART_BoundFt=4,         ///<Floating boundary particles.
    PART_Fluid=8,           ///<Fluid particles.
    PART_BoundFt_Fluid=12   ///<Both floating and fluid particles.
}TpParticle;                

///Types of execution with or without OpenMP.
typedef enum{ 
    OMPM_Single,            ///<Single execution with one core of CPU without OpenMP.
    OMPM_Dynamic,           ///<Multiple-core execution using OpenMP with dynamic load balancing.
    OMPM_Static             ///<Multiple-core execution using OpenMP with staticc load balancing.
}TpOmpMode;                 

///Returns the name of the type of interaction with or without OpenMP in text format.
inline const char* GetNameOmpMode(TpOmpMode mode){
  switch(mode){
    case OMPM_Single:    return("Single");
    case OMPM_Dynamic:   return("Dynamic");
    case OMPM_Static:    return("Static");
  }
  return("???");
}

///Modes of cells division.
typedef enum{ 
   CELLMODE_None=0
  ,CELLMODE_2H=1            ///<Cells of size 2h.
  ,CELLMODE_H=2             ///<Cells of size h.
  ,CELLMODE_Hneigs=3        ///<Cells of size h using a list of 25 ranges per cell.
}TpCellMode; 

///Returns the name of the cellmode in text format.
inline const char* GetNameCellMode(TpCellMode cellmode){
  switch(cellmode){
    case CELLMODE_2H:      return("2H");
    case CELLMODE_H:       return("H");
    case CELLMODE_Hneigs:  return("H-neigs");
  }
  return("???");
}

inline bool CellModeIsGPU(TpCellMode m){
  return(m==CELLMODE_2H||m==CELLMODE_H||m==CELLMODE_Hneigs);
}
inline bool CellModeIsCPU(TpCellMode m){
  return(m==CELLMODE_2H||m==CELLMODE_H);
}


#endif






