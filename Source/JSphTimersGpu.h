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

/// \file JSphTimersGpu.h \brief Measures time intervals during GPU execution.

#ifndef _JSphTimersGpu_
#define _JSphTimersGpu_

//#define TmcStart(x,y) ;
//#define TmcStop(x,y) ;
#define TmgStart(x,y) _TmgStart(x,y)
#define TmgStop(x,y) _TmgStop(x,y)

#include "JTimerCuda.h"

/// Structure with information of the timer and time value in GPU.
typedef struct{
  JTimerCuda timer;
  bool active;
  float time;
}StSphTimerGpu; 

typedef enum{
   TMG_NlLimits=0
  ,TMG_NlPreSort=1
  ,TMG_NlRadixSort=2
  ,TMG_NlSortData=3
  ,TMG_NlCellBegin=4
  ,TMG_NlOutCheck=5
  ,TMG_NlCellNv=6
  ,TMG_CfMatrixKgc=7
  ,TMG_CfPreForces=8
  ,TMG_CfForcesFluid=9
  ,TMG_CfForcesBound=10
  ,TMG_CfShepard=11
  ,TMG_SuCalcDt=12
  ,TMG_SuComputeStep=13
  ,TMG_SuFloating=14
  ,TMG_SuMotion=15
  ,TMG_SuDownData=16
  ,TMG_SuSavePart=17
}CsTypeTimerGPU;
#define TMG_COUNT 18
typedef StSphTimerGpu TimersGpu[TMG_COUNT];

//==============================================================================
/// Returns the name of the timer.
//==============================================================================
inline const char* TmgGetName(CsTypeTimerGPU ct){
  switch(ct){
    case TMG_NlLimits:         return("NL-Limits");
    case TMG_NlPreSort:        return("NL-PreSort");
    case TMG_NlRadixSort:      return("NL-RadixSort");
    case TMG_NlSortData:       return("NL-SortData");
    case TMG_NlCellBegin:      return("NL-CellBegin");
    case TMG_NlCellNv:         return("NL-CellNv");
    case TMG_NlOutCheck:       return("NL-OutCheck");
    case TMG_CfMatrixKgc:      return("CF-MatrixKgc");
    case TMG_CfPreForces:      return("CF-PreForces");
    case TMG_CfForcesFluid:    return("CF-ForcesFluid");
    case TMG_CfForcesBound:    return("CF-ForcesBound");
    case TMG_CfShepard:        return("CF-Shepard");
    case TMG_SuCalcDt:         return("SU-CalcDt");
    case TMG_SuComputeStep:    return("SU-ComputeStep");
    case TMG_SuFloating:       return("SU-Floating");
    case TMG_SuMotion:         return("SU-Motion");
    case TMG_SuDownData:       return("SU-DownData");
    case TMG_SuSavePart:       return("SU-SavePart");
  }
  return("???");
}
//==============================================================================
/// Returns the number of timers.
//==============================================================================
inline unsigned TmgGetCount(){ return(TMG_COUNT); }

//==============================================================================
/// Creates timers to measure time intervals.
//==============================================================================
inline void TmgCreation(TimersGpu vtimer,bool active){
  for(unsigned c=0;c<TMG_COUNT;c++){
    StSphTimerGpu* t=vtimer+c;
    t->timer.Reset();
    t->active=active;
    t->time=0;
  }
}

//==============================================================================
/// Destroys the timers used to measure times.
//==============================================================================
inline void TmgDestruction(TimersGpu vtimer){ TmgCreation(vtimer,false); }

//==============================================================================
/// Marks start of timer.
//==============================================================================
inline void _TmgStart(TimersGpu vtimer,CsTypeTimerGPU ct){ if(vtimer[ct].active)vtimer[ct].timer.Start(); }

//==============================================================================
/// Marks end of timer and accumulates time.
//==============================================================================
inline void _TmgStop(TimersGpu vtimer,CsTypeTimerGPU ct){
  StSphTimerGpu* t=vtimer+unsigned(ct);
  if(t->active){
    t->timer.Stop();
    t->time+=t->timer.GetElapsedTimeF();
  }
}

//==============================================================================
/// Returns the time accumulated by the timer.
//==============================================================================
inline float TmgGetValue(const TimersGpu vtimer,CsTypeTimerGPU ct){ return(vtimer[ct].time); }


//==============================================================================
/// Indicates whether a timer is active or not.
//==============================================================================
inline bool TmgIsActive(const TimersGpu vtimer,CsTypeTimerGPU ct){ return(unsigned(ct)<TmgGetCount()&&vtimer[ct].active); }

#endif






