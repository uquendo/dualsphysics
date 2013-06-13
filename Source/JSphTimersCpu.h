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

/// \file JSphTimersCpu.h \brief Measures time intervals during CPU execution.

#ifndef _JSphTimersCpu_
#define _JSphTimersCpu_

//#define TmcStart(x,y) ;
//#define TmcStop(x,y) ;
#define TmcStart(x,y) _TmcStart(x,y)
#define TmcStop(x,y) _TmcStop(x,y)

#include "JTimer.h" 

/// Structure with information of the timer and time value in CPU.
typedef struct{
  JTimer timer; 
  bool active;
  float time;
}StSphTimerCpu; 

typedef enum{
   TMC_NlDivideBound=0
  ,TMC_NlSortDataBound=1
  ,TMC_NlDivideFluid=2
  ,TMC_NlSortDataFluid=3
  ,TMC_CfMatrixKgc=4
  ,TMC_CfPreForces=5
  ,TMC_CfForces=6
  ,TMC_CfShepard=7
  ,TMC_SuCalcDt=8
  ,TMC_SuComputeStep=9
  ,TMC_SuFloating=10
  ,TMC_SuMotion=11
  ,TMC_SuSavePart=12
}CsTypeTimerCPU;

#define TMC_COUNT 13

typedef StSphTimerCpu TimersCpu[TMC_COUNT];

//==============================================================================
/// Returns the name of the timer.
//==============================================================================
inline const char* TmcGetName(CsTypeTimerCPU ct){
  switch(ct){
    case TMC_NlDivideBound:    return("NL-DivideBound");
    case TMC_NlSortDataBound:  return("NL-SortDataBound");
    case TMC_NlDivideFluid:    return("NL-DivideFluid");
    case TMC_NlSortDataFluid:  return("NL-SortDataFluid");
    case TMC_CfMatrixKgc:      return("CF-MatrixKgc");
    case TMC_CfPreForces:      return("CF-PreForces");
    case TMC_CfForces:         return("CF-Forces");
    case TMC_CfShepard:        return("CF-Shepard");
    case TMC_SuCalcDt:         return("SU-CalcDt");
    case TMC_SuComputeStep:    return("SU-ComputeStep");
    case TMC_SuFloating:       return("SU-Floating");
    case TMC_SuMotion:         return("SU-Motion");
    case TMC_SuSavePart:       return("SU-SavePart");
  }
  return("???");
}
//==============================================================================
/// Returns the number of timers.
//==============================================================================
inline unsigned TmcGetCount(){ return(TMC_COUNT); }

//==============================================================================
/// Creates timers to measure time intervals.
//==============================================================================
inline void TmcCreation(TimersCpu vtimer,bool active){
  for(unsigned c=0;c<TMC_COUNT;c++){
    StSphTimerCpu* t=vtimer+c;
    t->timer.Reset();
    t->active=active;
    t->time=0;
  }
}

//==============================================================================
/// Destroys the timers used to measure times.
//==============================================================================
inline void TmcDestruction(TimersCpu vtimer){ TmcCreation(vtimer,false); }

//==============================================================================
/// Marks start of timer.
//==============================================================================
inline void _TmcStart(TimersCpu vtimer,CsTypeTimerCPU ct){ if(vtimer[ct].active)vtimer[ct].timer.Start(); }

//==============================================================================
/// Marks end of timer and accumulates time.
//==============================================================================
inline void _TmcStop(TimersCpu vtimer,CsTypeTimerCPU ct){
  StSphTimerCpu* t=vtimer+unsigned(ct);
  if(t->active){
    t->timer.Stop();
    t->time+=t->timer.GetElapsedTimeF();
  }
}

//==============================================================================
/// Returns the time accumulated by the timer.
//==============================================================================
inline float TmcGetValue(const TimersCpu vtimer,CsTypeTimerCPU ct){ return(vtimer[ct].time); }

//==============================================================================
/// Indicates whether a timer is active or not.
//==============================================================================
inline bool TmcIsActive(const TimersCpu vtimer,CsTypeTimerCPU ct){ return(unsigned(ct)<TmcGetCount()&&vtimer[ct].active); }

#endif







