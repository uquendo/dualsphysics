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

/// \file JSphGpu.h \brief Declares the class \ref JSphGpu

#ifndef _JSphGpu_
#define _JSphGpu_

#include "TypesDef.h"
#include "JSphTimersCpu.h"
#include "CudaSphApi.h"
#include "JSph.h"
#include "JDivideCpu.h"

#include <string>

//##############################################################################
//# JSphGpu
//##############################################################################
/// \brief Defines the attributes and functions used only in GPU simulations.

class JSphGpu : public JSph
{
protected:
  std::string PtxasFile;      ///<File with information of number of registers to optimise the execution.
  std::string BlockSizes;     ///<Stores configuration of BlockSizes.
  StDeviceContext Dc;         ///<Structure with variables to be used on the CUDA files.
  StDeviceCte CteVars;        ///<Structure with variables to be stored in the constant memory of GPU.
  int dev;                    ///<GpuId.

  void ConfigDevice();

  void AllocMemory();
  void Reset();
  void DivideConfig(tfloat3 &posmin,tfloat3 &difmax);
  void GetVarsDevice();
  void SaveData();
  void ShowTimers();
  static unsigned OptimizeBlockSize(unsigned compute,unsigned nreg);
  unsigned BlockSizeConfig(const string& opname,unsigned compute,unsigned regs);

  void RunMotion(float stepdt);
  void RunFloating(float stepdt);
  void InitVars();

public:
  JSphGpu();
  ~JSphGpu();
  void Run(JCfgRun *cfg,JLog *log);

  unsigned GetMemoryCpu()const{ return(JSph::GetMemoryCpu()+Dc.memcpu); }
  unsigned GetMemoryGpu()const{ return(Dc.memgpu+Dc.memgpu2); }

  unsigned TimerGetCount()const{ return(TmgGetCount()); }
  bool TimerIsActive(unsigned ct)const{ return(TmgIsActive(Dc.timers,(CsTypeTimerGPU)ct)); }
  float TimerGetValue(unsigned ct)const{ return(TmgGetValue(Dc.timers,(CsTypeTimerGPU)ct)); }
  string TimerGetName(unsigned ct)const{ return(TmgGetName((CsTypeTimerGPU)ct)); }
};

#endif






