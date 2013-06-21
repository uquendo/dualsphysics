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

/// \file JCfgRun.h \brief Declares the class \ref JCfgRun.

#ifndef _JCfgRun_
#define _JCfgRun_

#pragma warning(disable : 4996) //Anula sprintf() deprecated.

#include "Types.h"
#include "Functions.h"
#include "JObject.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>


//##############################################################################
//# JCfgRun
//##############################################################################
/// \brief Defines the class responsible of collecting the execution parameters by command line.

class JCfgRun : protected JObject
{
public:
protected:
  bool SvDef;
  int DirsDef;

public:
  bool PrintInfo;
  bool Cpu;
  bool Gpu;
  int GpuId;
  
  int OmpThreads;
  TpOmpMode OmpMode;

  TpStep TStep;
  int VerletSteps;
  TpKernel TKernel;
  byte Kgc;   
  TpVisco TVisco; 
  float Visco;
  float TimeMax,TimePart;
  int ShepardSteps;
  int DBCSteps;
  bool SvDt,SvRes,SvTimers;
  bool Sv_Binx2,Sv_Csv,Sv_Ascii,Sv_Vtk,Sv_Flw;
  std::string CaseName,RunName,DirOut;
  std::string PartBeginDir;
  unsigned PartBegin,PartBeginFirst;
  float Incz;
  bool RhopOut;                   ///<Indicates whether \ref RhopOut density correction is actived or not.
  float RhopOutMin,RhopOutMax;    ///<Limits for \ref RhopOut density correction.
  TpCellMode  CellMode;           ///<Modes of cells division.
  std::string PtxasFile;          ///<File with ptxas information.

  JCfgRun();
  void Reset();
  void VisuInfo()const;
  void VisuConfig()const;
  void LoadArgv(int argc,char** argv);
  void LoadFile(std::string fname,int lv);
  void LoadOpts(std::string *optlis,int optn,int lv,std::string file);
  void ErrorParm(const std::string &opt,int optc,int lv,const std::string &file)const;
};

#endif







