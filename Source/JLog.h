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

/// \file JLog.h \brief Declares the class \ref JLog.

#ifndef _JLog_
#define _JLog_

#include "JObject.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

//##############################################################################
//# JLog
//##############################################################################
/// \brief Allows creating a log file with info of the execution.

class JLog : protected JObject
{
public:
  typedef enum{ Out_ScrFile=3,Out_File=2,Out_Screen=1,Out_None=0 }TpMode_Out;
protected:
  std::ofstream *Pf;
  bool Ok;
public:
  JLog();
  ~JLog();
  void Reset();
  void Init(const std::string &fname);
  void Print(const std::string &tx,TpMode_Out mode=Out_ScrFile);
  bool IsOk()const{ return(Ok); }
};

#endif




