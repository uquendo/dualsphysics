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

/// \file JLog.cpp \brief Implements the class \ref JLog.

#include "JLog.h"

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JLog::JLog(){
  ClassName="JLog";
  Pf=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JLog::~JLog(){
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JLog::Reset(){
  Ok=false;
  if(Pf){
    if(Pf->is_open())Pf->close();
    delete Pf; Pf=NULL;
  }
}

//==============================================================================
/// Initialisation file of logs.
//==============================================================================
void JLog::Init(const std::string &fname){
  Reset();
  Pf=new ofstream; 
  Pf->open(fname.c_str());
  if(Pf)Ok=true;
  else RunException("Init","Cannot open the file.",fname);
}
  
//==============================================================================
/// Visualises and/or stores information of the execution.
//==============================================================================
void JLog::Print(const std::string &tx,TpMode_Out mode){
  if(mode&Out_Screen)printf("%s\n",tx.c_str());
  if((mode&Out_File)&&Pf)(*Pf) << tx << endl;
}




