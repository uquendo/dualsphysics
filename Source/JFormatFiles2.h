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

/// \file JFormatFiles2.h \brief Declares the class \ref JFormatFiles2.

#ifndef _JFormatFiles2_
#define _JFormatFiles2_

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

#include <string>
#include <cstring>
#include "TypesDef.h"

//##############################################################################
//# JFormatFiles2
//##############################################################################
/// \brief Provides functions to store particle data in formats VTK, CSV, ASCII.

class JFormatFiles2
{
public:

  //==============================================================================
  /// Throws a simple exception.
  //==============================================================================
  static void RunException(std::string method,std::string msg);
  
  //==============================================================================
  /// Throws an exception related to a file.
  //==============================================================================  
  static void RunException(std::string method,std::string msg,std::string file);

  //==============================================================================
  /// Stores data in CSV format.
  //==============================================================================
  static void ParticlesToCsv(std::string fname,unsigned np,unsigned nfixed,unsigned nmoving,unsigned nfloat,unsigned nfluid,unsigned nfluidout,float timestep,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press,const float *mass,const unsigned *id,const byte *type,const byte *mk,const tfloat3 *ace,const tfloat3 *vor);
  
  //==============================================================================
  /// Stores data in ASCII format.
  //==============================================================================  
  static void ParticlesToAscii(std::string fname,unsigned np,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press,const float *mass,const unsigned *id,const byte *type,const byte *mk,const tfloat3 *ace,const tfloat3 *vor);
  
  //==============================================================================
  /// Stores data in CSV format (splits positive and negative part of Ace).
  //==============================================================================  
  static void ParticlesToCsv2(std::string fname,unsigned np,unsigned nfixed,unsigned nmoving,unsigned nfloat,unsigned nfluid,unsigned nfluidout,float timestep,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press,const float *mass,const unsigned *id,const byte *type,const byte *mk,const tfloat3 *acepos,const tfloat3 *aceneg,const tfloat3 *vor);

  //==============================================================================
  /// Stores data in ASCII format (splits positive and negative part of Ace).
  //==============================================================================
  static void ParticlesToAscii2(std::string fname,unsigned np,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press,const float *mass,const unsigned *id,const byte *type,const byte *mk,const tfloat3 *acepos,const tfloat3 *aceneg,const tfloat3 *vor);

  //==============================================================================
  /// Stores data in VTK format.
  //============================================================================== 
  static void ParticlesToVtk(std::string fname,unsigned np,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press,const float *mass,const unsigned *id,const byte *type,const byte *mk,const tfloat3 *ace,const tfloat3 *vor,int domain=0);

  //==============================================================================
  /// Stores data in VTK format for variables of type float.
  //==============================================================================
  static void ParticlesToVtkFloat(std::string fname,unsigned np,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press,const float *mass,const float *id,const float *type,const float *mk,const tfloat3 *ace,const tfloat3 *vor);

  //==============================================================================
  /// Stores information of points in  VTK format.
  //==============================================================================
  static void SaveVtkPointsVar(std::string fname,unsigned np,const tfloat3 *pos,const std::string &varname,const float *var);

  //==============================================================================
  /// Stores information of points in CSV file for variables of type float.
  //==============================================================================
  static void SaveCsvPointsVar(const std::string &fname,int part,float timestep,unsigned np,const tfloat3* pos,const float* data,bool first=false);

  //==============================================================================
  /// Stores information of points in CSV file for variables of type float3.
  //==============================================================================
  static void SaveCsvPointsVar3(const std::string &fname,int part,float timestep,unsigned np,const tfloat3* pos,const tfloat3* data,bool first=false);

  //==============================================================================
  /// Stores information of points in ASCII file for variables of type float.
  //==============================================================================  
  static void SaveAscPointsVar(const std::string &fname,float timestep,unsigned np,const tfloat3* pos,const float* data,bool first=false);

  //==============================================================================
  /// Stores information of points in ASCII file for variables of type float3.
  //============================================================================== 
  static void SaveAscPointsVar3(const std::string &fname,float timestep,unsigned np,const tfloat3* pos,const tfloat3* data,bool first=false);

  //==============================================================================
  /// Stores time and position for predefined motion.
  //==============================================================================
  static void SaveMotionPredef(const std::string &fname,unsigned np,const float *time,const tfloat3 *pos);

  //==============================================================================
  /// Stores time and position for predefined motion.
  //==============================================================================
  static void SaveMotionPredef(const std::string &fname,unsigned np,const float *time,const float *pos);

  //==============================================================================
  /// Generates a VTK file with map cells.
  //==============================================================================
  static void SaveVtkCells(const std::string &fname,const tfloat3 &posmin,const tuint3 &cells,float scell);

};


#endif




