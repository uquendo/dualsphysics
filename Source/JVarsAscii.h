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

/// \file JVarsAscii.h \brief Declares the class \ref JVarsAscii.

#ifndef _JVarsAscii_
#define _JVarsAscii_

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

#include "JObject.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>


//##############################################################################
//# JVarsAscii
//##############################################################################
/// \brief Reads variables from a text file in ASCII format.

class JVarsAscii : protected JObject
{
protected:
  static const unsigned MAXVARS=200;
  unsigned Count;
  std::string *Names,*Values;

  unsigned AddValue(const std::string &name);

public:
  JVarsAscii();
  ~JVarsAscii();
  void Reset();
  void LoadFile(const std::string &file);
  void SaveFile(const std::string &file)const;
  unsigned GetCount()const{ return(Count); }
  std::string GetName(unsigned pos)const{ return(pos<Count? Names[pos]: ""); }
  std::string GetValue(unsigned pos)const{ return(pos<Count? Values[pos]: ""); }
  int FindName(const std::string &name)const;

  bool Exists(const std::string &name,unsigned value=0)const{ return(ValuesCount(name)!=0); }
  unsigned ValuesCount(const std::string &name)const{ return(ValuesCount(FindName(name))); }
  unsigned ValuesCount(unsigned pos)const;

  std::string ValueString(const std::string &name,unsigned value)const;
  int ValueInt(const std::string &name,unsigned value=0)const;
  double ValueDouble(const std::string &name,unsigned value=0)const;
  float ValueFloat(const std::string &name,unsigned value=0)const{ return(float(ValueDouble(name,value))); }
  bool ValueBool(const std::string &name,unsigned value=0)const{ return(ValueInt(name,value)!=0); }

  unsigned AddValueString(const std::string &name,const std::string &v1="",const std::string &v2="",const std::string &v3="");
  unsigned AddValueInt(const std::string &name,int v1);
  unsigned AddValueInt(const std::string &name,int v1,int v2);
  unsigned AddValueInt(const std::string &name,int v1,int v2,int v3);
  unsigned AddValueDouble(const std::string &name,double v1,const char fmt[]="%g");
  unsigned AddValueDouble(const std::string &name,double v1,double v2,const char fmt[]="%g");
  unsigned AddValueDouble(const std::string &name,double v1,double v2,double v3,const char fmt[]="%g");

  static std::string StrTrim(const std::string &cad);
  static std::string StrUpper(const std::string &cad);
};

#endif




