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

/// \file JException.cpp \brief Implements the class \ref JException.

#include "JException.h"
#include <cstdio>

//==============================================================================
/// Constructor of objects.
/// \param classname Name of the class that throws the exception.
/// \param method Name of the method that throws the exception.
/// \param text Text of the exception.
/// \param file Name of the file related to the exception.
//==============================================================================
JException::JException(const std::string &classname,const std::string &method,const std::string &text,const std::string &file){
  ExName="JException";
  ClassName=classname;
  Method=method;
  Text=text;
  File=file;
}

//==============================================================================
/// Returns the complete text message with the information of the exception. 
//==============================================================================
std::string JException::ToStr()const{
  std::string tx;
  char cad[512];
  sprintf(cad,"Exception (%s::%s)\n",ClassName.c_str(),Method.c_str());  tx=cad;
  if(!Text.empty()){ sprintf(cad,"Text: %s\n",Text.c_str()); tx=tx+cad; }
  if(!File.empty()){ sprintf(cad,"File: %s\n",File.c_str()); tx=tx+cad; }
  return(tx);
}

//==============================================================================
/// Visualises the exception message in console.
//==============================================================================
void JException::Print()const{
  printf("\n*** %s\n",ToStr().c_str());
}




