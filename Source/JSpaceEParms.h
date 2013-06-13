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

/// \file JSpaceEParms.h \brief Declares the class \ref JSpaceEParms.

#ifndef _JSpaceEParms_
#define _JSpaceEParms_


#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

#include <string>
#include <vector>
#include "JObject.h"

class JXml;
class TiXmlElement;

//##############################################################################
//# JSpaceEParms
//##############################################################################
/// \brief Manages the info of execution parameters from the input XML file.

class JSpaceEParms : protected JObject
{
public:
  typedef struct{
    std::string key;
    std::string value;
    std::string comment;
  }JSpaceEParmsItem;
private:
  typedef std::vector<JSpaceEParmsItem> VecList;
  typedef std::vector<JSpaceEParmsItem>::iterator VecListIt;

  VecList List;

  JSpaceEParmsItem* GetItemPointer(const std::string &key);
  std::string GetValueNum(const std::string &key,int num);
  void ReadXml(JXml *sxml,TiXmlElement* lis);
  void WriteXml(JXml *sxml,TiXmlElement* lis)const;
public:
  JSpaceEParms();
  ~JSpaceEParms();
  void Reset();
  void Add(const std::string &key,const std::string &value,const std::string &comment);
  void SetValue(const std::string &key,const std::string &value);
  void SetComment(const std::string &key,const std::string &comment);
  bool Exists(const std::string &key){ return(GetItemPointer(key)!=NULL); }

  int GetValueNumInt(const std::string &key,int num,bool optional=false,int valdef=0);
  double GetValueNumDouble(const std::string &key,int num,bool optional=false,double valdef=0);
  float GetValueNumFloat(const std::string &key,int num,bool optional=false,float valdef=0){ return(float(GetValueNumDouble(key,num,optional,valdef))); }
  
  int GetValueInt(const std::string &key,bool optional=false,int valdef=0){ return(GetValueNumInt(key,0,optional,valdef)); }
  double GetValueDouble(const std::string &key,bool optional=false,double valdef=0){ return(GetValueNumDouble(key,0,optional,valdef)); }
  float GetValueFloat(const std::string &key,bool optional=false,float valdef=0){ return(GetValueNumFloat(key,0,optional,valdef)); }

  unsigned Count()const{ return(unsigned(List.size())); }
  std::string ToString(unsigned pos)const;
  JSpaceEParmsItem GetParm(unsigned pos)const;
  void LoadFileXml(const std::string &file,const std::string &path);
  void SaveFileXml(const std::string &file,const std::string &path,bool newfile=true)const;
  void LoadXml(JXml *sxml,const std::string &place);
  void SaveXml(JXml *sxml,const std::string &place)const;
};

#endif




