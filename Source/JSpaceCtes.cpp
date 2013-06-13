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

/// \file JSpaceCtes.cpp \brief Implements the class \ref JSpaceCtes.

#include "JSpaceCtes.h"
#include "JXml.h"

//==============================================================================
/// Constructor.
//==============================================================================
JSpaceCtes::JSpaceCtes(){
  ClassName="JSpaceCtes";
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSpaceCtes::Reset(){
  SetLatticeBound(true);
  SetLatticeFluid(true);
  Gravity=TFloat3(0);
  CFLnumber=0; HSwlAuto=true; HSwl=0;
  CoefSound=0; Coefficient=0; Gamma=0; Rhop0=0; Eps=0;
  ResetCalc();
}

//==============================================================================
/// Resets the computed variables.
//==============================================================================
void JSpaceCtes::ResetCalc(){
  Dp=0; H=0; B=0;
  MassBound=0; MassFluid=0;
}

//==============================================================================
/// Loads values by default.
//==============================================================================
void JSpaceCtes::LoadDefault(){
  Reset();
  SetLatticeBound(true);
  SetLatticeFluid(true);
  SetGravity(TFloat3(0,0,-9.81f));
  SetCFLnumber(0.2f);
  SetHSwlAuto(true);
  SetHSwl(0);
  SetCoefSound(10);
  SetCoefficient(0.866025f);
  SetGamma(7);
  SetRhop0(1000);
  SetEps(0.5f);
}

//==============================================================================
/// Reads constants for definition of the case of xml node.
//==============================================================================
void JSpaceCtes::ReadXmlDef(JXml *sxml,TiXmlElement* node){
  TiXmlElement* lattice=sxml->GetFirstElement(node,"lattice");
  SetLatticeBound(sxml->GetAttributeInt(lattice,"bound")==1);
  SetLatticeFluid(sxml->GetAttributeInt(lattice,"fluid")==1);
  SetGravity(sxml->ReadElementFloat3(node,"gravity"));
  SetCFLnumber(sxml->ReadElementFloat(node,"cflnumber","value"));
  TiXmlElement* hswl=sxml->GetFirstElement(node,"hswl");
  SetHSwl(sxml->GetAttributeFloat(hswl,"value"));
  SetHSwlAuto(sxml->GetAttributeBool(hswl,"auto"));
  SetCoefSound(sxml->ReadElementFloat(node,"coefsound","value"));
  SetCoefficient(sxml->ReadElementFloat(node,"coefficient","value"));
  SetGamma(sxml->ReadElementFloat(node,"gamma","value"));
  SetRhop0(sxml->ReadElementFloat(node,"rhop0","value"));
  SetEps(sxml->ReadElementFloat(node,"eps","value"));
}

//==============================================================================
/// Writes constants for definition of the case of xml node.
//==============================================================================
void JSpaceCtes::WriteXmlDef(JXml *sxml,TiXmlElement* node)const{
  TiXmlElement lattice("lattice");
  JXml::AddAttribute(&lattice,"bound",GetLatticeBound());
  JXml::AddAttribute(&lattice,"fluid",GetLatticeFluid());
  node->InsertEndChild(lattice);
  JXml::AddElementFloat3(node,"gravity",GetGravity());
  JXml::AddElementAttrib(node,"cflnumber","value",GetCFLnumber());
  TiXmlElement hswl("hswl");
  JXml::AddAttribute(&hswl,"value",GetHSwl()); 
  JXml::AddAttribute(&hswl,"auto",GetHSwlAuto());
  node->InsertEndChild(hswl);
  JXml::AddElementAttrib(node,"coefsound","value",GetCoefSound());
  JXml::AddElementAttrib(node,"coefficient","value",GetCoefficient());
  JXml::AddElementAttrib(node,"gamma","value",GetGamma());
  JXml::AddElementAttrib(node,"rhop0","value",GetRhop0());
  JXml::AddElementAttrib(node,"eps","value",GetEps());
}

//==============================================================================
/// Reads constants for execution of the case of xml node.
//==============================================================================
void JSpaceCtes::ReadXmlRun(JXml *sxml,TiXmlElement* node){
  SetGravity(sxml->ReadElementFloat3(node,"gravity"));
  SetCFLnumber(sxml->ReadElementFloat(node,"cflnumber","value"));
  SetGamma(sxml->ReadElementFloat(node,"gamma","value"));
  SetRhop0(sxml->ReadElementFloat(node,"rhop0","value"));
  SetEps(sxml->ReadElementFloat(node,"eps","value"));
  SetDp(sxml->ReadElementFloat(node,"dp","value"));
  SetH(sxml->ReadElementFloat(node,"h","value"));
  SetB(sxml->ReadElementFloat(node,"b","value"));
  SetMassBound(sxml->ReadElementFloat(node,"massbound","value"));
  SetMassFluid(sxml->ReadElementFloat(node,"massfluid","value"));
}

//==============================================================================
/// Writes constants for execution of the case of xml node.
//==============================================================================
void JSpaceCtes::WriteXmlRun(JXml *sxml,TiXmlElement* node)const{
  JXml::AddElementFloat3(node,"gravity",GetGravity());
  JXml::AddElementAttrib(node,"cflnumber","value",GetCFLnumber());
  JXml::AddElementAttrib(node,"gamma","value",GetGamma());
  JXml::AddElementAttrib(node,"rhop0","value",GetRhop0());
  JXml::AddElementAttrib(node,"eps","value",GetEps());
  JXml::AddElementAttrib(node,"dp","value",GetDp());
  JXml::AddElementAttrib(node,"h","value",GetH(),"%.7E");
  JXml::AddElementAttrib(node,"b","value",GetB(),"%.7E");
  JXml::AddElementAttrib(node,"massbound","value",GetMassBound(),"%.7E");
  JXml::AddElementAttrib(node,"massfluid","value",GetMassFluid(),"%.7E");
}

//==============================================================================
/// Loads constants for execution of the case of xml node.
//==============================================================================
void JSpaceCtes::LoadXmlDef(JXml *sxml,const std::string &place){
  Reset();
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)RunException("LoadXmlDef",std::string("The item is not found \'")+place+"\'.");
  ReadXmlDef(sxml,node->ToElement());
}

//==============================================================================
/// Stores constants for execution of the case of xml node.
//==============================================================================
void JSpaceCtes::SaveXmlDef(JXml *sxml,const std::string &place)const{
  WriteXmlDef(sxml,sxml->GetNode(place,true)->ToElement());
}

//==============================================================================
/// Loads constants for execution of the case of xml node.
//==============================================================================
void JSpaceCtes::LoadXmlRun(JXml *sxml,const std::string &place){
  Reset();
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)RunException("LoadXmlRun",std::string("The item is not found \'")+place+"\'.");
  ReadXmlRun(sxml,node->ToElement());
}

//==============================================================================
/// Stores constants for execution of the case of xml node.
//==============================================================================
void JSpaceCtes::SaveXmlRun(JXml *sxml,const std::string &place)const{
  WriteXmlRun(sxml,sxml->GetNode(place,true)->ToElement());
}




