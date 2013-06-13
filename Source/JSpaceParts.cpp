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

/// \file JSpaceParts.cpp \brief Implements the class \ref JSpaceParts.

#include "JSpaceParts.h"
#include "JXml.h"


//##############################################################################
std::string JSpacePartBlock::GetNameXml()const{
  std::string ret;
  for(unsigned c=16;c<ClassName.length();c++)ret=ret+char(tolower(ClassName[c]));
  return(ret);
}
//==============================================================================
void JSpacePartBlock::ReadXml(JXml *sxml,TiXmlElement* ele){
  MkType=sxml->GetAttributeByte(ele,(Bound? "mkbound": "mkfluid")); 
  Begin=sxml->GetAttributeUnsigned(ele,"begin");
  Count=sxml->GetAttributeUnsigned(ele,"count");
}
//==============================================================================
TiXmlElement* JSpacePartBlock::WriteXml(JXml *sxml,TiXmlElement* ele)const{
  TiXmlElement item(GetNameXml().c_str());
  JXml::AddAttribute(&item,(Bound? "mkbound": "mkfluid"),MkType);
  JXml::AddAttribute(&item,"mk",Mk);
  JXml::AddAttribute(&item,"begin",Begin);
  JXml::AddAttribute(&item,"count",Count);
  return(ele->InsertEndChild(item)->ToElement());
}
//##############################################################################
void JSpacePartBlock_Moving::ReadXml(JXml *sxml,TiXmlElement* ele){
  JSpacePartBlock::ReadXml(sxml,ele);
  RefMotion=sxml->GetAttributeUnsigned(ele,"refmotion"); 
}
//==============================================================================
TiXmlElement* JSpacePartBlock_Moving::WriteXml(JXml *sxml,TiXmlElement* ele)const{
  ele=JSpacePartBlock::WriteXml(sxml,ele);
  JXml::AddAttribute(ele,"refmotion",RefMotion);
  return(ele);
}
//##############################################################################
void JSpacePartBlock_Floating::ReadXml(JXml *sxml,TiXmlElement* ele){
  JSpacePartBlock::ReadXml(sxml,ele);
  Massbody=sxml->ReadElementFloat(ele,"massbody","value");
  Center=sxml->ReadElementFloat3(ele,"center");
  Inertia=sxml->ReadElementFloat3(ele,"inertia");
  Velini=(sxml->GetFirstElement(ele,"velini",true)!=NULL? sxml->ReadElementFloat3(ele,"velini"): TFloat3(0));
  Omegaini=(sxml->GetFirstElement(ele,"omegaini",true)!=NULL? sxml->ReadElementFloat3(ele,"omegaini"): TFloat3(0));
}
//==============================================================================
TiXmlElement* JSpacePartBlock_Floating::WriteXml(JXml *sxml,TiXmlElement* ele)const{
  ele=JSpacePartBlock::WriteXml(sxml,ele);
  sxml->AddElementAttrib(ele,"massbody","value",Massbody);
  sxml->AddElementFloat3(ele,"center",Center);
  sxml->AddElementFloat3(ele,"inertia",Inertia);
  if(Velini!=TFloat3(0))sxml->AddElementFloat3(ele,"velini",Velini);
  if(Omegaini!=TFloat3(0))sxml->AddElementFloat3(ele,"omegaini",Omegaini);
  return(ele);
}

//##############################################################################
//##############################################################################
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSpaceParts::JSpaceParts(){
  ClassName="JSpaceParts";
  Reset();
}
//==============================================================================
/// Destructor.
//==============================================================================
JSpaceParts::~JSpaceParts(){
  Reset();
}
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSpaceParts::Reset(){
  for(unsigned c=0;c<Blocks.size();c++)delete Blocks[c];
  Blocks.clear();
  Begin=0;
  LastType=PT_Fixed;
  SetMkFirst(0,0);
}
//==============================================================================
/// Returns the number of particles of a given type.
//==============================================================================
unsigned JSpaceParts::Count(TpParticles type)const{
  unsigned n=0;
  for(unsigned c=0;c<Blocks.size();c++)if(Blocks[c]->Type==type)n+=Blocks[c]->GetCount();
  return(n);
}
//==============================================================================
/// Returns the total number of particles.
//==============================================================================
unsigned JSpaceParts::Count()const{
  unsigned n=0;
  for(unsigned c=0;c<Blocks.size();c++)n+=Blocks[c]->GetCount();
  return(n);
}
//==============================================================================
/// Returns the number of particles of the requested type.
//==============================================================================
unsigned JSpaceParts::CountBlocks(TpParticles type)const{
  unsigned n=0;
  for(unsigned c=0;c<Blocks.size();c++)if(Blocks[c]->Type==type)n++;
  return(n);
}

//==============================================================================
/// Returns the requested block.
//==============================================================================
const JSpacePartBlock& JSpaceParts::GetBlock(unsigned pos)const{
  if(pos>=CountBlocks())RunException("GetBlock","The requested particles block is missing.");
  return(*(Blocks[pos]));
}

//==============================================================================
/// Returns the block with the given mk.
//==============================================================================
JSpacePartBlock* JSpaceParts::GetByMkType(bool bound,byte mktype)const{
  JSpacePartBlock* bk=NULL;
  for(unsigned c=0;c<Blocks.size()&&!bk;c++)if((Blocks[c]->Type!=PT_Fluid)==bound&&Blocks[c]->GetMkType()==mktype)bk=Blocks[c];
  return(bk);
}

//==============================================================================
/// Checks and adds a new block of particles.
//==============================================================================
void JSpaceParts::Add(JSpacePartBlock* block){
  if(GetByMkType(block->Bound,block->GetMkType()))RunException("Add","Cannot add a block with a existing mk.");
  if(block->Type<LastType)RunException("Add","The block type is invalid after the last type added.");
  block->ConfigMk(block->Type==PT_Fluid? MkFluidFirst: MkBoundFirst);
  Blocks.push_back(block);
  Begin+=block->GetCount();
  LastType=block->Type;
}

//==============================================================================
/// Loads data in xml format from a file.
//==============================================================================
void JSpaceParts::LoadFileXml(const std::string &file,const std::string &path){
  JXml jxml;
  jxml.LoadFile(file);
  LoadXml(&jxml,path);
}

//==============================================================================
/// Stores data in xml format into a file.
//==============================================================================
void JSpaceParts::SaveFileXml(const std::string &file,const std::string &path,bool newfile)const{
  JXml jxml;
  if(!newfile)jxml.LoadFile(file);
  SaveXml(&jxml,path);
  jxml.SaveFile(file);
}

//==============================================================================
/// Loads initial conditions from the object xml.
//==============================================================================
void JSpaceParts::LoadXml(JXml *sxml,const std::string &place){
  Reset();
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)RunException("LoadXml",std::string("Cannot find the element \'")+place+"\'.");
  ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Stores initial conditions in the object xml.
//==============================================================================
void JSpaceParts::SaveXml(JXml *sxml,const std::string &place)const{
  WriteXml(sxml,sxml->GetNode(place,true)->ToElement());
}

//==============================================================================
/// Reads the list of initial conditions in xml format.
//==============================================================================
void JSpaceParts::ReadXml(JXml *sxml,TiXmlElement* lis){
  const char met[]="ReadXml";
  unsigned np=sxml->GetAttributeUnsigned(lis,"np");
  unsigned nb=sxml->GetAttributeUnsigned(lis,"nb");
  unsigned nbf=sxml->GetAttributeUnsigned(lis,"nbf");
  byte mkboundfirst=sxml->GetAttributeByte(lis,"mkboundfirst");
  byte mkfluidfirst=sxml->GetAttributeByte(lis,"mkfluidfirst");
  SetMkFirst(mkboundfirst,mkfluidfirst);
  TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    std::string cmd=ele->Value();
    if(cmd.length()&&cmd[0]!='_'){
      if(cmd=="fixed")         Add(new JSpacePartBlock_Fixed(sxml,ele));
      else if(cmd=="moving")   Add(new JSpacePartBlock_Moving(sxml,ele));
      else if(cmd=="floating") Add(new JSpacePartBlock_Floating(sxml,ele));
      else if(cmd=="fluid")    Add(new JSpacePartBlock_Fluid(sxml,ele));
      else sxml->ErrReadElement(ele,cmd,false);
    }
    ele=ele->NextSiblingElement();
  }
  Begin=Count();
  if(np!=Count()||nb!=np-Count(PT_Fluid)||nbf!=Count(PT_Fixed))RunException(met,"The amount of particles does not match the header.");
}

//==============================================================================
/// Writes the list in xml format.
//==============================================================================
void JSpaceParts::WriteXml(JXml *sxml,TiXmlElement* lis)const{
  unsigned np=Count();
  JXml::AddAttribute(lis,"np",np);
  JXml::AddAttribute(lis,"nb",np-Count(PT_Fluid));
  JXml::AddAttribute(lis,"nbf",Count(PT_Fixed));
  JXml::AddAttribute(lis,"mkboundfirst",GetMkBoundFirst());
  JXml::AddAttribute(lis,"mkfluidfirst",GetMkFluidFirst());
  for(unsigned c=0;c<Blocks.size();c++)Blocks[c]->WriteXml(sxml,lis);
}

//==============================================================================
/// Adjusts the value of Mk according to boundfirst and fluidfirst.
//==============================================================================
void JSpaceParts::SetMkFirst(byte boundfirst,byte fluidfirst){
  MkBoundFirst=boundfirst; MkFluidFirst=fluidfirst;
  for(unsigned c=0;c<Blocks.size();c++)Blocks[c]->ConfigMk(Blocks[c]->Type==PT_Fluid? MkFluidFirst: MkBoundFirst);
}




