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

/// \file JSpaceParts.h \brief Declares the class \ref JSpaceParts.

#ifndef _JSpaceParts_
#define _JSpaceParts_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"

class JXml;
class TiXmlElement;

///<Types of particles.
typedef enum{ PT_Fixed=1,PT_Moving=2,PT_Floating=4,PT_Fluid=5 }TpParticles; 

//==============================================================================
//##############################################################################
//==============================================================================
class JSpacePartBlock : public JObject
{
private:
  byte Mk;                   ///<Absolute label.
  byte MkType;               ///<Label of block fluid or bound.
  unsigned Begin;            ///<Id of the first particle of the block.
  unsigned Count;            ///<Number of particles.
public:
  const TpParticles Type;    ///<Type of particle.
  const bool Bound;          ///<Indicates whether a particle is boundary or not.

  JSpacePartBlock(TpParticles type,const char* name,bool bound,byte mktype=0,unsigned begin=0,unsigned count=0):Type(type),Bound(bound),MkType(mktype),Begin(begin),Count(count){ 
    ClassName=std::string("JSpacePartBlock_")+name;
  } 
  virtual ~JSpacePartBlock(){}
  void ConfigMk(byte mkfirst){ Mk=MkType+mkfirst; }
  std::string GetNameXml()const;
  unsigned GetBegin()const{ return(Begin); }
  unsigned GetCount()const{ return(Count); }
  unsigned GetMkType()const{ return(MkType); }
  unsigned GetMk()const{ return(Mk); }
  virtual void ReadXml(JXml *sxml,TiXmlElement* ele);
  virtual TiXmlElement* WriteXml(JXml *sxml,TiXmlElement* ele)const;
};

//##############################################################################
class JSpacePartBlock_Fixed : public JSpacePartBlock
{
public:
  JSpacePartBlock_Fixed(byte mktype,unsigned begin,unsigned count):JSpacePartBlock(PT_Fixed,"Fixed",true,mktype,begin,count){}
  JSpacePartBlock_Fixed(JXml *sxml,TiXmlElement* ele):JSpacePartBlock(PT_Fixed,"Fixed",true){ ReadXml(sxml,ele); }
};  
//##############################################################################
class JSpacePartBlock_Moving : public JSpacePartBlock
{
private:
  unsigned RefMotion;
public:
  JSpacePartBlock_Moving(byte mktype,unsigned begin,unsigned count,unsigned refmotion):JSpacePartBlock(PT_Moving,"Moving",true,mktype,begin,count),RefMotion(refmotion){}
  JSpacePartBlock_Moving(JXml *sxml,TiXmlElement* ele):JSpacePartBlock(PT_Moving,"Moving",true){ ReadXml(sxml,ele); }
  unsigned GetRefMotion()const{ return(RefMotion); }
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  TiXmlElement* WriteXml(JXml *sxml,TiXmlElement* ele)const;
};  
//##############################################################################
class JSpacePartBlock_Floating : public JSpacePartBlock
{
private:
  float Massbody;
  tfloat3 Center;
  tfloat3 Inertia;
  tfloat3 Velini;
  tfloat3 Omegaini;
public:
  JSpacePartBlock_Floating(byte mktype,unsigned begin,unsigned count,float massbody,const tfloat3& center,const tfloat3& inertia,const tfloat3& velini,const tfloat3& omegaini):JSpacePartBlock(PT_Floating,"Floating",true,mktype,begin,count),Massbody(massbody),Center(center),Inertia(inertia),Velini(velini),Omegaini(omegaini){}
  JSpacePartBlock_Floating(JXml *sxml,TiXmlElement* ele):JSpacePartBlock(PT_Floating,"Floating",true){ ReadXml(sxml,ele); }
  float GetMassbody()const{ return(Massbody); }
  tfloat3 GetCenter()const{ return(Center); }
  tfloat3 GetInertia()const{ return(Inertia); }
  tfloat3 GetVelini()const{ return(Velini); }
  tfloat3 GetOmegaini()const{ return(Omegaini); }
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  TiXmlElement* WriteXml(JXml *sxml,TiXmlElement* ele)const;
};  
//##############################################################################
class JSpacePartBlock_Fluid : public JSpacePartBlock
{
public:
  JSpacePartBlock_Fluid(byte mktype,unsigned begin,unsigned count):JSpacePartBlock(PT_Fluid,"Fluid",false,mktype,begin,count){}
  JSpacePartBlock_Fluid(JXml *sxml,TiXmlElement* ele):JSpacePartBlock(PT_Fluid,"Fluid",false){ ReadXml(sxml,ele); }
};  

//##############################################################################
//# JSpaceParts
//##############################################################################
/// \brief Manages the info of particles from the input XML file.

class JSpaceParts  : protected JObject
{
private:
  std::vector<JSpacePartBlock*> Blocks;
  unsigned Begin;
  TpParticles LastType;
  byte MkBoundFirst,MkFluidFirst;
  
  unsigned GetBegin()const{ return(Begin); }
  JSpacePartBlock* GetByMkType(bool bound,byte mktype)const;
  void Add(JSpacePartBlock* block);
  void ReadXml(JXml *sxml,TiXmlElement* lis);
  void WriteXml(JXml *sxml,TiXmlElement* lis)const;

public:
  JSpaceParts();
  ~JSpaceParts();
  void Reset();
  unsigned Count(TpParticles type)const;
  unsigned Count()const;
  unsigned CountBlocks()const{ return(unsigned(Blocks.size())); }
  unsigned CountBlocks(TpParticles type)const;
  const JSpacePartBlock& GetBlock(unsigned pos)const;

  void LoadFileXml(const std::string &file,const std::string &path);
  void SaveFileXml(const std::string &file,const std::string &path,bool newfile=true)const;
  void LoadXml(JXml *sxml,const std::string &place);
  void SaveXml(JXml *sxml,const std::string &place)const;

  void SetMkFirst(byte boundfirst,byte fluidfirst);
  byte GetMkBoundFirst()const{ return(MkBoundFirst); }
  byte GetMkFluidFirst()const{ return(MkFluidFirst); }

  void AddFixed(byte mktype,unsigned count){ Add(new JSpacePartBlock_Fixed(mktype,GetBegin(),count)); }
  void AddMoving(byte mktype,unsigned count,unsigned refmotion){ Add(new JSpacePartBlock_Moving(mktype,GetBegin(),count,refmotion)); }
  void AddFloating(byte mktype,unsigned count,float massbody,const tfloat3& center,const tfloat3& inertia,const tfloat3& velini,const tfloat3& omegaini){ Add(new JSpacePartBlock_Floating(mktype,GetBegin(),count,massbody,center,inertia,velini,omegaini)); }
  void AddFluid(byte mktype,unsigned count){ Add(new JSpacePartBlock_Fluid(mktype,GetBegin(),count)); }
};


#endif




