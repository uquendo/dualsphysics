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

/// \file JSpaceCtes.h \brief Declares the class \ref JSpaceCtes.

#ifndef _JSpaceCtes_
#define _JSpaceCtes_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"

class JXml;
class TiXmlElement;

//#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

//##############################################################################
//# JSpaceCtes
//##############################################################################
/// \brief Manages the info of constants from the input XML file.

class JSpaceCtes : protected JObject 
{
private:
  int LatticeBound;       ///<Lattice to create boundary particles on its nodes.
  int LatticeFluid;       ///<Lattice to create fluid particles on its nodes.
  tfloat3 Gravity;        ///<Gravity acceleration.
  float CFLnumber;        ///<CFL number (0.1-0.5).
  bool HSwlAuto;          ///<Activates the automatic computation of H_Swl.
  float HSwl;             ///<Maximum height of the volume of fluid.
  float CoefSound;        ///<Coefficient of speed of sound.
  float Coefficient;      ///<Coefficient to calculate the smoothing length H.
  float Gamma;            ///<Politropic constant. (1-7).
  float Rhop0;            ///<Density of reference.
  float Eps;              ///<Epsilon constant for XSPH variant.
  //-Computed values:
  float Dp;               ///<Inter-particle distance.
  float H;                ///<Smoothing length.
  float B;                ///<Constant that sets a limit for the maximum change in density.
  float MassBound;        ///<Mass of a boundary particle.
  float MassFluid;        ///<Mass of a fluid particle.

  void ReadXmlDef(JXml *sxml,TiXmlElement* ele);
  void WriteXmlDef(JXml *sxml,TiXmlElement* ele)const;
  void ReadXmlRun(JXml *sxml,TiXmlElement* ele);
  void WriteXmlRun(JXml *sxml,TiXmlElement* ele)const;
public:
  
  JSpaceCtes();
  void Reset();
  void ResetCalc();
  void LoadDefault();
  void LoadXmlDef(JXml *sxml,const std::string &place);
  void SaveXmlDef(JXml *sxml,const std::string &place)const;
  void LoadXmlRun(JXml *sxml,const std::string &place);
  void SaveXmlRun(JXml *sxml,const std::string &place)const;

  int GetLatticeBound()const{ return(LatticeBound); }
  int GetLatticeFluid()const{ return(LatticeFluid); }
  tfloat3 GetGravity()const{ return(Gravity); }
  float GetCFLnumber()const{ return(CFLnumber); }
  bool GetHSwlAuto()const{ return(HSwlAuto); }
  float GetHSwl()const{ return(HSwl); }
  float GetCoefSound()const{ return(CoefSound); }
  float GetCoefficient()const{ return(Coefficient); }
  float GetGamma()const{ return(Gamma); }
  float GetRhop0()const{ return(Rhop0); }
  float GetEps()const{ return(Eps); }

  void SetLatticeBound(bool simple){ LatticeBound=(simple? 1: 2); }
  void SetLatticeFluid(bool simple){ LatticeFluid=(simple? 1: 2); }
  void SetGravity(const tfloat3& g){ Gravity=g; }
  void SetCFLnumber(float v){ 
    if(v<0.1f||v>0.5f)RunException("SetCFLnumber","Value out of allowed range (0.1-0.5).");
    CFLnumber=v;
  }
  void SetHSwlAuto(bool on){ HSwlAuto=on; }
  void SetHSwl(float v){ HSwl=v; }
  void SetCoefSound(float v){ CoefSound=v; }
  void SetCoefficient(float v){ Coefficient=v; }
  void SetGamma(float v){ Gamma=v; }
  void SetRhop0(float v){ Rhop0=v; }
  void SetEps(float v){ Eps=v; }

  float GetDp()const{ return(Dp); }
  float GetH()const{ return(H); }
  float GetB()const{ return(B); }
  float GetMassBound()const{ return(MassBound); }
  float GetMassFluid()const{ return(MassFluid); }

  void SetDp(float v){ Dp=v; }
  void SetH(float v){ H=v; }
  void SetB(float v){ B=v; }
  void SetMassBound(float v){ MassBound=v; }
  void SetMassFluid(float v){ MassFluid=v; }
};

#endif




