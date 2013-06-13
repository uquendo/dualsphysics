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

/// \file JDivideCpu.h \brief Declares the class \ref JDivideCpu.

#ifndef _JDivideCpu_
#define _JDivideCpu_


#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

#include "TypesDef.h"
#include "JObject.h"
#include "JLog.h"
#include <cmath>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
using namespace std;

//##############################################################################
//# JDivideCpu
//##############################################################################
/// \brief Defines the class responsible of computing the Neighbour List

class JDivideCpu : protected JObject
{
public:

/// Structure with the information used to compute the Neighbour List.
  typedef struct {
    float scell;             ///<Cells size. 
    float border;            ///<Limit applied to posmin and posmax.
    float incz;              ///<Allowable increase in Z+.
    tfloat3 posmin,posmax;   ///<Minimum and maximum position of the particles.
    tfloat3 difmax;          ///<Difference between posmax and posmin.
    unsigned ncellx;         ///<Number of cells in axis X.
    unsigned ncelly;         ///<Number of cells in axis Y.
    unsigned ncellz;         ///<Number of cells in axis Z.
    unsigned ncells;         ///<Total number of cells. 
    float poszmax_incz;      ///<New value of posmax in Z after applying incz.
    float difzmax_incz;      ///<New value of difzmax in Z after applying incz.
    unsigned ncellz_incz;    ///<New number of cells in Z after applying incz. 
    unsigned ncells_incz;    ///<New number of total cells after applying incz.     
  }StDivideInfo;

/// Structure with the information of cells.
  typedef struct {
    unsigned cells_void;
    unsigned min_incell;
    unsigned min_count;
    unsigned max_incell;
    unsigned max_count;
    float media;
  }StCellInfo;

/// Structure with the information of excluded particles.
  typedef struct{ 
    float timeout;
    unsigned id;
    tfloat3 pos;
    tfloat3 vel;
    float rhop;
    unsigned p;
  }StParticleOut;


protected:
  unsigned *CellPart;         ///<Cell of each particle. [Np]
  unsigned *Parts;            ///<All particles ordered according to the cell they belong. [Np]
  unsigned *PartsInCell;      ///<Number of particles per cell at the last Neighbour List. [NctTot]
  unsigned *BeginBound;       ///<Position of the first boundary particle of each cell. [NctTot+1]
  unsigned *BeginFluid;       ///<Position of the first fluid particle of each cell. [NctTot+1]

  //-Variables to reordered particles.
  byte *VSort;                ///<Memory to reorder particles. [sizeof(tfloat3)*Np]
  int *VSortInt;              ///<To order arrays of type int (points to VSort).
  float *VSortFloat;          ///<To order arrays of type float (points to VSort).
  tfloat3 *VSortFloat3;       ///<To order arrays of type tfloat3 (points to VSort).
  tsymatrix3f *VSortFloat6;   ///<To order arrays of type tsymatrix3f (points to VSort).
  
  //-Variables to manage particles out.
  unsigned OutSize;           ///<Number of particles for which there is memory space in \ref Out.
  unsigned OutCount;          ///<Number of excluded particles stored in \ref Out, that can be larger than \ref NfluidOut.
  StParticleOut* Out;         ///<List of excluded particles by the exclusion order.
  
  int *PartsOut;              ///<Excluded particles at the last Neighbour List step. [Np-Nbound]

  tfloat3 PosMin,DifMax;
  float IncZ;

  unsigned Np;                ///<Number of total particles. 
  unsigned Nfixed;            ///<Number of fixed boundary particles. 
  unsigned Nmoving;           ///<Number of moving boundary particles. 
  unsigned Nfloat;            ///<Number of floating boundary particles. 
  unsigned Nbound;            ///<Number of boundary particles ( \ref Nfixed + \ref Nmoving + \ref Nfloat ). 
  unsigned Npb;              ///<Number of particles of the boundary block ( \ref Nbound - \ref Nfloat ) or ( \ref Nfixed + \ref Nmoving).
  unsigned Nfluid;            ///<Number of fluid particles. 
  unsigned NfluidOut;         ///<Number of fluid particles excluded from the total \ref Np.
  unsigned NOutLast;          ///<Number of fluid particles excluded at the last Neighbour List.


  unsigned Ndivb;             ///<Number of times that \ref DivBoundary was called.
  unsigned Ndivf;             ///<Number of times that \ref DivFluid was called.

  tint3 Ncells;               ///<Number of cells for each direction.
  unsigned Nct;               ///<Number of total cells.
  unsigned Nsheet;            ///<Number of cells in X direction * Y direction.
  unsigned NctTot,NctTotMax;
  float SizeCell;             ///<Cells size.
  float OvSizeCell;           ///<Inverse value of Cells size.
  unsigned BoundMaxCellZ,FluidMaxCellZ,LastMaxCellZ;

  bool CellInfoOk;
  StCellInfo CellInfo;       ///<Collects information about cells.
  string DirOut;             ///<Name of the output directory.

  bool RhopOut;                ///<Indicates whether \ref RhopOut density correction is actived or not.
  float RhopOutMin,RhopOutMax; ///<Limits for \ref RhopOut density correction.
  unsigned RhopOutCount;       ///<Number of excluded particles by \ref RhopOut.

  void CalcCells();
  void VisuLimits()const;
  void ClearMemory();
  void CheckLimits(unsigned pini,unsigned n,const tfloat3 *pos,tfloat3 &pmin,tfloat3 &pmax)const;

  void OutResize(unsigned size);


public:
  JLog *Log;
  bool ShowCellsInfo;

  JDivideCpu(JLog *log,string dirout);
  ~JDivideCpu();
  void Reset();
  void OutClear(){ OutCount=0; }
  
  void GetDataOut(unsigned* id,tfloat3* pos,tfloat3* vel,float* rhop,bool clearout);
  
  void Config(float sizecell,float incz,unsigned np,unsigned nfixed,unsigned nmoving,unsigned nfloat,unsigned nfluid,const tfloat3 *pos,float dp);
  void ConfigRhopOut(bool rhopout,float rhopoutmin,float rhopoutmax){ RhopOut=rhopout; RhopOutMin=rhopoutmin; RhopOutMax=rhopoutmax; }
  void ConfigMemory();
  unsigned int GetMemoryAlloc()const;

  void DivBoundary(const tfloat3 *pos);
  bool DivFluid(const unsigned* idp,const tfloat3* pos,const tfloat3* vel,const float* rhop,float time);
   
  void RevDivide(bool sorted,tfloat3 *pos,int *idp);

  void SortBoundary(unsigned *vec);
  void SortBoundary(float *vec);
  void SortBoundary(tfloat3 *vec);
  void SortFluid(unsigned *vec);
  void SortFluid(float *vec);
  void SortFluid(tfloat3 *vec);
  void SortOnlyFluid(tsymatrix3f *vec);    
  
  bool CellNoEmpty(int box,byte kind)const;
  int CellBegin(int box,byte kind)const;
  int CellSize(int box,byte kind)const;

  void SetLimit(tfloat3 posmin,tfloat3 posmax);
  tfloat3 GetLimitMin()const{ return(PosMin); }      ///<Returns minimum limits of the domain.
  float GetLimitIncZ()const{ return(IncZ); }         ///<Returns the maximum increase in Z direction.

  tint3 GetNcells()const{ return(Ncells); }          ///<Returns the number of cells.
  unsigned GetNp()const{ return(Np); }               ///<Returns the total number of particles.
  unsigned GetNpOk()const{ return(Np-NfluidOut); }   ///<Returns the total number of particles without the excluded ones.
  unsigned GetNfluidOut()const{ return(NfluidOut); } ///<Returns the number of excluded particles from the total.
  unsigned GetOutCount()const{ return(OutCount); }   ///<Returns the number of excluded particles and stored in \ref Out.
  unsigned GetRhopOutCount(){ return(RhopOutCount); }///<Returns the number of excluded particles by \ref RhopOut.

  tfloat3 GetPosMin()const{ return(PosMin); }

  StCellInfo GetCellInfo();
  StCellInfo GetCellInfoData(unsigned n,unsigned nct,unsigned *beginb,unsigned *beginf)const;

  const unsigned* GetBeginBound()const{ return(BeginBound); }
  const unsigned* GetBeginFluid()const{ return(BeginFluid); }

  unsigned GetNctTotMax()const{ return(NctTotMax); } 
  unsigned GetNmoving()const{ return(Nmoving); }

  static StDivideInfo CalcDivideInfo(float scell,float border,float incz,unsigned np,const tfloat3* pos);

};

#endif







