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

/// \file JPartData.h \brief Declares the class \ref JPartData.

#ifndef _JPartData_
#define _JPartData_

#include <string>
#include "JObject.h"
#include "TypesDef.h"

//#define DBG_JPartData

//##############################################################################
//# JPartData
//##############################################################################
/// \brief Allows reading/writing files with data of particles in formats binx2,ascii...

class JPartData : protected JObject
{
public:
  ///Type of format files.
  typedef enum{ FmtNull=0,FmtBin=1,FmtBi2=2,FmtBi2Out=3,FmtAscii=4,FmtFlw=5}TpFmtFile; 
  typedef struct{  
    float dp,h,b,rhop0,gamma;
    float massbound,massfluid;
  }StConfig;

  ///Manages the information of the excluded particles.
  typedef struct{ 
    unsigned part;
    float timeout;
    unsigned id;
    tfloat3 pos;
    tfloat3 vel;
    float rhop;
  }StParticleOut; //-sizeof(40)

private:
  ///Structure that describes the header of binary format files.
  typedef struct{
    char titu[16];               ///<Title of the file "#File BINX-000".
    byte fmt;                    ///<File format.
    byte bitorder;               ///<1:BigEndian 0:LittleEndian.
    byte full;                   ///<1:All data of header and particles, 0:Without... [id,pos,vel] of fixed particles
    byte data2d;                 ///<1:Data for a 2D case, 0:3D Case.
  }StHeadFmtBin;//-sizeof(20)  

  //-Structures to be used with the format BINX-011:
  typedef struct{  
    byte fmt;                    ///<File format.
    byte bitorder;               ///<1:BigEndian 0:LittleEndian.
    bool full;                   ///<1:All data of header and particles, 0:Without... [id,pos,vel] of fixed particles
    float h,dp,b,rhop0,gamma;
    float massbound,massfluid;
    unsigned np,nb,nbf,npok;
    float time;
  }StInfoFileBin;
  typedef struct{//-They must be all of 4 bytes due to conversion ByteOrder ...
    float h,dp,b,rhop0,gamma;
    float massbound,massfluid;
    int np,nb,nbf,npok;
    float timestep;
  }StHeadDatFullBin;
  typedef struct{ //-They must be all of 4 bytes due to conversion ByteOrder ...
    int np,nb,nbf,npok;
    float timestep;
  }StHeadDatBin;

  //-Structures to be used with the format BINX-020:
  typedef struct{ 
    byte fmt;                    ///<File format.
    byte bitorder;               ///<1:BigEndian 0:LittleEndian.
    bool full;                   ///<1:All data of header and particles, 0:Without... [id,pos,vel] of fixed particles
    byte data2d;                 ///<1:Data for a 2D case, 0:3D Case.
    float dp,h,b,rhop0,gamma;
    float massbound,massfluid;
    unsigned np,nfixed,nmoving,nfloat,nfluidout;
    float time;
  }StInfoFileBi2;
  typedef struct{  //-They must be all of 4 bytes due to conversion ByteOrder ...
    float h,dp,b,rhop0,gamma;
    float massbound,massfluid;
    int np,nfixed,nmoving,nfloat,nfluidout;
    float timestep;
  }StHeadDatFullBi2;//-sizeof(52)  
  typedef struct{ //-They must be all of 4 bytes due to conversion ByteOrder ...
    int np,nfixed,nmoving,nfloat,nfluidout;
    float timestep;
  }StHeadDatBi2;
  
  std::string MaskFileName[5];   ///<Masks to create filenames in different formats.

  TpFmtFile FmtData;             ///<Current format of data. 
  // - FmtNull: Undetermined format...
  // - FmtAscii: There are no excluded particles and data of all particles 
  //   are in the main arrays, sorted by id.
  // - FmtBin: Data of no excluded particles are only stored in the main arrays (beginning), 
  //   particles can be disorderd but always boundaries first and then fluid particles.
  //   Data of excluded particles are stored in Out[].
  //   Using SortDataOut() data of particles Out can be located in the end of the main arrays.
  // - FmtBi2: Data of no excluded particles are only stored in the main arrays (beginning),
  //   particles are sorted by id but keeping the positions of the excluded particles.
  //   Data of excluded particles can be in Out[] or not.
  //   Using SortDataOut() data of particles Out can be located in the end of the main arrays.

  unsigned ByteOrderFile;        ///<Order of bytes of the last loaded binary file. 
  float PartTime;                ///<Instant of simulation.
  unsigned PartNumber;           ///<Number of files of data.

  bool DataConfig;               ///<Indicates whether the configuration data is loaded.
  bool DataBound;                ///<Indicates whether the boundary data is loaded.

  StConfig ConfigInfo;           ///<Basic parameters of the simulation configuration.
  
  float TimeByPart;              ///<Time per part to calculate PartTime in ASCII format. 
  
  bool Data2D;                   ///<Indicates whether data comes from a 2D case. 

  unsigned Np;                   ///<Number of total particles. 
  unsigned Nfixed;               ///<Number of fixed boundary particles. 
  unsigned Nmoving;              ///<Number of moving boundary particles. 
  unsigned Nfloat;               ///<Number of floating boundary particles.
  unsigned Nbound;               ///<Number of boundary particles ( \ref Nfixed + \ref Nmoving + \ref Nfloat ). 
  unsigned Nfluid;               ///<Number of fluid particles (including the excluded ones).  
  unsigned NfluidOut;            ///<Number of fluid particles excluded for the current \ref PartNumber.
  unsigned Nprobe;          ///<Number of probe particles.


  unsigned *Id;                  ///<ID of particles.
  tfloat3 *Pos;                  ///<Position of particles.
  tfloat3 *Vel;                  ///<Velocity of particles.
  float *Rhop;                   ///<Density of particles.

  tfloat3 *ProbeVel;         ///<Velocity of the probe particles (X,Y,Z).
  float *ProbeRhop;          ///<Density of the probe particles.

  int *OrderOut;                 ///<Only for excluded particles, indicates the position of particle \ref Out (-1 if it is not an excluded particle).

  unsigned OutSize;              ///<Number of particles that can be allocated in \ref Out.
  unsigned OutCount;             ///<Number of excluded particles stored in \ref Out, that can be larger than \ref NfluidOut.
  StParticleOut* Out;            ///<List of excluded particles by order of exclusion.
  unsigned OutCountSaved;        ///<Number of excluded particles and stored in the last call to \ref SaveFileBi2Out.


  void SetConfigInfo(float dp,float h,float b,float rhop0,float gamma,float massbound,float massfluid){
    ConfigInfo.dp=dp; ConfigInfo.h=h; ConfigInfo.b=b; ConfigInfo.rhop0=rhop0; ConfigInfo.gamma=gamma; ConfigInfo.massbound=massbound; ConfigInfo.massfluid=massfluid;
  }
  bool CheckConfigInfo(float dp,float h,float b,float rhop0,float gamma,float massbound,float massfluid)const{
    return(ConfigInfo.dp==dp&&ConfigInfo.h==h&&ConfigInfo.b==b&&ConfigInfo.rhop0==rhop0&&ConfigInfo.gamma==gamma&&ConfigInfo.massbound==massbound&&ConfigInfo.massfluid==massfluid);
  }
  StInfoFileBin GetInfoFileBin(std::ifstream &pf)const;
  StInfoFileBi2 GetInfoFileBi2(std::ifstream &pf,bool fileout)const;
  void SetOutSize(unsigned size);
  void FluidOutClear(bool clearout);

public:

  JPartData();
  ~JPartData();
  void Reset();
  void Config(TpFmtFile fmt,unsigned np,unsigned nbound,unsigned nfluid,unsigned nfixed,unsigned nmoving,unsigned nfloat,float dp,float h,float b,float rhop0,float gamma,float massbound,float massfluid,bool data2d=false,unsigned nprobe=0);

  unsigned SetDataUnsorted(unsigned part,float timepart,bool outreset,unsigned npok,unsigned nout,unsigned* id,tfloat3* pos,tfloat3* vel,float* rhop,tfloat3* probevel=NULL,float* proberhop=NULL);

  unsigned GetDataSort(unsigned size,unsigned *id,tfloat3 *pos,tfloat3 *vel,float *rhop,bool full)const;
  void GetDataPointers(unsigned* &id,tfloat3* &pos,tfloat3* &vel,float* &rhop){ id=Id; pos=Pos; vel=Vel; rhop=Rhop; }
  void GetDataPointers(unsigned* &id,tfloat3* &pos,tfloat3* &vel,float* &rhop,tfloat3* &probevel,float* &proberhop){ id=Id; pos=Pos; vel=Vel; rhop=Rhop; probevel=ProbeVel; proberhop=ProbeRhop; }
  
  void SortDataOut();

  TpFmtFile GetFmtData()const{ return(FmtData); }

  void LoadFile(TpFmtFile fmt,unsigned part,const std::string &dir="");
  void SaveFile(TpFmtFile fmt,const std::string &dir="",bool updatebi2out=false);

  void LoadFileAscii(unsigned part,const std::string &file);
  void SaveFileAscii(std::string file="")const;

  void LoadFileBin(unsigned part,const std::string &file);
  void SaveFileBin(bool full=false,std::string file="")const;

  void LoadFileBi2(unsigned part,const std::string &file);
  void SaveFileBi2(bool full=false,std::string file="")const;
  void LoadFileBi2Out(std::string file);
  void SaveFileBi2Out(std::string file,bool updatebi2out=false);

  void LoadFileInfo(TpFmtFile fmt,const std::string &file);
  void SaveFileInfo(const std::string &file)const;
  void SaveFileInfoDat(const std::string &file)const;

  void SaveFileFlw(std::string file="")const;

  void SetMaskFileName(TpFmtFile fmt,const std::string &mask){ if(1<=int(fmt)&&int(fmt)<=5)MaskFileName[int(fmt)-1]=mask; }
  std::string GetMaskFileName(TpFmtFile fmt)const{ return(1<=int(fmt)&&int(fmt)<=5? MaskFileName[int(fmt)-1]: std::string("")); }
  std::string GetFileName(TpFmtFile fmt,unsigned part,const std::string &dir)const; 
  void SetPartNumber(unsigned part){ PartNumber=part; }
  unsigned GetPartNumber()const{ return(PartNumber); }
  void SetPartTime(float parttime){ PartTime=parttime; }
  float GetPartTime()const{ return(PartTime); }
  void SetTimeByPart(float timebypart){ TimeByPart=timebypart; }
  float GetTimeByPart()const{ return(TimeByPart); }

  unsigned GetByteOrderFile()const{ return(ByteOrderFile); }
  StConfig GetConfigInfo()const{ return(ConfigInfo); }
  bool DataConfigOk()const{ return(DataConfig); }
  bool DataBoundOk()const{ return(DataBound); }

  bool GetData2D()const{ return(Data2D); }

  unsigned GetNp()const{        return(Np); }
  unsigned GetNbound()const{    return(Nbound); }
  unsigned GetNfluid()const{    return(Nfluid); }
  unsigned GetNfixed()const{    return(Nfixed); }
  unsigned GetNmoving()const{   return(Nmoving); }
  unsigned GetNfloat()const{    return(Nfloat); }
  unsigned GetNfluidOut()const{ return(NfluidOut); }
  unsigned GetNprobe()const{return(Nprobe); }

  unsigned GetMemoryAlloc()const;

  StInfoFileBin GetInfoFromFileBin(const std::string &file)const;
  StInfoFileBi2 GetInfoFromFileBi2(const std::string &file)const;  

  void DgCompare(const JPartData& pdat);

};


#endif





