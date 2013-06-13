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

/// \file JPartData.cpp \brief Implements the class \ref JPartData

#include "JPartData.h"
#include "JVarsAscii.h"

#include <fstream>
#include <cmath>
#include <cstring>
//#include <iostream>
//#include <sstream>
//#include <cstdlib>
//#include <ctime>
//using namespace std;
#include "Functions.h"

using std::string;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::fstream;
using std::min;
using std::max;
using std::endl;


//==============================================================================
/// Constructor.
//==============================================================================
JPartData::JPartData(){
  ClassName="JPartData";
  Id=NULL; Pos=NULL; Vel=NULL; Rhop=NULL; OrderOut=NULL;
  Out=NULL;
  SetMaskFileName(FmtBin,"PartBinx_%04d.bin");
  SetMaskFileName(FmtBi2,"Part%04d.bi2");
  SetMaskFileName(FmtBi2Out,"PartOut.bi2");
  SetMaskFileName(FmtAscii,"PART_%04d");
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartData::~JPartData(){
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartData::Reset(){
  delete[] Id;        Id=NULL; 
  delete[] Pos;       Pos=NULL;
  delete[] Vel;       Vel=NULL;
  delete[] Rhop;      Rhop=NULL;
  delete[] OrderOut;  OrderOut=NULL;
  delete[] Out;       Out=NULL;
  OutSize=0; OutCount=0; OutCountSaved=0;
  FmtData=FmtNull;
  SetPartTime(0); SetPartNumber(0);
  TimeByPart=1;
  DataConfig=false; DataBound=false;
  Data2D=false;
  Np=0;
  Nbound=0; Nfluid=0;
  NfluidOut=0;
  Nfixed=0; Nmoving=0; Nfloat=0;
  memset(&ConfigInfo,0,sizeof(StConfig));
}

//==============================================================================
/// Sets the basic setup.
//==============================================================================
void JPartData::Config(TpFmtFile fmt,unsigned np,unsigned nbound,unsigned nfluid,unsigned nfixed,unsigned nmoving,unsigned nfloat,float dp,float h,float b,float rhop0,float gamma,float massbound,float massfluid,bool data2d){
  const char met[]="Config";
  Reset();
  if(nbound+nfluid!=np||nfixed+nmoving+nfloat!=nbound)RunException(met,"Error in the number of particles of each type.");
  FmtData=fmt;
  Data2D=(fmt==FmtBi2||fmt==FmtBi2Out||fmt==FmtAscii? data2d: false);
  SetConfigInfo(dp,h,b,rhop0,gamma,massbound,massfluid);
  Np=np; Nbound=nbound; Nfluid=nfluid;
  Nfixed=nfixed; Nmoving=nmoving; Nfloat=nfloat;
  try{
    Id=new unsigned[Np];        
    Pos=new tfloat3[Np];       
    Vel=new tfloat3[Np];  
    Rhop=new float[Np];       
    OrderOut=new int[Nfluid]; 
  }
  catch(const std::bad_alloc){
    RunException(met,"Cannot allocate the requested memory.");
  }
  if(fmt==FmtBi2)for(unsigned c=0;c<Np;c++)Id[c]=c;
  FluidOutClear(true);
  DataConfig=true;
}

//==============================================================================
/// Returns the amount of allocated memory.
//==============================================================================
unsigned JPartData::GetMemoryAlloc()const{
  unsigned m=0;
  if(Id)m+=Np*sizeof(unsigned);
  if(Pos)m+=Np*sizeof(tfloat3);
  if(Vel)m+=Np*sizeof(tfloat3);
  if(Rhop)m+=Np*sizeof(float);
  if(OrderOut)m+=Nfluid*sizeof(int);
  if(Out)m+=OutSize*sizeof(StParticleOut);
  return(m);
}

//==============================================================================
/// Modifies particle data. 
/// Returns true if any particle was not excluded, when previous data were.
//==============================================================================
unsigned JPartData::SetDataUnsorted(unsigned part,float timepart,bool outreset,unsigned npok,unsigned nout,unsigned* id,tfloat3* pos,tfloat3* vel,float* rhop){
  const char met[]="SetDataUnsorted";
  unsigned errout=0;
  if(!DataConfigOk())RunException(met,"Initial configuration missing.");
  if(NfluidOut!=OutCount)RunException(met,"number of excluded particles does not match with stored number.");
  if(npok+nout>Np)RunException(met,"The number of partibles is larger than the maximum set.");
  if(!id||!pos||!vel||!rhop)RunException(met,"Particle data missing.");
  SetPartNumber(part);
  SetPartTime(timepart);
  if(outreset)FluidOutClear(true);
  if(GetFmtData()==FmtBi2){
    const unsigned npokold=Np-NfluidOut;
    //-Generates index of particles in id [] from particle identifier.
    int *rid=new int[Np];
    for(unsigned c=0;c<Np;c++)rid[c]=-1;
    const unsigned nn=npok+nout;
    for(unsigned c=0;c<nn;c++)rid[id[c]]=c;
    //-Modifies boundary particles.
    for(unsigned p=0;p<Nbound;p++){
      int pid=rid[p]; 
      if(pid>=0){ Pos[p]=pos[pid]; Vel[p]=vel[pid]; Rhop[p]=rhop[pid]; }
    }
#ifdef DBG_JPartData
    //-Registers new partices out (ORDERING).
    for(unsigned cp=Nbound;cp<Np;cp++){
      int c=rid[cp]; 
      if(c>=int(npok)){
        unsigned p=id[c];
        if(OrderOut[p-Nbound]==-1){//-No excluded particles.
          if(OutCount>=int(OutSize))SetOutSize(min(max(unsigned(1000),unsigned(OutSize*2)),Nfluid));
          StParticleOut* out=Out+OutCount;
          out->part=GetPartNumber();
          out->id=p;
          out->pos=pos[c];
          out->vel=vel[c];
          out->rhop=rhop[c];
          out->timeout=GetPartTime();
          OrderOut[p-Nbound]=OutCount;
          OutCount++; NfluidOut++;
        }
      }
    }   
#endif
    //-Registers new partices out (WITHOUT ORDERING).
    for(unsigned c=npok;c<nn;c++){
      unsigned p=id[c];
      if(OrderOut[p-Nbound]==-1){//-No excluded particles.
        if(OutCount>=int(OutSize))SetOutSize(min(max(unsigned(1000),unsigned(OutSize*2)),Nfluid));
        StParticleOut* out=Out+OutCount;
        out->part=GetPartNumber();
        out->id=p;
        out->pos=pos[c];
        out->vel=vel[c];
        out->rhop=rhop[c];
        out->timeout=GetPartTime();
        OrderOut[p-Nbound]=OutCount;
        OutCount++; NfluidOut++;
      }
    }
    //-Relocates and modifies fluid particles.
    unsigned pok=Nbound,pokold=Nbound;
    for(unsigned p=Nbound;p<Np;p++){
      const int pid=rid[p];
      if(OrderOut[p-Nbound]==-1){//-No excluded particles.
        if(pid>=0){
          Id[pok]=p; Pos[pok]=pos[pid]; Vel[pok]=vel[pid]; Rhop[pok]=rhop[pid]; pok++;
        }
        else{
          for(;pokold<npokold&&Id[pokold]<p;pokold++);
          if(pokold>=npokold||Id[pokold]!=p)RunException(met,"Cannot find the data of a supposed non-excluded particle.");
          if(pokold!=pok){ Id[pok]=p; Pos[pok]=Pos[pokold]; Vel[pok]=Vel[pokold]; Rhop[pok]=Rhop[pokold]; }
          pok++;
        } 
      } 
      else if(pid>=0&&pid<int(npok))errout++;//-Excluded particle that now it was not...
    }
    delete[] rid;
    if(pok!=Np-NfluidOut)RunException(met,"errr pok \n");
#ifdef DBG_JPartData
    for(unsigned p=0;p<Nbound;p++)if(Id[p]!=p)RunException(met,"Error in Id[] for boundary.");
    unsigned xid=Nbound;
    for(unsigned p=Nbound;p<Np-NfluidOut;p++){
      for(;OrderOut[xid-Nbound]>=0&&xid<Np;xid++);
      if(Id[p]!=xid)RunException(met,"Error in Id[] for fluids.");
      else xid++;
    } 
#endif
  }
  else RunException(met,"This method is not implmentd for the format BINX-20.");
  DataBound=true;
  return(errout);
}

//==============================================================================
/// Returns the complete filename according to format, part and directory.
//==============================================================================
std::string JPartData::GetFileName(TpFmtFile fmt,unsigned part,const std::string &dir)const{
  std::string file;
  std::string mask=GetMaskFileName(fmt);
  if(!mask.empty()){
    char cad[1024];
    sprintf(cad,mask.c_str(),part);
    file=fun::GetDirWithSlash(dir)+cad;
  }
  return(file);
}

//==============================================================================
/// Loads data from a file.
//==============================================================================
void JPartData::LoadFile(TpFmtFile fmt,unsigned part,const std::string &dir){
  switch(fmt){
    case FmtBin:    LoadFileBin(part,GetFileName(fmt,part,dir));    break;
    case FmtBi2:    LoadFileBi2(part,GetFileName(fmt,part,dir));    break;
    case FmtBi2Out: LoadFileBi2Out(GetFileName(fmt,part,dir));      break;
    case FmtAscii:  LoadFileAscii(part,GetFileName(fmt,part,dir));  break;
  }
}

//==============================================================================
/// Stores data in a file.
//==============================================================================
void JPartData::SaveFile(TpFmtFile fmt,const std::string &dir,bool updatebi2out){
  switch(fmt){
    case FmtBin:    SaveFileBin(false,GetFileName(fmt,GetPartNumber(),dir));            break;
    case FmtBi2:    SaveFileBi2(false,GetFileName(fmt,GetPartNumber(),dir));            break;
    case FmtBi2Out: SaveFileBi2Out(GetFileName(fmt,GetPartNumber(),dir),updatebi2out);  break;
    case FmtAscii:  SaveFileAscii(GetFileName(fmt,GetPartNumber(),dir));                break;
  }
}

//==============================================================================
/// Loads data from a file with ASCII format.
//==============================================================================
void JPartData::LoadFileAscii(unsigned part,const std::string &file){
  const char met[]="LoadFileAscii";
  if(FmtData!=FmtAscii)Reset();
  //-Gets basic configuration.
  if(!DataConfigOk()){
    string dirfile=fun::GetDirWithSlash(fun::GetDirParent(file));
    string fbase=(!part? file: GetFileName(FmtAscii,0,dirfile));
    string fileinfo=fbase+".inf";
    LoadFileInfo(FmtAscii,fileinfo);
  }
  //-Gets general data from a file.
  SetPartTime(GetTimeByPart()*part);
  SetPartNumber(part);
  DataBound=true;
  ifstream pf;
  pf.open(file.c_str());
  if(pf){
    float fpx,fpy,fpz,fvx,fvy,fvz,frhop,fpres,fmass;
    if(Data2D)for(unsigned c=0;c<Np;c++){ //-2D data.
      pf >> fpx >> fpz >> fvx >> fvz >> frhop >> fpres >> fmass;
      Id[c]=c;  Pos[c]=TFloat3(fpx,0,fpz);  Vel[c]=TFloat3(fvx,0,fvz);  Rhop[c]=frhop;
    } 
    else for(unsigned c=0;c<Np;c++){ //-3D data. 
      pf >> fpx >> fpy >> fpz >> fvx >> fvy >> fvz >> frhop >> fpres >> fmass;
      Id[c]=c;  Pos[c]=TFloat3(fpx,fpy,fpz);  Vel[c]=TFloat3(fvx,fvy,fvz);  Rhop[c]=frhop;
    } 
    if(pf.fail())RunException(met,"Error reading data file.",file);
    pf.close();
    //-Clears the list of paticles out.
    if(NfluidOut||OutCount)FluidOutClear(true);
  }
  else RunException(met,"Cannot open the file.",file);
}
//==============================================================================
/// Stores data in a file with ASCII format.
//==============================================================================
void JPartData::SaveFileAscii(std::string file)const{
  const char met[]="SaveFileAscii";
  if(!DataConfigOk())RunException(met,"Initial configuration missing.");
  if(!DataBoundOk())RunException(met,"Boundary particles data missing.");
  if(NfluidOut>OutCount)RunException(met,"Excluded particles data missing.");
  if(file.empty())file=GetFileName(FmtAscii,GetPartNumber(),"");
  ofstream pf;
  pf.open(file.c_str());
  if(pf){
    char cad[4096];
    char fmt3d[]=" %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E";
    char fmt2d[]=" %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E %15.7E";
    int *rid=NULL;
    if(FmtData==FmtBin||FmtData==FmtBi2)rid=new int[Np];
    else if(FmtData!=FmtAscii)RunException(met,"Unrecognised data format.");
    if(rid){
      unsigned npok=Np-NfluidOut;
      for(unsigned c=0;c<npok;c++)rid[Id[c]]=c;
      for(unsigned c=0;c<NfluidOut;c++)rid[Out[c].id]=-int(c+1);
    }
    for(unsigned p=0;p<Np;p++){
      int c=(rid? rid[p]: int(p));
      if(c>=0){//-Particles in.
        float pres=ConfigInfo.b*(pow(Rhop[c]/ConfigInfo.rhop0,ConfigInfo.gamma)-1.0f);
        if(Data2D)sprintf(cad,fmt2d,Pos[c].x,Pos[c].z,Vel[c].x,Vel[c].z,Rhop[c],pres,(p<Nbound? ConfigInfo.massbound: ConfigInfo.massfluid));
        else sprintf(cad,fmt3d,Pos[c].x,Pos[c].y,Pos[c].z,Vel[c].x,Vel[c].y,Vel[c].z,Rhop[c],pres,(p<Nbound? ConfigInfo.massbound: ConfigInfo.massfluid));
      }
      else{//-Particles out.
        StParticleOut* out=Out+(-c-1);
        float pres=ConfigInfo.b*(pow(out->rhop/ConfigInfo.rhop0,ConfigInfo.gamma)-1.0f);
        if(Data2D)sprintf(cad,fmt2d,out->pos.x,out->pos.z,out->vel.x,out->vel.z,out->rhop,pres,ConfigInfo.massfluid);
        else sprintf(cad,fmt3d,out->pos.x,out->pos.y,out->pos.z,out->vel.x,out->vel.y,out->vel.z,out->rhop,pres,ConfigInfo.massfluid);      
      }
      pf << cad << endl;
    }
    if(pf.fail())RunException(met,"File writing failure.",file);
    pf.close();
    delete[] rid;
  }
  else RunException(met,"Cannot open the file.",file);
  //-Generates .inf and .data file with meta-information of particles.
  if(GetPartNumber()==0){
    int pos=(int)file.find_last_of(".");
    SaveFileInfo(file.substr(0,pos)+".inf");
    SaveFileInfoDat(file.substr(0,pos)+".dat");
  }
}

//==============================================================================
/// Reads and returns information from file BINX-011.
//==============================================================================
JPartData::StInfoFileBin JPartData::GetInfoFileBin(ifstream &pf)const{
  StInfoFileBin ret;
  memset(&ret,0,sizeof(StInfoFileBin));
  StHeadFmtBin hfmt;
  pf.read((char*)&hfmt,sizeof(StHeadFmtBin));
  if(!strcmp(hfmt.titu,"#File BINX-011 ")){
    ret.fmt=11;
    ret.bitorder=hfmt.bitorder;
    ret.full=hfmt.full!=0;
    bool rendi=(hfmt.bitorder!=byte(fun::GetByteOrder()));
    if(hfmt.full){
      StHeadDatFullBin hdat;
      pf.read((char*)&hdat,sizeof(StHeadDatFullBin));
      if(rendi)fun::ReverseByteOrder((int*)&hdat,sizeof(StHeadDatFullBin)/4);//-Conversion Big/LittleEndian
      ret.dp=hdat.dp;
      ret.h=hdat.h;
      ret.b=hdat.b;
      ret.rhop0=hdat.rhop0;
      ret.gamma=hdat.gamma;
      ret.massbound=hdat.massbound;
      ret.massfluid=hdat.massfluid;
      ret.np=hdat.np;
      ret.nb=hdat.nb;
      ret.nbf=hdat.nbf;
      ret.npok=hdat.npok;
      ret.time=hdat.timestep;
    }
    else{
      StHeadDatBin hd;
      pf.read((char*)&hd,sizeof(StHeadDatBin));
      if(rendi)fun::ReverseByteOrder((int*)&hd,sizeof(StHeadDatBin)/4);//-Conversion Big/LittleEndian
      ret.np=hd.np;
      ret.nb=hd.nb;
      ret.nbf=hd.nbf;
      ret.npok=hd.npok;
      ret.time=hd.timestep;
    }
  }
  return(ret);
}
//==============================================================================
/// Loads data from a file with BINX-11 format. 
//==============================================================================
void JPartData::LoadFileBin(unsigned part,const std::string &file){
  const char* met="LoadFileBin";
  if(!part||FmtData!=FmtBin)Reset();
  ifstream pf;
  pf.open(file.c_str(),ios::binary|ios::in);
  if(pf){
    //-Checks file format and configures object.
    StInfoFileBin inf=GetInfoFileBin(pf);
    if(inf.fmt!=11)RunException(met,"Invalid file format.",file);
    if(!DataConfigOk())Config(FmtBin,inf.np,inf.nb,inf.np-inf.nb,inf.nbf,inf.nb-inf.nbf,0,inf.dp,inf.h,inf.b,inf.rhop0,inf.gamma,inf.massbound,inf.massfluid,false);
    else{  
      if(Np!=inf.np||Nbound!=inf.nb||Nfixed!=inf.nbf||Nmoving!=inf.nb-inf.nbf||Nfloat!=0)RunException(met,"The number of particles does not match the loaded data.",file);
      if(inf.full&&!CheckConfigInfo(inf.dp,inf.h,inf.b,inf.rhop0,inf.gamma,inf.massbound,inf.massfluid))RunException(met,"The configuration parameter does not match the loaded data.",file);
    }
    if(!DataBoundOk()&&!inf.full)RunException(met,"Fixed boundary data missing.",file);
    //-Gets general data of a file.
    ByteOrderFile=unsigned(inf.bitorder? fun::BigEndian: fun::LittleEndian);
    SetPartTime(inf.time);
    SetPartNumber(part);
    //-Gets particle data of a file.
    bool rendi=(inf.bitorder!=byte(fun::GetByteOrder()));
    if(inf.full){
      if(Id)pf.read((char*)Id,sizeof(int)*inf.np);       else pf.seekg(sizeof(int)*inf.np,ios::cur);
      if(Pos)pf.read((char*)Pos,sizeof(tfloat3)*inf.np); else pf.seekg(sizeof(tfloat3)*inf.np,ios::cur);
      if(Vel)pf.read((char*)Vel,sizeof(tfloat3)*inf.np); else pf.seekg(sizeof(tfloat3)*inf.np,ios::cur);
      if(Rhop)pf.read((char*)Rhop,sizeof(float)*inf.np); else pf.seekg(sizeof(float)*inf.np,ios::cur);
      if(rendi){
        if(Id)fun::ReverseByteOrder((int*)Id,inf.np);
        if(Pos)fun::ReverseByteOrder((int*)Pos,inf.np*3);
        if(Vel)fun::ReverseByteOrder((int*)Vel,inf.np*3);
        if(Rhop)fun::ReverseByteOrder((int*)Rhop,inf.np);
      }
      DataBound=true;
    }
    else{
      unsigned nf=inf.np-inf.nb;
      if(Id)pf.read((char*)(Id+inf.nb),sizeof(int)*nf);       else pf.seekg(sizeof(int)*nf,ios::cur);
      if(Pos)pf.read((char*)(Pos+inf.nb),sizeof(tfloat3)*nf); else pf.seekg(sizeof(tfloat3)*nf,ios::cur);
      if(Vel)pf.read((char*)(Vel+inf.nb),sizeof(tfloat3)*nf); else pf.seekg(sizeof(tfloat3)*nf,ios::cur);
      if(Rhop)pf.read((char*)Rhop,sizeof(float)*inf.np);      else pf.seekg(sizeof(float)*inf.np,ios::cur);
      if(rendi){
        if(Id)fun::ReverseByteOrder((int*)(Id+inf.nb),nf);
        if(Pos)fun::ReverseByteOrder((int*)(Pos+inf.nb),nf*3);
        if(Vel)fun::ReverseByteOrder((int*)(Vel+inf.nb),nf*3);
        if(Rhop)fun::ReverseByteOrder((int*)Rhop,inf.np);
      }
    }
    pf.close();
    //-Checks if any particle cease to be excluded and excludes it again
    if(1){
      for(unsigned c=Nbound;c<inf.npok;c++)if(OrderOut[Id[c]-Nbound]>=0){
        printf("The particle %d was OUT and now it is not...\n",Id[c]);
        inf.npok--;
        //-Particle to the end to be converted again in Out.
        unsigned id0=Id[inf.npok];  Id[inf.npok]=Id[c];     Id[c]=id0;
        tfloat3 pos0=Pos[inf.npok];  Pos[inf.npok]=Pos[c];   Pos[c]=pos0;
        tfloat3 vel0=Vel[inf.npok];  Vel[inf.npok]=Vel[c];   Vel[c]=vel0;
        float rhop0=Rhop[inf.npok]; Rhop[inf.npok]=Rhop[c]; Rhop[c]=rhop0;
      }
    }
    //-Processes excluded particles.
    if(inf.npok!=inf.np){
      int nout=inf.np-inf.npok;
      if(int(OutSize)<nout)SetOutSize(min(max(unsigned(1000),unsigned(nout*2)),Nfluid));
      for(unsigned c=inf.np-1;c>=inf.npok;c--){
        unsigned id=Id[c];
        if(OrderOut[id-Nbound]<0){
          StParticleOut* out=Out+NfluidOut;
          out->part=part;
          out->id=id;
          out->pos=Pos[c];
          out->vel=Vel[c];
          out->rhop=Rhop[c];
          out->timeout=GetPartTime();
          OrderOut[id-Nbound]=NfluidOut;
          NfluidOut++; OutCount=NfluidOut;
        }
      }
    }
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Places data of excluded particles at the end of the main arrays.
//==============================================================================
void JPartData::SortDataOut(){
  const char met[]="SortDataOut";
  if(FmtData!=FmtBin&&FmtData!=FmtBi2)RunException(met,"Method applicable only to formats BINX-11 and BINX-20.");
  if(!DataConfigOk())RunException(met,"Initial confoguration missing.");
  if(!DataBoundOk())RunException(met,"Boundary particle data missing.");
  if(NfluidOut>OutCount)RunException(met,"Excluded particle data missing.");
  for(unsigned c=0;c<NfluidOut;c++){//-Coloca datos de las particulas excluidas en vectores principales.
    StParticleOut* out=Out+c;
    unsigned idp=Np-c-1;
    Id[idp]=out->id; Pos[idp]=out->pos; Vel[idp]=out->vel; Rhop[idp]=out->rhop;
  }
}

//==============================================================================
/// Stores data in a file with format BINX-11.
//==============================================================================
void JPartData::SaveFileBin(bool full,std::string file)const{
  const char met[]="SaveFileBin";
  if(!DataConfigOk())RunException(met,"Initial confoguration missing.");
  if(!DataBoundOk())RunException(met,"Boundary particle data missing.");
  if(NfluidOut>OutCount)RunException(met,"Excluded particle data missing.");
  if(Nfloat)RunException(met,"Invalid format with floating bodies.");
  if(file.empty())file=GetFileName(FmtBin,GetPartNumber(),"");
  ofstream pf;
  pf.open(file.c_str(),ios::binary|ios::out);
  if(pf){
    StHeadFmtBin hfmt;
    strcpy(hfmt.titu,"#File BINX-011 ");
    hfmt.fmt=11;
    hfmt.bitorder=byte(fun::GetByteOrder());
    hfmt.full=(Nmoving>0||full||!GetPartNumber()? 1: 0);
    hfmt.data2d=0;
    pf.write((char*)&hfmt,sizeof(StHeadFmtBin));
    unsigned *id;
    tfloat3 *pos,*vel;
    float* rhop;
    if(FmtData==FmtAscii||FmtData==FmtBin||FmtData==FmtBi2){ 
      id=Id; pos=Pos; vel=Vel; rhop=Rhop;
      if(FmtData==FmtBin||FmtData==FmtBi2)for(unsigned c=0;c<NfluidOut;c++){//-Places data of excluded particles in the main arrays.
        StParticleOut* out=Out+c;
        unsigned idp=Np-c-1;
        Id[idp]=out->id; Pos[idp]=out->pos; Vel[idp]=out->vel; Rhop[idp]=out->rhop;
      }
    } 
    else RunException(met,"Unrecognised file format.");
    if(hfmt.full){
      StHeadDatFullBin hdata;
      hdata.h=ConfigInfo.h; hdata.dp=ConfigInfo.dp;  hdata.b=ConfigInfo.b;
      hdata.rhop0=ConfigInfo.rhop0;  hdata.gamma=ConfigInfo.gamma;
      hdata.massbound=ConfigInfo.massbound;
      hdata.massfluid=ConfigInfo.massfluid;
      hdata.np=Np;  hdata.nb=Nbound;  hdata.nbf=Nfixed;
      hdata.npok=Np-NfluidOut;  hdata.timestep=GetPartTime();
      pf.write((char*)&hdata,sizeof(StHeadDatFullBin));
      //-Writing data.
      pf.write((char*)id,sizeof(int)*hdata.np);
      pf.write((char*)pos,sizeof(tfloat3)*hdata.np);
      pf.write((char*)vel,sizeof(tfloat3)*hdata.np);
      pf.write((char*)rhop,sizeof(float)*hdata.np);
    }
    else{
      StHeadDatBin hd;
      hd.np=Np; hd.nb=Nbound; hd.nbf=Nfixed; hd.npok=Np-NfluidOut; hd.timestep=GetPartTime();
      pf.write((char*)&hd,sizeof(StHeadDatBin));
      //-Writing data.
      int nf=hd.np-hd.nb;
      pf.write((char*)(id+hd.nb),sizeof(int)*nf);
      pf.write((char*)(pos+hd.nb),sizeof(tfloat3)*nf);
      pf.write((char*)(vel+hd.nb),sizeof(tfloat3)*nf);
      pf.write((char*)rhop,sizeof(float)*hd.np);
    }
    if(pf.fail())RunException(met,"File writing failure.",file);
    pf.close();
    if(id!=Id)delete[] id;
    if(pos!=Pos)delete[] pos;
    if(vel!=Vel)delete[] vel;
    if(rhop!=Rhop)delete[] rhop;
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Reads and returns the infomration of a file BINX-020.
//==============================================================================
JPartData::StInfoFileBi2 JPartData::GetInfoFileBi2(ifstream &pf,bool fileout)const{
  StInfoFileBi2 ret;
  memset(&ret,0,sizeof(StInfoFileBi2));
  StHeadFmtBin hfmt;
  pf.read((char*)&hfmt,sizeof(StHeadFmtBin));
  if(!strcmp(hfmt.titu,"#File BINX-020 ")||(fileout&&hfmt.full==2)||(!fileout&&hfmt.full<=1)){
    ret.fmt=20;
    ret.bitorder=hfmt.bitorder;
    ret.full=(hfmt.full!=0);
    ret.data2d=(hfmt.data2d!=0);
    bool rendi=(hfmt.bitorder!=byte(fun::GetByteOrder()));
    if(hfmt.full){
      StHeadDatFullBi2 hdat;
      pf.read((char*)&hdat,sizeof(StHeadDatFullBi2));
      if(rendi)fun::ReverseByteOrder((int*)&hdat,sizeof(StHeadDatFullBi2)/4);//-Conversion Big/LittleEndian
      ret.dp=hdat.dp;
      ret.h=hdat.h;
      ret.b=hdat.b;
      ret.rhop0=hdat.rhop0;
      ret.gamma=hdat.gamma;
      ret.massbound=hdat.massbound;
      ret.massfluid=hdat.massfluid;
      ret.np=hdat.np;
      ret.nfixed=hdat.nfixed;
      ret.nmoving=hdat.nmoving;
      ret.nfloat=hdat.nfloat;
      ret.nfluidout=hdat.nfluidout;
      ret.time=hdat.timestep;
    }
    else{
      StHeadDatBi2 hd;
      pf.read((char*)&hd,sizeof(StHeadDatBi2));
      if(rendi)fun::ReverseByteOrder((int*)&hd,sizeof(StHeadDatBi2)/4);//-Conversion Big/LittleEndian
      ret.np=hd.np;
      ret.nfixed=hd.nfixed;
      ret.nmoving=hd.nmoving;
      ret.nfloat=hd.nfloat;
      ret.nfluidout=hd.nfluidout;
      ret.time=hd.timestep;
    }
  }
  return(ret);
}
//==============================================================================
/// Loads data from a file with format BINX-020.
//==============================================================================
void JPartData::LoadFileBi2(unsigned part,const std::string &file){
  const char* met="LoadFileBi2";
  if(FmtData!=FmtBi2)Reset();
  ifstream pf;
  pf.open(file.c_str(),ios::binary|ios::in);
  if(pf){
    //-Checks file format and configures object.
    StInfoFileBi2 inf=GetInfoFileBi2(pf,false);
    if(inf.fmt!=20)RunException(met,"Invalid file format.",file);
    unsigned nbound=inf.nfixed+inf.nmoving+inf.nfloat;
    unsigned nfluid=inf.np-nbound;
    if(!DataConfigOk())Config(FmtBi2,inf.np,nbound,nfluid,inf.nfixed,inf.nmoving,inf.nfloat,inf.dp,inf.h,inf.b,inf.rhop0,inf.gamma,inf.massbound,inf.massfluid,inf.data2d!=0);
    else{  
      if(Np!=inf.np||Nbound!=nbound||Nfixed!=inf.nfixed||Nmoving!=inf.nmoving||Nfloat!=inf.nfloat)RunException(met,"The particle number does not match with loaded data.",file);
      if(inf.full&&!CheckConfigInfo(inf.dp,inf.h,inf.b,inf.rhop0,inf.gamma,inf.massbound,inf.massfluid))RunException(met,"Configuration parameters does not match with loaded data.",file);
    }
    if(!DataBoundOk()&&!inf.full)RunException(met,"Fixed boundary data missing.",file);
    //-Obtains general data of the file.
    ByteOrderFile=unsigned(inf.bitorder? fun::BigEndian: fun::LittleEndian);
    SetPartTime(inf.time);
    SetPartNumber(part);
    //-Obtains particle data of the file.
    unsigned npok=inf.np-inf.nfluidout;
    bool rendi=(inf.bitorder!=byte(fun::GetByteOrder()));
    if(inf.full){
      if(Pos)pf.read((char*)Pos,sizeof(tfloat3)*npok); else pf.seekg(sizeof(tfloat3)*npok,ios::cur);
      if(Vel)pf.read((char*)Vel,sizeof(tfloat3)*npok); else pf.seekg(sizeof(tfloat3)*npok,ios::cur);
      if(Rhop)pf.read((char*)Rhop,sizeof(float)*npok); else pf.seekg(sizeof(float)*npok,ios::cur);
      if(rendi){
        if(Pos)fun::ReverseByteOrder((int*)Pos,npok*3);
        if(Vel)fun::ReverseByteOrder((int*)Vel,npok*3);
        if(Rhop)fun::ReverseByteOrder((int*)Rhop,npok);
      }
      DataBound=true;
    }
    else{
      unsigned nmov=npok-inf.nfixed;
      if(Pos)pf.read((char*)(Pos+inf.nfixed),sizeof(tfloat3)*nmov); else pf.seekg(sizeof(tfloat3)*nmov,ios::cur);
      if(Vel)pf.read((char*)(Vel+inf.nfixed),sizeof(tfloat3)*nmov); else pf.seekg(sizeof(tfloat3)*nmov,ios::cur);
      if(Rhop)pf.read((char*)Rhop,sizeof(float)*npok);              else pf.seekg(sizeof(float)*npok,ios::cur);
      if(rendi){
        if(Pos)fun::ReverseByteOrder((int*)(Pos+inf.nfixed),nmov*3);
        if(Vel)fun::ReverseByteOrder((int*)(Vel+inf.nfixed),nmov*3);
        if(Rhop)fun::ReverseByteOrder((int*)Rhop,npok);
      }
    }
    //-Recovers list of excluded particles.
    if(inf.nfluidout){
      //-Marks in orderout the excluded particles.
      byte *orderout=new byte[nfluid];
      memset(orderout,0,nfluid);
      bool idmode=inf.nfluidout<unsigned(nfluid/32);
      if(idmode){
        unsigned *idlis=new unsigned[inf.nfluidout];
        pf.read((char*)(idlis),sizeof(unsigned)*inf.nfluidout);
        if(rendi)fun::ReverseByteOrder((int*)idlis,inf.nfluidout); 
        for(unsigned c=0;c<inf.nfluidout;c++)orderout[idlis[c]-nbound]=1;
        delete[] idlis;
      }
      else{
        int size=int((nfluid+7)/8);
        byte *idbool=new byte[size];
        pf.read((char*)(idbool),size);
        for(int cbyte=0;cbyte<size;cbyte++)if(idbool[cbyte])
          for(int cbit=0;cbit<8;cbit++)if(idbool[cbyte]&(1<<cbit))orderout[cbyte*8+cbit]=1;
        delete[] idbool;
      }
      //-Stablishes id of valid fluid particles starting from orderout .
      int pos=nbound,cid=nbound;
      for(unsigned c=0;c<nfluid;c++){
        if(!orderout[c]){ Id[pos]=cid; cid++; pos++; }
        else cid++;
      }
      #ifdef DBG_JPartData
        if(pos!=npok||unsigned(cid)>Np)RunException(met,"Fails when adjusting Id[] for fluids.");
      #endif
      //-Updates content of OrderOut[]
      if(NfluidOut)FluidOutClear(false); //-Clears list of excluded particles.
      if(OutCount>=inf.nfluidout)for(unsigned c=0;c<inf.nfluidout;c++)OrderOut[Out[c].id-Nbound]=c; 
      else for(unsigned c=0;c<Nfluid;c++)if(orderout[c])OrderOut[c]=0;
      delete[] orderout;
    }
    else if(NfluidOut){
      FluidOutClear(false); //-Clears list of excluded particles.
      for(unsigned c=Nbound;c<Np;c++)Id[c]=c;
    }
    NfluidOut=inf.nfluidout;
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Stores data in a file with format BINX-020.
//==============================================================================
void JPartData::SaveFileBi2(bool full,std::string file)const{
  const char met[]="SaveFileBi2";
  if(!DataConfigOk())RunException(met,"Initial configuration missing.");
  if(!DataBoundOk())RunException(met,"boundary pasticles data missing.");
  if(file.empty())file=GetFileName(FmtBi2,GetPartNumber(),"");
  ofstream pf;
  pf.open(file.c_str(),ios::binary|ios::out);
  if(pf){
    StHeadFmtBin hfmt;
    strcpy(hfmt.titu,"#File BINX-020 ");
    hfmt.fmt=20;
    hfmt.bitorder=byte(fun::GetByteOrder());
    hfmt.full=(full||!GetPartNumber()? 1: 0);
    hfmt.data2d=GetData2D();
    pf.write((char*)&hfmt,sizeof(StHeadFmtBin));
    //-Prepares data to be stored.
    unsigned npok=Np-NfluidOut;
    tfloat3 *pos,*vel;
    float* rhop;
    if(FmtData==FmtBi2){ pos=Pos; vel=Vel; rhop=Rhop; }
    else if(FmtData==FmtAscii||FmtData==FmtBin){
      pos=new tfloat3[npok];
      vel=new tfloat3[npok];
      rhop=new float[npok];
      GetDataSort(npok,NULL,pos,vel,rhop,false);
    }
    else RunException(met,"Unrecognised data format.");
    //-Stores particle data in file. 
    if(hfmt.full){
      StHeadDatFullBi2 hdata;
      hdata.h=ConfigInfo.h; hdata.dp=ConfigInfo.dp;  hdata.b=ConfigInfo.b;
      hdata.rhop0=ConfigInfo.rhop0;  hdata.gamma=ConfigInfo.gamma;
      hdata.massbound=ConfigInfo.massbound;
      hdata.massfluid=ConfigInfo.massfluid;
      hdata.np=Np; hdata.nfixed=Nfixed; hdata.nmoving=Nmoving; hdata.nfloat=Nfloat;
      hdata.nfluidout=NfluidOut;  hdata.timestep=GetPartTime();
      pf.write((char*)&hdata,sizeof(StHeadDatFullBi2));
      //-Writing of data.
      pf.write((char*)pos,sizeof(tfloat3)*npok);
      pf.write((char*)vel,sizeof(tfloat3)*npok);
      pf.write((char*)rhop,sizeof(float)*npok);
    }
    else{
      StHeadDatBi2 hd;
      hd.np=Np; hd.nfixed=Nfixed; hd.nmoving=Nmoving; hd.nfloat=Nfloat;
      hd.nfluidout=NfluidOut;  hd.timestep=GetPartTime();
      pf.write((char*)&hd,sizeof(StHeadDatBi2));
      //-Writing of data.
      int nmov=npok-hd.nfixed;
      pf.write((char*)(pos+hd.nfixed),sizeof(tfloat3)*nmov);
      pf.write((char*)(vel+hd.nfixed),sizeof(tfloat3)*nmov);
      pf.write((char*)rhop,sizeof(float)*npok);
    }
    //-Stores list of excluded particles.
    if(NfluidOut){
      bool idmode=NfluidOut<unsigned(Nfluid/32);
      if(idmode){
        unsigned *idlis=new unsigned[NfluidOut];
        for(unsigned c=0;c<NfluidOut;c++)idlis[c]=Out[c].id;
        pf.write((char*)(idlis),sizeof(unsigned)*NfluidOut);
        delete[] idlis;
      }
      else{
        int size=int((Nfluid+7)/8);
        byte *idbool=new byte[size];
        memset(idbool,0,size);
        for(unsigned c=0;c<NfluidOut;c++){
          int idp=int(Out[c].id)-Nbound;
          idbool[int(idp/8)]|=byte(1<<(idp%8));
        }
        pf.write((char*)(idbool),size);
        delete[] idbool;
      }
    }
    if(pf.fail())RunException(met,"File writing failure.",file);
    pf.close();
    if(pos!=Pos)delete[] pos;
    if(vel!=Vel)delete[] vel;
    if(rhop!=Rhop)delete[] rhop;
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Loads data of excluded particles in file with format BINX-20.
//==============================================================================
void JPartData::LoadFileBi2Out(std::string file){
  const char* met="LoadFileBi2Out";
  if(FmtData!=FmtBi2)Reset();
  ifstream pf;
  pf.open(file.c_str(),ios::binary|ios::in);
  if(pf){
    //-Checks file format and configure object.
    StInfoFileBi2 inf=GetInfoFileBi2(pf,true);
    if(inf.fmt!=20)RunException(met,"Invalid file format.",file);
    unsigned nbound=inf.nfixed+inf.nmoving+inf.nfloat;
    unsigned nfluid=inf.np-nbound;
    if(!DataConfigOk())Config(FmtBi2,inf.np,nbound,nfluid,inf.nfixed,inf.nmoving,inf.nfloat,inf.dp,inf.h,inf.b,inf.rhop0,inf.gamma,inf.massbound,inf.massfluid,inf.data2d!=0);
    else{  
      if(Np!=inf.np||Nbound!=nbound||Nfixed!=inf.nfixed||Nmoving!=inf.nmoving||Nfloat!=inf.nfloat)RunException(met,"Particle number does not match the loaded data.",file);
      if(inf.full&&!CheckConfigInfo(inf.dp,inf.h,inf.b,inf.rhop0,inf.gamma,inf.massbound,inf.massfluid))RunException(met,"Configuration parameters does not match loaded data.",file);
    }
    //-Prepares Out[] and loads particles out.
    if(OutSize<inf.nfluidout){ OutCount=0; SetOutSize(inf.nfluidout); }
    OutCount=inf.nfluidout;
    pf.read((char*)(Out),sizeof(StParticleOut)*OutCount);
    if(inf.bitorder!=byte(fun::GetByteOrder()))fun::ReverseByteOrder((int*)Out,int(sizeof(StParticleOut)/4)*OutCount);
    pf.close();
    //-Updates content of OrderOut[].
    if(OutCount>=NfluidOut){
      unsigned nfluidout=NfluidOut;
      if(NfluidOut)FluidOutClear(false); //-Clears list of particles out.
      NfluidOut=nfluidout;
      for(unsigned c=0;c<NfluidOut;c++)OrderOut[Out[c].id-Nbound]=c; 
    }
  }
  else RunException(met,"Cannot open the file.",file);
}
//==============================================================================
/// Stores data of excluded particles in file with format BINX-20.
//==============================================================================
void JPartData::SaveFileBi2Out(std::string file,bool updatebi2out){
  const char met[]="SaveFileBi2Out";
  if(!DataConfigOk())RunException(met,"initial configuration missing.");
  if(!DataBoundOk())RunException(met,"Boundary particles data missing.");
  if(NfluidOut>OutCount)RunException(met,"Excluded particles data missing.");
  if(!updatebi2out||NfluidOut>OutCountSaved){
    if(file.empty())file=GetFileName(FmtBi2Out,GetPartNumber(),"");
    fstream pf;
    if(updatebi2out)pf.open(file.c_str(),ios::binary|ios::out|ios::in);
    else pf.open(file.c_str(),ios::binary|ios::out);
    if(pf){
      //-Writing data.
      pf.seekp(sizeof(StHeadFmtBin)+sizeof(StHeadDatFullBi2),ios::beg);
      if(updatebi2out){
        pf.seekp(sizeof(StParticleOut)*OutCountSaved,ios::cur);
        pf.write((char*)(Out+OutCountSaved),sizeof(StParticleOut)*(NfluidOut-OutCountSaved));
      }
      else pf.write((char*)(Out),sizeof(StParticleOut)*NfluidOut);
      pf.flush();
      //-Writing headers. 
      pf.seekp(0,ios::beg);
      StHeadFmtBin hfmt;
      strcpy(hfmt.titu,"#File BINX-020 ");
      hfmt.fmt=20;
      hfmt.bitorder=byte(fun::GetByteOrder());
      hfmt.full=2;
      hfmt.data2d=GetData2D();
      pf.write((char*)&hfmt,sizeof(StHeadFmtBin));
      StHeadDatFullBi2 hdata;
      hdata.h=ConfigInfo.h; hdata.dp=ConfigInfo.dp;  hdata.b=ConfigInfo.b;
      hdata.rhop0=ConfigInfo.rhop0;  hdata.gamma=ConfigInfo.gamma;
      hdata.massbound=ConfigInfo.massbound;
      hdata.massfluid=ConfigInfo.massfluid;
      hdata.np=Np; hdata.nfixed=Nfixed; hdata.nmoving=Nmoving; hdata.nfloat=Nfloat;
      hdata.nfluidout=NfluidOut;  hdata.timestep=GetPartTime();
      pf.write((char*)&hdata,sizeof(StHeadDatFullBi2));
      pf.flush();
      if(pf.fail())RunException(met,"File writing failure.",file);
      pf.close();
      OutCountSaved=NfluidOut;
    }
    else RunException(met,"Cannot open the file.",file);
  }
}

//==============================================================================
/// Returns information of a file BINX-011 without loading data.
//==============================================================================
JPartData::StInfoFileBin JPartData::GetInfoFromFileBin(const std::string &file)const{
  const char met[]="GetInfoFromFileBi2";
  StInfoFileBin inf;
  ifstream pf;
  pf.open(file.c_str(),ios::binary|ios::in);
  if(pf){
    inf=GetInfoFileBin(pf);
    if(inf.fmt!=11)RunException(met,"Invalid file format.",file);
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
  return(inf);
}

//==============================================================================
/// Returns information of a file BINX-20 without loading data.
//==============================================================================
JPartData::StInfoFileBi2 JPartData::GetInfoFromFileBi2(const std::string &file)const{
  const char met[]="GetInfoFromFileBi2";
  StInfoFileBi2 inf;
  ifstream pf;
  pf.open(file.c_str(),ios::binary|ios::in);
  if(pf){
    inf=GetInfoFileBi2(pf,false);
    if(inf.fmt!=20)RunException(met,"Invalid file format.",file);
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
  return(inf);
}

//==============================================================================
/// Loads file .inf with meta-info of the particles
//==============================================================================
void JPartData::LoadFileInfo(TpFmtFile fmt,const std::string &file){
  const char met[]="LoadFileInfo";
  Reset();
  JVarsAscii vars;
  vars.LoadFile(file);
  unsigned np,nbound,nfixed,nmoving,nfloat,nfluid;
  if(!vars.Exists("np"))RunException(met,"Missing variable \'np\'.",file);          else np=(unsigned)vars.ValueInt("np");
  if(!vars.Exists("nbound"))RunException(met,"Missing variable \'nbound\'.",file);  else nbound=(unsigned)vars.ValueInt("nbound");
  if(!vars.Exists("nfixed"))RunException(met,"Missing variable \'nfixed\'.",file);  else nfixed=(unsigned)vars.ValueInt("nfixed");
  nmoving=(unsigned)(vars.Exists("nmoving")? vars.ValueInt("nmoving"): 0);
  nfloat=(unsigned)(vars.Exists("nfloat")? vars.ValueInt("nfloat"): 0);
  if(!vars.Exists("nfluid"))RunException(met,"Missing variable \'nfluid\'.",file);  else nfluid=(unsigned)vars.ValueInt("nfluid");
  float dp,h,b,rhop0,gamma,massbound,massfluid;
  if(!vars.Exists("dp"))RunException(met,"Missing variable \'dp\'.",file);                else dp=vars.ValueFloat("dp");
  if(!vars.Exists("h"))RunException(met,"Missing variable \'h\'.",file);                  else h=vars.ValueFloat("h");
  if(!vars.Exists("b"))RunException(met,"Missing variable \'b\'.",file);                  else b=vars.ValueFloat("b");
  if(!vars.Exists("rhop0"))RunException(met,"Missing variable \'rhop0\'.",file);          else rhop0=vars.ValueFloat("rhop0");
  if(!vars.Exists("gamma"))RunException(met,"Missing variable \'gamma\'.",file);          else gamma=vars.ValueFloat("gamma");
  if(!vars.Exists("massbound"))RunException(met,"Missing variable \'massbound\'.",file);  else massbound=vars.ValueFloat("massbound");
  if(!vars.Exists("massfluid"))RunException(met,"Missing variable \'massfluid\'.",file);  else massfluid=vars.ValueFloat("massfluid");
  int data2d=0;
  if(vars.Exists("data2d"))data2d=vars.ValueInt("data2d");
  Config(fmt,np,nbound,nfluid,nfixed,nmoving,nfloat,dp,h,b,rhop0,gamma,massbound,massfluid,data2d==1);
  if(vars.Exists("timebypart"))SetTimeByPart(vars.ValueFloat("timebypart"));
}

//==============================================================================
/// Stores file .inf with meta-info of the particles
//==============================================================================
void JPartData::SaveFileInfo(const std::string &file)const{
  JVarsAscii vars;
  vars.AddValueInt("data2d",(Data2D? 1: 0));
  vars.AddValueInt("np",Np);
  vars.AddValueInt("nbound",Nbound);
  vars.AddValueInt("nfixed",Nfixed);
  vars.AddValueInt("nmoving",Nmoving);
  vars.AddValueInt("nfloat",Nfloat);
  vars.AddValueInt("nfluid",Nfluid);
  vars.AddValueDouble("dp",ConfigInfo.dp);
  vars.AddValueDouble("h",ConfigInfo.h,"%15.7E");
  vars.AddValueDouble("b",ConfigInfo.b,"%15.7E");
  vars.AddValueDouble("rhop0",ConfigInfo.rhop0);
  vars.AddValueDouble("gamma",ConfigInfo.gamma);
  vars.AddValueDouble("massbound",ConfigInfo.massbound,"%15.7E");
  vars.AddValueDouble("massfluid",ConfigInfo.massfluid,"%15.7E");
  vars.AddValueDouble("timebypart",GetTimeByPart());
  vars.SaveFile(file);
}

//==============================================================================
/// Stores file .dat with meta-info of the particles.
//==============================================================================
void JPartData::SaveFileInfoDat(const std::string &file)const{
  const char met[]="SaveFileInfoDat";
  ofstream pf;
  pf.open(file.c_str());
  if(pf){
    char cad[4096];
    sprintf(cad,"%9u %9u %9u %9u %9u %9u",Np,Nbound,Nfixed,Nmoving,Nfloat,Nfluid); pf << cad;
    sprintf(cad," %15.7E %15.7E %15.7E",ConfigInfo.dp,ConfigInfo.h,ConfigInfo.b); pf << cad;
    sprintf(cad," %15.7E %15.7E %15.7E %15.7E",ConfigInfo.rhop0,ConfigInfo.gamma,ConfigInfo.massbound,ConfigInfo.massfluid); pf << cad;
    sprintf(cad," %15.7E",GetTimeByPart()); pf << cad;
    pf << endl;
    if(pf.fail())RunException(met,"File writing failure.",file);
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Sorts paticle data in the arrays and returns the number of copied particles.
/// - size: Indicates the number of reserved positions in the arrays.
/// - full: Indicates that data of excluded particles must be included.
/// - id,vel,rhop: They can be NULL if they are not necessary.
//==============================================================================
unsigned JPartData::GetDataSort(unsigned size,unsigned *id,tfloat3 *pos,tfloat3 *vel,float *rhop,bool full)const{
  const char met[]="GetDataSort";
  unsigned rsize=(full? Np: Np-NfluidOut);
  if(rsize>size)RunException(met,"Target vector dimension is insufficient.");
  if(!DataBoundOk())RunException(met,"Boundary particle data is missing.");
  if(full&&NfluidOut>OutCount)RunException(met,"Excluded pasticle data is missing.");
  if(FmtData==FmtAscii){
    if(id)memcpy(id,Id,rsize*sizeof(unsigned));
    memcpy(pos,Pos,rsize*sizeof(tfloat3));
    if(vel)memcpy(vel,Vel,rsize*sizeof(tfloat3));
    if(rhop)memcpy(rhop,Rhop,rsize*sizeof(float));
  }
  else if(FmtData==FmtBin){
    unsigned npok=Np-NfluidOut;
    unsigned *idx=NULL;
    //-Creates indexes of adjacent positioning for no excluded fluid particles.
    if(NfluidOut){
      idx=new unsigned[Nfluid];
      unsigned cidx=Nbound;
      for(unsigned c=0;c<Nfluid;c++)if(OrderOut[c]<0){ idx[c]=cidx; cidx++; }
    }
    //-Places all particles not excluded.
    #ifdef DBG_JPartData
      unsigned nb=0,nf=0;
    #endif 
    for(unsigned c=0;c<npok;c++){
      unsigned idp=Id[c];
      if(idp>=Nbound&&idx)idp=idx[idp-Nbound];
      if(id)id[idp]=Id[c];
      pos[idp]=Pos[c]; 
      if(vel)vel[idp]=Vel[c]; 
      if(rhop)rhop[idp]=Rhop[c];
      #ifdef DBG_JPartData
        if(idp<Nbound)nb++; else nf++;
      #endif 
    }
    #ifdef DBG_JPartData
      if(nb!=Nbound||nf!=(Nfluid-NfluidOut))RunException(met,"Error placing no excluded partices.");
    #endif 
    delete[] idx;
  }
  else if(FmtData==FmtBi2){
    unsigned npok=Np-NfluidOut;
    if(id)memcpy(id,Id,npok*sizeof(unsigned));
    memcpy(pos,Pos,npok*sizeof(tfloat3));
    if(vel)memcpy(vel,Vel,npok*sizeof(tfloat3));
    if(rhop)memcpy(rhop,Rhop,npok*sizeof(float));
  }
  else RunException(met,"Unrecognised file format.");
  //-Copies data of excluded particles.
  if(full&&(FmtData==FmtBin||FmtData==FmtBi2)){
    for(unsigned c=0;c<NfluidOut;c++){
      StParticleOut* out=Out+c;
      unsigned idp=Np-c-1;
      if(id)id[idp]=out->id; 
      pos[idp]=out->pos; 
      if(vel)vel[idp]=out->vel; 
      if(rhop)rhop[idp]=out->rhop;
    }
  }
  return(rsize);
}

//==============================================================================
/// Changes the allocated memory space for excluded particles.
//==============================================================================
void JPartData::SetOutSize(unsigned size){
  StParticleOut* out=new StParticleOut[size];
  OutCount=min(OutCount,size);
  if(OutCount)memcpy(out,Out,OutCount*sizeof(StParticleOut));
  delete[] Out;
  Out=out;
  OutSize=size;
}

//==============================================================================
/// Clears list of fluid excluded particles.
//==============================================================================
void JPartData::FluidOutClear(bool clearout){
  for(unsigned c=0;c<Nfluid;c++)OrderOut[c]=-1;
  NfluidOut=0;
  OutCountSaved=0;
  if(clearout)OutCount=0;
}


//==============================================================================
/// Compares data between 2 objects (debug).
//==============================================================================
void JPartData::DgCompare(const JPartData& pdat){
  const char met[]="DgCompare";
  if(Np!=pdat.GetNp())RunException(met,"Variable Np does not match.");
  if(Nbound!=pdat.GetNbound())RunException(met,"Variable Nbound does not match.");
  if(Nfluid!=pdat.GetNfluid())RunException(met,"Variable Nfluid does not match.");
  if(Nfixed!=pdat.GetNfixed())RunException(met,"Variable Nfixed does not match.");
  if(Nmoving!=pdat.GetNmoving())RunException(met,"Variable Nmoving does not match.");
  if(Nfloat!=pdat.GetNfloat())RunException(met,"Variable Nfloat does not match.");
  if(NfluidOut!=pdat.GetNfluidOut())RunException(met,"Variable NfluidOut does not match.");
  StConfig f1=GetConfigInfo();
  StConfig f2=pdat.GetConfigInfo();
  if(f1.dp!=f2.dp)RunException(met,"Variable dp does not match.");
  if(f1.h!=f2.h)RunException(met,"Variable h does not match.");
  if(f1.b!=f2.b)RunException(met,"Variable b does not match.");
  if(f1.rhop0!=f2.rhop0)RunException(met,"Variable rhop0 does not match.");
  if(f1.gamma!=f2.gamma)RunException(met,"Variable gamma does not match.");
  if(f1.massbound!=f2.massbound)RunException(met,"Variable massbound does not match.");
  if(f1.massfluid!=f2.massfluid)RunException(met,"Variable massfluid does not match.");
  unsigned *id1=new unsigned[Np],*id2=new unsigned[Np];
  tfloat3 *pos1=new tfloat3[Np],*pos2=new tfloat3[Np];
  tfloat3 *vel1=new tfloat3[Np],*vel2=new tfloat3[Np];
  float *rhop1=new float[Np],*rhop2=new float[Np];
  unsigned size1=GetDataSort(Np,id1,pos1,vel1,rhop1,true);
  unsigned size2=pdat.GetDataSort(Np,id2,pos2,vel2,rhop2,true);
  if(size1!=size2)RunException(met,"Variable size does not match.");
  for(unsigned c=0;c<size1;c++){
    if(id1[c]!=id2[c])RunException(met,"Variable id[] does not match.");
    if(pos1[c]!=pos2[c])RunException(met,"Variable pos[] does not match.");
    if(vel1[c]!=vel2[c])RunException(met,"Variable vel[] does not match.");
  }
  delete[] id1;
  delete[] pos1;
  delete[] vel1;
  delete[] id2;
  delete[] pos2;
  delete[] vel2;
}





