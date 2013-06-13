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

/// \file JDivideCpu.cpp \brief Implements the class \ref JDivideCpu.

#include "JDivideCpu.h"
#include "Functions.h"

#include <cstring>

//==============================================================================
/// Constructor of objects.
//==============================================================================
JDivideCpu::JDivideCpu(JLog *log,string dirout){
  ClassName="JDivideCpu";
  Log=log;
  DirOut=dirout;
  CellPart=NULL; Parts=NULL; PartsOut=NULL;
  PartsInCell=NULL; BeginBound=NULL; BeginFluid=NULL;
  VSort=NULL;
  Out=NULL; OutSize=0;
  ShowCellsInfo=true;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDivideCpu::~JDivideCpu(){
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDivideCpu::Reset(){
  delete[] CellPart; CellPart=NULL;
  delete[] Parts;    Parts=NULL;
  delete[] VSort;    VSort=NULL;
  delete[] Out;      Out=NULL; OutSize=0;
  ClearMemory();
  Np=0; Nbound=0; Nfluid=0; Nfixed=0; Nmoving=0; Nfloat=0; NfluidOut=0; Npb=0;
  NOutLast=0;
  Ndivb=0; Ndivf=0;
  ConfigRhopOut(false,0,0); RhopOutCount=0;
  CellInfoOk=false;
  OutClear();
}

//==============================================================================
/// Changes the allocated memory space for particles out.
//==============================================================================
void JDivideCpu::OutResize(unsigned size){
  StParticleOut* out=new StParticleOut[size];
  OutCount=min(OutCount,size);
  if(OutCount)memcpy(out,Out,OutCount*sizeof(StParticleOut));
  delete[] Out;
  Out=out;
  OutSize=size;
}

//==============================================================================
/// Calculates cells according to the dimensions of the case.
//==============================================================================
void JDivideCpu::CalcCells(){
  Ncells.x=int(ceil(DifMax.x/SizeCell));
  Ncells.y=int(ceil(DifMax.y/SizeCell));
  float modz=SizeCell*0.05f;       //-Margin of errors due to precision when using real values.
  Ncells.z=int(ceil((DifMax.z+modz)/SizeCell));
  Nsheet=Ncells.x*Ncells.y; Nct=Nsheet*Ncells.z; NctTot=Nct+1; NctTotMax=NctTot;
  LastMaxCellZ=0;
}

//==============================================================================
/// Shows the current limits of the simulation.
//==============================================================================
void JDivideCpu::VisuLimits()const{
  char tx[256];
  sprintf(tx,"Domain limits: (%.4f,%.4f,%.4f)-(%.4f,%.4f,%.4f)",PosMin.x,PosMin.y,PosMin.z,PosMin.x+DifMax.x,PosMin.y+DifMax.y,PosMin.z+DifMax.z);
  Log->Print(tx);
  sprintf(tx,"Dimensions: %f x %f x %f",DifMax.x,DifMax.y,DifMax.z);
  Log->Print(tx);
  if(ShowCellsInfo){
    sprintf(tx,"Cells of the domain: %d x %d x %d",Ncells.x,Ncells.y,Ncells.z);
    Log->Print(tx);
    Log->Print(fun::VarStr("Nsheet",Nsheet));
    Log->Print(fun::VarStr("Nct",Nct));
  }
}

//==============================================================================
// Modifies the limits of the domain.
//==============================================================================
void JDivideCpu::SetLimit(tfloat3 posmin,tfloat3 posmax){
  PosMin=posmin;
  DifMax.x=posmax.x-posmin.x; DifMax.y=posmax.y-posmin.y; DifMax.z=posmax.z-posmin.z;
  Log->Print("**Limits adjustment");
  CalcCells();
  VisuLimits();
}

//==============================================================================
/// Releases allocated memory during the Neighbour List step.
//==============================================================================
void JDivideCpu::ClearMemory(){
  delete[] PartsInCell; PartsInCell=NULL;
  delete[] BeginBound;  BeginBound=NULL;
  delete[] BeginFluid;  BeginFluid=NULL;
}

//==============================================================================
/// Allocates the memory necessary for each configuration.
//==============================================================================
void JDivideCpu::ConfigMemory(){
  const char met[]="ConfigMemory";
  ClearMemory();
  try{
    PartsInCell=new unsigned[NctTotMax];
    BeginBound=new unsigned[NctTotMax+1];
    BeginFluid=new unsigned[NctTotMax+1];
  }
  catch(const std::bad_alloc){
    RunException(met,"Could not allocate the requested memory.");
  }
  memset(BeginBound,0,sizeof(unsigned)*(NctTotMax+1));
  memset(BeginFluid,0,sizeof(unsigned)*(NctTotMax+1));
}

//==============================================================================
/// Returns the dynamic allocated memory.
//==============================================================================
unsigned int JDivideCpu::GetMemoryAlloc()const{
  unsigned int mem=0;
  mem+=sizeof(unsigned)*(NctTotMax + NctTotMax+1 + NctTotMax+1); //- PartsInCell+BeginBound+BeginFluid
  mem+=sizeof(unsigned)*(Np + Np + (Np-Npb));              //- CellPart+Parts+PartsOut
  mem+=(unsigned)max(sizeof(tfloat3)*Np,sizeof(tsymatrix3f)*(Np-Npb));   //- VSort
  mem+=sizeof(StParticleOut)*OutSize;                      //- Out
  return(mem);
}

//==============================================================================
/// Obtains limits of the domain.
//==============================================================================
void JDivideCpu::CheckLimits(unsigned pini,unsigned n,const tfloat3 *pos,tfloat3 &pmin,tfloat3 &pmax)const{
  if(n){ pmin.x=pos[pini].x; pmin.y=pos[pini].y; pmin.z=pos[pini].z; }
  else{ pmin.x=0; pmin.y=0; pmin.z=0; }
  pmax=pmin;
  const unsigned pfin=pini+n;
  for(unsigned p=pini;p<pfin;p++){
    const tfloat3 *ps=(pos+p);
    if(pmin.x>ps->x)pmin.x=ps->x;
    if(pmin.y>ps->y)pmin.y=ps->y;
    if(pmin.z>ps->z)pmin.z=ps->z;
    if(pmax.x<ps->x)pmax.x=ps->x;
    if(pmax.y<ps->y)pmax.y=ps->y;
    if(pmax.z<ps->z)pmax.z=ps->z;
  }
}

//==============================================================================
/// Configuration of the Neighbour List: Memory allocation, domain and cells.
//==============================================================================
void JDivideCpu::Config(float sizecell,float incz,unsigned np,unsigned nfixed,unsigned nmoving,unsigned nfloat,unsigned nfluid,const tfloat3 *pos,float dp){
  const char met[]="Config";
  Reset();
  string aux="**Initialising ";
  Log->Print(aux+ClassName);
  if(np<1)RunException(met,"Invalid number of particles.");
  if(sizecell<=0)RunException(met,"Invalid value for cell size.");
  Np=np; Nfixed=nfixed; Nmoving=nmoving; Nfloat=nfloat; Nfluid=nfluid;
  Nbound=Nfixed+Nmoving+Nfloat;
  if(Np!=Nbound+Nfluid)RunException(met,"Error in number of particles.");
  Npb=Nbound-Nfloat;
  SizeCell=sizecell; OvSizeCell=float(1.0/SizeCell); IncZ=incz;
  //-Allocates memory dependent on Np.
  try{
    CellPart=new unsigned[Np]; 
    Parts=new unsigned[Np];
    VSort=(byte*)new byte[max(sizeof(tfloat3)*Np,sizeof(tsymatrix3f)*(Np-Npb))];
  }
  catch(const std::bad_alloc){
    RunException(met,"Could not allocate the requested memory.");
  }
  VSortInt=(int*)VSort;
  VSortFloat=(float*)VSort;
  VSortFloat3=(tfloat3*)VSort;
  VSortFloat6=(tsymatrix3f*)VSort;
  //-Calculates limits of the domain.
  tfloat3 pmin,pmax;
  CheckLimits(0,Np,pos,pmin,pmax);
  float mod=sizecell*0.05f;  
  pmin.x-=mod; pmin.y-=mod; pmin.z-=mod;
  pmax.x+=mod; pmax.y+=mod; pmax.z+=mod;
  DifMax.x=pmax.x-pmin.x; DifMax.y=pmax.y-pmin.y; DifMax.z=pmax.z-pmin.z;
  CalcCells();
  if(ShowCellsInfo){
    char tx[128];
    sprintf(tx,"Cells of the initial domain: %d x %d x %d  (%d)",Ncells.x,Ncells.y,Ncells.z,(Ncells.x*Ncells.y*Ncells.z));
    Log->Print(tx);
  }
  pmax.z=(pmax.z-pmin.z)*(IncZ+1.f)+pmin.z;
  PosMin=pmin;
  DifMax.x=pmax.x-pmin.x; DifMax.y=pmax.y-pmin.y; DifMax.z=pmax.z-pmin.z;
  CalcCells();
  VisuLimits();
  CellInfoOk=false;
}

//==============================================================================
/// Computes the Neighbour List of boundary particles (Npb).
//==============================================================================
void JDivideCpu::DivBoundary(const tfloat3 *pos){
  const char met[]="DivBoundary";
  BoundMaxCellZ=0;
  //-Initialise amount and position of particles per cell.
  memset(PartsInCell,0,sizeof(unsigned)*NctTotMax);
  //-Calculates the cell of each particle.
  for(unsigned p=0;p<Npb;p++){
    float dx=pos[p].x-PosMin.x,dy=pos[p].y-PosMin.y,dz=pos[p].z-PosMin.z;
    if(dx<0||dy<0||dz<0||dx>=DifMax.x||dy>=DifMax.y||dz>=DifMax.z)RunException(met,"A boundary particle was found outside the domain.");
    unsigned cx=unsigned(dx*OvSizeCell),cy=unsigned(dy*OvSizeCell),cz=unsigned(dz*OvSizeCell);
    if(BoundMaxCellZ<cz)BoundMaxCellZ=cz;
    int box=cx+cy*Ncells.x+cz*Nsheet;
    CellPart[p]=box;
    PartsInCell[box]++;//-Number of particles of each cell.
  }
  //-Calculates the first boundary particle of each cell. 
  BeginBound[0]=0;
  for(unsigned box=0;box<NctTotMax;box++){
    BeginBound[box+1]=BeginBound[box]+PartsInCell[box];
    PartsInCell[box]=0;
  }
  //-Addresses particles in their cells.  
  for(unsigned p=0;p<Npb;p++){
    unsigned box=CellPart[p];
    Parts[BeginBound[box]+PartsInCell[box]]=p;
    PartsInCell[box]++;
  }
  CellInfoOk=false;
  Ndivb++;
}

//==============================================================================
/// Computes the Neighbour List of fluid particles.
//==============================================================================
bool JDivideCpu::DivFluid(const unsigned* idp,const tfloat3* pos,const tfloat3* vel,const float* rhop,float time){
  const char met[]="DivFluid";
  bool modifdata=false;
  //-Checks limits and calculates the cell of each particle.
  FluidMaxCellZ=BoundMaxCellZ;
  unsigned outcountpre=OutCount;
  const unsigned pfin=Np-NfluidOut;
  for(unsigned p=Npb;p<pfin;p++){
    float dx=pos[p].x-PosMin.x,dy=pos[p].y-PosMin.y,dz=pos[p].z-PosMin.z;
    bool partfloating=(Nfloat&&idp[p]<Nbound? true: false);                                      //-Particle floating. Floating particles are in the same block as fluid ones.
    bool partin=(dx>=0&&dy>=0&&dz>=0&&dx<DifMax.x&&dy<DifMax.y&&dz<DifMax.z);                    //-Particle within the domain.
    bool partrhopout=(partin&&!partfloating&&RhopOut&&(rhop[p]<RhopOutMin||rhop[p]>RhopOutMax)); //-Particle excluded by \ref RhopOut.

    if(partin&&!partrhopout){ //-Within the domain and valid density.
      unsigned cx=unsigned(dx*OvSizeCell),cy=unsigned(dy*OvSizeCell),cz=unsigned(dz*OvSizeCell);
      if(FluidMaxCellZ<cz)FluidMaxCellZ=cz;
      CellPart[p]=cx+cy*Ncells.x+cz*Nsheet; 
    }
    else{ //-Out of the domain and no-valid density.
      if(partfloating){
        char tx[1024];
        sprintf(tx,"particle floating out> p:%u id:%u pos:(%f,%f,%f)\n",p,idp[p],pos[p].x,pos[p].y,pos[p].z);
        Log->Print(tx);
        RunException(met,"A floating body particle was found outside the domain.");
      }
      if(OutCount>=int(OutSize))OutResize(min(max(unsigned(1000),unsigned(OutSize*2)),Nfluid));
      StParticleOut* out=Out+OutCount;
      out->p=p;
      out->id=idp[p];
      out->pos=pos[p];
      out->vel=vel[p];
      out->rhop=rhop[p];
      out->timeout=time;
      OutCount++;
      if(partrhopout)RhopOutCount++;
    }
  }
  //-Adjust doamin in Z axis.
  if(LastMaxCellZ!=FluidMaxCellZ){
    Ncells.z=FluidMaxCellZ+1;
    Nct=Nsheet*Ncells.z;
    NctTot=Nct+1;
    LastMaxCellZ=FluidMaxCellZ;
    modifdata=true;
  }
  //-Assigns cell Nct to the new Particles-Out.
  if(outcountpre<OutCount){
    NOutLast=OutCount-outcountpre;
    for(unsigned c=outcountpre;c<OutCount;c++)CellPart[Out[c].p]=Nct;
    NfluidOut+=(OutCount-outcountpre);
    modifdata=true;
  }
  else NOutLast=0;
  //-Counts particles per cell.
  memset(PartsInCell,0,sizeof(int)*NctTot);
  for(unsigned p=Npb;p<pfin;p++)PartsInCell[CellPart[p]]++;
  //-Calculates the first fluid particle of each cell. 
  BeginFluid[0]=Npb;
  for(unsigned box=0;box<NctTot;box++){
    BeginFluid[box+1]=BeginFluid[box]+PartsInCell[box];
    PartsInCell[box]=0;
  }
  //-Addresses particles in their cells.  
  for(unsigned p=Npb;p<pfin;p++){
    unsigned box=CellPart[p];
    Parts[BeginFluid[box]+PartsInCell[box]]=p;
    PartsInCell[box]++;
  }
  CellInfoOk=false;
  Ndivf++;
  return(modifdata);
}

//==============================================================================
/// Stores the information of the excluded particles in the end of the arrays. 
//==============================================================================
void JDivideCpu::GetDataOut(unsigned* id,tfloat3* pos,tfloat3* vel,float* rhop,bool clearout){
  StParticleOut* out=Out;
  const unsigned pfin=Np-NfluidOut+OutCount;
  for(unsigned p=Np-NfluidOut;p<pfin;p++){
    id[p]=out->id;
    pos[p]=out->pos;
    vel[p]=out->vel;
    rhop[p]=out->rhop;
    out++;
  }
  if(clearout)OutClear();
}

//==============================================================================
/// Returns the info of the Neighbour List starting from the particles position.
//==============================================================================
JDivideCpu::StDivideInfo JDivideCpu::CalcDivideInfo(float scell,float border,float incz,unsigned np,const tfloat3* pos){
  StDivideInfo info;
  memset(&info,0,sizeof(StDivideInfo));
  info.scell=scell;
  info.border=border;
  info.incz=incz;
  if(np){
    //-Calculates minimum and maximum position. 
    tfloat3 pmin,pmax;
    pmin=pos[0]; pmax=pmin;
    for(unsigned p=1;p<np;p++){
      const tfloat3 *ps=(pos+p);
      if(pmin.x>ps->x)pmin.x=ps->x;
      if(pmin.y>ps->y)pmin.y=ps->y;
      if(pmin.z>ps->z)pmin.z=ps->z;
      if(pmax.x<ps->x)pmax.x=ps->x;
      if(pmax.y<ps->y)pmax.y=ps->y;
      if(pmax.z<ps->z)pmax.z=ps->z;
    }
    pmin.x-=border; pmin.y-=border; pmin.z-=border;
    pmax.x+=border; pmax.y+=border; pmax.z+=border;
    info.posmin=pmin;
    info.posmax=pmax;
    info.difmax=TFloat3(pmax.x-pmin.x,pmax.y-pmin.y,pmax.z-pmin.z);
    info.ncellx=unsigned(ceil(info.difmax.x/scell));
    info.ncelly=unsigned(ceil(info.difmax.y/scell));
    info.ncellz=unsigned(ceil(info.difmax.z/scell));
    info.ncells=info.ncellx*info.ncelly*info.ncellz;
    info.poszmax_incz=info.difmax.z*(incz+1.f)+info.posmin.z;
    info.difzmax_incz=info.poszmax_incz-info.posmin.z;
    info.ncellz_incz=unsigned(ceil(info.difzmax_incz/scell));
    info.ncells_incz=info.ncellx*info.ncelly*info.ncellz_incz;
  }
  return(info);
}

//==============================================================================
/// Reorders data of boundary particles.
//==============================================================================
void JDivideCpu::SortBoundary(unsigned *vec){
  for(unsigned p=0;p<Npb;p++)VSortInt[p]=vec[Parts[p]];
  memcpy(vec,VSortInt,sizeof(unsigned)*Npb);
}
//==============================================================================
void JDivideCpu::SortBoundary(float *vec){
  for(unsigned p=0;p<Npb;p++)VSortFloat[p]=vec[Parts[p]];
  memcpy(vec,VSortFloat,sizeof(float)*Npb);
}
//==============================================================================
void JDivideCpu::SortBoundary(tfloat3 *vec){
  for(unsigned p=0;p<Npb;p++)VSortFloat3[p]=vec[Parts[p]];
  memcpy(vec,VSortFloat3,sizeof(tfloat3)*Npb);
}

//==============================================================================
/// Reorder data of fluid particles.
//==============================================================================
void JDivideCpu::SortFluid(unsigned *vec){
  const unsigned noutpre=NfluidOut-NOutLast;
  const unsigned pfin=Np-noutpre;
  for(unsigned p=Npb;p<pfin;p++)VSortInt[p]=vec[Parts[p]];
  memcpy(vec+Npb,VSortInt+Npb,sizeof(unsigned)*(Np-Npb-noutpre));
}
//==============================================================================
void JDivideCpu::SortFluid(float *vec){
  const unsigned noutpre=NfluidOut-NOutLast;
  const unsigned pfin=Np-noutpre;
  for(unsigned p=Npb;p<pfin;p++)VSortFloat[p]=vec[Parts[p]];
  memcpy(vec+Npb,VSortFloat+Npb,sizeof(float)*(Np-Npb-noutpre));
}
//==============================================================================
void JDivideCpu::SortFluid(tfloat3 *vec){
  const unsigned noutpre=NfluidOut-NOutLast;
  const unsigned pfin=Np-noutpre;
  for(unsigned p=Npb;p<pfin;p++)VSortFloat3[p]=vec[Parts[p]];
  memcpy(vec+Npb,VSortFloat3+Npb,sizeof(tfloat3)*(Np-Npb-noutpre));
}
//==============================================================================
void JDivideCpu::SortOnlyFluid(tsymatrix3f *vec){     
  const unsigned noutpre=NfluidOut-NOutLast;
  const unsigned pfin=Np-noutpre;
  for(unsigned p=Npb;p<pfin;p++)VSortFloat6[p-Npb]=vec[Parts[p]-Npb];
  memcpy(vec,VSortFloat6,sizeof(tsymatrix3f)*(Np-Npb-noutpre));
}

//==============================================================================
/// Indicates whether the cell is empty or not.
//==============================================================================
bool JDivideCpu::CellNoEmpty(int box,byte kind)const{
  return(kind==1? BeginBound[box]<BeginBound[box+1]: BeginFluid[box]<BeginFluid[box+1]);
}

//==============================================================================
/// Returns the first particle of a cell.
//==============================================================================
int JDivideCpu::CellBegin(int box,byte kind)const{
  return(kind==1? BeginBound[box]: BeginFluid[box]);
}

//==============================================================================
/// Returns the number of particles of a cell.
//==============================================================================
int JDivideCpu::CellSize(int box,byte kind)const{
  return(kind==1? int(BeginBound[box+1]-BeginBound[box]): int(BeginFluid[box+1]-BeginFluid[box]));
}


//==============================================================================
/// Collects information about the cells.
//==============================================================================
JDivideCpu::StCellInfo JDivideCpu::GetCellInfo(){
  const char met[]="GetCellInfo";
  if(!CellInfoOk){
    CellInfo=GetCellInfoData(Np,Nct,BeginBound,BeginFluid);
    CellInfoOk=true;
  }
  return(CellInfo);
}

//==============================================================================
/// Collects information about the cells of a given list.
//==============================================================================
JDivideCpu::StCellInfo JDivideCpu::GetCellInfoData(unsigned n,unsigned nct,unsigned *beginb,unsigned *beginf)const{
  const char met[]="GetCellInfoData";
  if(!nct)RunException(met,"There are no cells for analysis.");
  if(beginf[nct]==0)RunException(met,"The data do not provide a full neighbour list.");
  unsigned vmin=-1,nmin=0,vmax=0,nmax=0,suma=0,nsuma=0,nvoid=0;
  for(unsigned box=0;box<nct;box++){
    unsigned count=(beginb[box+1]-beginb[box])+(beginf[box+1]-beginf[box]);
    if(count>0){ 
      if(vmin>count||!nmin){ vmin=count; nmin=1; }
      else if(vmin==count)nmin++;
      if(vmax<count||!nmax){ vmax=count; nmax=1; }
      else if(vmax==count)nmax++;
      suma+=count; nsuma++; 
    }
    else nvoid++;
  }
  StCellInfo celinfo;
  celinfo.cells_void=nvoid;
  celinfo.min_incell=vmin;
  celinfo.min_count=nmin;
  celinfo.max_incell=vmax;
  celinfo.max_count=nmax;
  celinfo.media=float(suma)/nsuma;
  return(celinfo);
}








