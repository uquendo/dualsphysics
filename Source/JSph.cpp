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

/// \file JSph.cpp \brief Implements the class \ref JSph

#include "JSph.h"
#include "Functions.h"
#include "JSphMotion.h"
#include "JXml.h"
#include "JSpaceCtes.h"
#include "JSpaceEParms.h"
#include "JSpaceParts.h"
#include "JFormatFiles2.h"
#include "JProbe.h"

//==============================================================================
/// Constructor.
//==============================================================================
JSph::JSph(){
  ClassName="JSph";
  TStep=STEP_None;
  Idp=NULL; Ridp=NULL; Pos=NULL; Vel=NULL; Rhop=NULL; 
  DtPf=NULL;
  FtObjs=NULL;
  Motion=new JSphMotion();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSph::~JSph(){
  Reset();
  delete Motion;
}

//==============================================================================
/// Returns the code version in text format.
//==============================================================================
std::string JSph::GetVersionStr(){
  char cad[128];
  sprintf(cad,"%1.2f",float(VersionMajor)/100);
  return(cad);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSph::Reset(){
  Simulate2D=false;
  Hdiv=0;
  OmpThreads=1;
  Log=NULL;
  CaseLoaded=false;
  VerletSteps=40;
  ShepardSteps=30;
  DBCSteps=1;
  memset(&CubicCte,0,sizeof(StCubicCte));
  memset(&WendlandCte,0,sizeof(StWendlandCte));
  CaseName="";
  RunName="";
  DirCase="";
  DirOut="";
  PartBeginDir=""; PartBegin=0; PartBeginFirst=0;
  Np=0; Nbound=0; Nfixed=0; Nmoving=0; Nfloat=0; Nfluid=0; H=0; NpOk=0;
  Npb=0;
  Part=0; Ndiv=0;
  IncZ=0.f;
  DtModif=0;
  ClearRidp();
  AllocMemory(0);
  Pdata.Reset();
  SvData=byte(SDAT_Binx2);
  SvSteps=false;
  SvDt=false;
  SvRes=false;
  SvTimers=false;
  BoundDatVer=0;
  if(DtPf){
    if(DtPf->is_open())DtPf->close();
    delete DtPf; DtPf=NULL;
  }
  Motion->Reset();
  MotionTimeMod=0;
  MotionObjCount=0;
  RhopOut=false; RhopOutMin=0; RhopOutMax=0; RhopOutCount=0;
  FtCount=0;
  AllocMemoryFloating(0);
}
//==============================================================================
/// Allocates memory of main data.
//==============================================================================
void JSph::AllocMemory(int np){
  MemCpuStatic=0;
  delete[] Idp;  Idp=NULL;
  delete[] Ridp; Ridp=NULL;
  delete[] Pos;  Pos=NULL;
  delete[] Vel;  Vel=NULL;
  delete[] Rhop; Rhop=NULL;
  if(np>0){
    try{
      Idp=new unsigned[np]; Ridp=new unsigned[np];  MemCpuStatic+=(sizeof(unsigned)*np)*2;
      Pos=new tfloat3[np];  Vel=new tfloat3[np];    MemCpuStatic+=(sizeof(tfloat3)*np)*2;
      Rhop=new float[np];                           MemCpuStatic+=(sizeof(float)*np);
    }
    catch(const std::bad_alloc){
      RunException("AllocMemory","Could not allocate the requested memory.");
    }
  }
}

//==============================================================================
/// Allocates memory of floating objectcs.
//==============================================================================
void JSph::AllocMemoryFloating(unsigned ftcount){
  delete[] FtObjs; FtObjs=NULL;
  if(ftcount)FtObjs=new StFloatingData[ftcount];  MemCpuStatic+=(sizeof(StFloatingData)*ftcount);
  //delete[] FtDist; FtDist=NULL;
  //if(ftnp)FtDist=new tfloat3[ftnp];                MemCpuStatic+=(sizeof(tfloat3)*ftnp);
}

//==============================================================================
/// Allocates memory for probes arrays needed to write flw.
//==============================================================================
void JSph::AllocMemoryProbes(unsigned long np){
  delete[] ProbePos;  ProbePos=NULL;
  delete[] ProbeVel;  ProbeVel=NULL;
  delete[] ProbeRhop; ProbeRhop=NULL;
  if(np>0){
    try{
      ProbePos=new tfloat3[np];  ProbeVel=new tfloat3[np];    MemCpuStatic+=(sizeof(tfloat3)*np)*2;
      ProbeRhop=new float[np];                           MemCpuStatic+=(sizeof(float)*np);
    }
    catch(const std::bad_alloc){
      RunException("AllocMemory","Could not allocate the requested memory.");
    }
  }
}

//==============================================================================
/// Loads the case to be executed.
//==============================================================================
void JSph::LoadCase(){
  const char* met="LoadCase";
  string filexml=DirCase+CaseName+".xml";
  if(!fun::FileExists(filexml))RunException(met,"Case configuration was not found.",filexml);
  JXml xml; xml.LoadFile(filexml);
  JSpaceCtes ctes;     ctes.LoadXmlRun(&xml,"case.execution.constants");
  JSpaceEParms eparms; eparms.LoadXml(&xml,"case.execution.parameters");
  JSpaceParts parts;    parts.LoadXml(&xml,"case.execution.particles");

  //-Execution parameters.
  switch(eparms.GetValueInt("StepAlgorithm",true,1)){
    case 1:  TStep=STEP_Verlet;     break;
    case 2:  TStep=STEP_Symplectic;  break;
    default: RunException(met,"Step algorithm is not valid.");
  }
  VerletSteps=eparms.GetValueInt("VerletSteps",true,40);
  switch(eparms.GetValueInt("Kernel",true,1)){
    case 1:  TKernel=KERNEL_Cubic;     break;
    case 2:  TKernel=KERNEL_Wendland;  break;
    default: RunException(met,"Kernel choice is not valid.");
  }
  ShepardSteps=eparms.GetValueInt("ShepardSteps",true,0);
  switch(eparms.GetValueInt("ViscoTreatment",true,1)){
    case 1:  TVisco=VISCO_Artificial;  break;
    case 2:  TVisco=VISCO_LaminarSPS;     break;
    default: RunException(met,"Viscosity treatment is not valid.");
  }
  Visco=eparms.GetValueFloat("Visco");
  Kgc=(eparms.GetValueInt("KernelGradientCorr",true,0)!=0);            
  DBCSteps=eparms.GetValueInt("DBCSteps",true,1);
  TimeMax=eparms.GetValueFloat("TimeMax");
  TimePart=eparms.GetValueFloat("TimeOut");
  Visco=eparms.GetValueFloat("Visco");
  DtIni=eparms.GetValueFloat("DtIni");
  DtMin=eparms.GetValueFloat("DtMin",true,0.00001f);
  IncZ=eparms.GetValueFloat("IncZ");

  //-Predefined constantes.
  H=ctes.GetH();
  CteB=ctes.GetB();
  Gamma=ctes.GetGamma();
  //Rhop0=ctes.GetRhop0();
  Eps=ctes.GetEps();
  CFLnumber=ctes.GetCFLnumber();
  Dp=ctes.GetDp();
  Gravity=ctes.GetGravity();
  MassFluid=ctes.GetMassFluid();
  MassBound=ctes.GetMassBound();

  //-Particle data.
  Np=parts.Count();
  Nfixed=parts.Count(PT_Fixed);
  Nmoving=parts.Count(PT_Moving);
  Nfloat=parts.Count(PT_Floating);
  Nfluid=parts.Count(PT_Fluid);
  Nbound=Np-Nfluid;
  Npb=Nbound-Nfloat;
  PartOutMax=unsigned(float(Nfluid)*eparms.GetValueFloat("PartsOutMax",true,1));

  //-Loads and configures MOTION.
  MotionObjCount=0;
  for(unsigned c=0;c<parts.CountBlocks();c++){
    const JSpacePartBlock &block=parts.GetBlock(c);
    if(block.Type==PT_Moving){
      if(MotionObjCount>=255)RunException(met,"The number of mobile objects exceeds the maximum.");
      //printf("block[%2d]=%d -> %d\n",c,block.GetBegin(),block.GetCount());
      MotionObjBegin[MotionObjCount]=block.GetBegin();
      MotionObjBegin[MotionObjCount+1]=MotionObjBegin[MotionObjCount]+block.GetCount();
      MotionObjCount++;
    }
  }
  if(int(MotionObjCount)<Motion->Init(&xml,"case.execution.motion",DirCase))RunException(met,"The number of mobile objects is lower than expected.");

  //-Loads floating objects.
  FtCount=parts.CountBlocks(PT_Floating);
  if(FtCount){
    AllocMemoryFloating(FtCount);
    unsigned cobj=0;
    for(unsigned c=0;c<parts.CountBlocks()&&cobj<FtCount;c++){
      const JSpacePartBlock &block=parts.GetBlock(c);
      if(block.Type==PT_Floating){
        const JSpacePartBlock_Floating &fblock=(const JSpacePartBlock_Floating &)block;
        StFloatingData* fobj=FtObjs+cobj;
        fobj->begin=fblock.GetBegin();
        fobj->count=fblock.GetCount();
        fobj->mass=fblock.GetMassbody();
        fobj->center=fblock.GetCenter();
        fobj->inertia=fblock.GetInertia();
        fobj->fvel=fblock.GetVelini();
        fobj->fomega=fblock.GetOmegaini();
        cobj++;
      }
    }
  }
  Log->Print("**Case configuration is loaded");

  //-Loads initial state of particles.
  AllocMemory(Np);
  Log->Print("Loading initial state of particles...");
  string filepart;
  string filepart1=DirCase+CaseName+".bi2";
  string filepart2=DirCase+CaseName+".bin";
  if(!fun::FileExists(filepart1)&&!fun::FileExists(filepart2))RunException(met,"Initial state of the case was not found.",filepart1);
  JPartData pdini;
  if(fun::FileExists(filepart1)){ pdini.LoadFileBi2(0,filepart1); filepart=filepart1; }
  else{ pdini.LoadFileBin(0,filepart2); filepart=filepart2; }
  JPartData::StConfig cf=pdini.GetConfigInfo();
  if(cf.h==0)RunException(met,"File data invalid",filepart);
  if(pdini.GetNp()!=Np||pdini.GetNbound()!=Nbound||pdini.GetNfixed()!=Nfixed||pdini.GetNmoving()!=Nmoving||pdini.GetNfloat()!=Nfloat)RunException(met,"Data file does not match the configuration of the case.",filepart);
  pdini.GetDataSort(Np,Idp,Pos,Vel,Rhop,true);
  Simulate2D=pdini.GetData2D();
  if(pdini.GetNfluidOut())RunException(met,"Particles OUT are not allowed in the initial case file.");
  NpOk=Np;

  //-Prepares Pdata for the maintenance of the particles.
  Pdata.Config(JPartData::FmtBi2,Np,Nbound,Nfluid,Nfixed,Nmoving,Nfloat,cf.dp,cf.h,cf.b,RHOPZERO,cf.gamma,cf.massbound,cf.massfluid,Simulate2D);

  Log->Print("**Initial state of particles is loaded");

  //-Loads initial state of probe particles if flw output requested
  if(SvData&SDAT_Flw){
        Log->Print("**Loading initial state of probe particles");
        string fileg3d=DirCase+CaseName+".g3d";
        if(!fun::FileExists(fileg3d)) RunException(met,"G3D file was not found.",fileg3d);
        JProbe probeini(0,ProbePos,ProbeVel,ProbeRhop);
        probeini.LoadCastNodesCountG3D(fileg3d);
        AllocMemoryProbes(probeini.GetNProbe());
        probeini.LoadFileG3D(fileg3d);
        //TODO: check and/or correct position mapping
  }

  CaseLoaded=true;
  BoundDatVer++;
}

//==============================================================================
/// Loads PART to restart simulation.
//==============================================================================
void JSph::LoadPartBegin(){
  const char* met="LoadPartBegin";
  if(!PartBeginDir.empty()){
    if(Nfloat)RunException(met,"Still can not use a startup file with floating body ...");  
    JPartData pdat;
    //-Reads file.
    pdat.LoadFile(JPartData::FmtBi2Out,0,PartBeginDir);
    pdat.LoadFile(JPartData::FmtBi2,0,PartBeginDir);
    pdat.LoadFile(JPartData::FmtBi2,PartBegin,PartBeginDir);
    //-Loads data from the boot file.
    if(pdat.GetNp()!=Np||pdat.GetNbound()!=Nbound||pdat.GetNfixed()!=Nfixed||pdat.GetNmoving()!=Nmoving||pdat.GetNfloat()!=Nfloat){
      string filepart=pdat.GetFileName(JPartData::FmtBi2,PartBegin,PartBeginDir);
      RunException(met,"Boot file data do not match the configuration of the case.",filepart);   
    }
    pdat.GetDataSort(Np,Idp,Pos,Vel,Rhop,true);
    //-Adjust paramaters to start.
    PartIni=PartBeginFirst;
    TimeStepIni=(!PartIni? 0: pdat.GetPartTime());
    //-Adjust motion for the instant of the loaded PART.
    if(Nmoving){
      MotionTimeMod=(!PartIni? pdat.GetPartTime(): 0);
      Motion->ProcesTime(0,TimeStepIni+MotionTimeMod);
    }
  }
  else{
    PartIni=0;
    TimeStepIni=0;
  }
}

//==============================================================================
/// Genarates output data files.
//==============================================================================
void JSph::SaveData(){
  const char met[]="SaveData";
  char cad[64],cad2[256];
  if(SvSteps){
    if(Part)sprintf(cad,"_%08d",Ndiv);  
    else sprintf(cad,"_--------");
    Pdata.SetMaskFileName(JPartData::FmtAscii,"PART_%08d");
    Pdata.SetMaskFileName(JPartData::FmtBin,"PartBinx_%08d.bin");
    Pdata.SetMaskFileName(JPartData::FmtBi2,"Part%08d.bi2");
  }
  else sprintf(cad,"_%04d",Part);   
  //-Ouput format.
  if(SvData&SDAT_Sphysics)Pdata.SaveFile(JPartData::FmtAscii,DirOut);
  if(SvData&SDAT_Binx2){
    Pdata.SaveFile(JPartData::FmtBi2,DirOut);
    Pdata.SaveFile(JPartData::FmtBi2Out,DirOut,Part!=PartIni);
  } 
  if((SvData&SDAT_Csv)||(SvData&SDAT_Vtk)){
    tfloat3 *pos,*vel;
    float *rhop;
    unsigned *id;
    Pdata.GetDataPointers(id,pos,vel,rhop);
    Pdata.SortDataOut();   
    if(SvData&SDAT_Csv)JFormatFiles2::ParticlesToCsv(DirOut+"PartCsv"+cad+".csv",Pdata.GetNp(),Pdata.GetNfixed(),Pdata.GetNmoving(),Pdata.GetNfloat(),Pdata.GetNfluid()-Pdata.GetNfluidOut(),Pdata.GetNfluidOut(),Pdata.GetPartTime(),pos,vel,rhop,NULL,NULL,id,NULL,NULL,NULL,NULL);
    if(SvData&SDAT_Vtk)JFormatFiles2::ParticlesToVtk(DirOut+"PartVtk"+cad+".vtk",Pdata.GetNp(),pos,vel,rhop,NULL,NULL,id,NULL,NULL,NULL,NULL);
  }
  if(SvData&SDAT_Flw){
     //TODO: TimeStep?
     Pdata.SaveFile(JPartData::FmtFlw,DirOut);
  }
  //-Time computation.
  if(Part>PartIni||Ndiv){   
    TimerPart.Stop();
    float tpart=TimerPart.GetElapsedTimeF()/1000.f;
    float tseg=tpart/(TimeStep-TimeStepM1);
    TimerSim.Stop();
    float tcalc=TimerSim.GetElapsedTimeF()/1000.f;
    float tleft=(tcalc/(TimeStep-TimeStepIni))*(TimeMax-(TimeStep-TimeStepIni));
    if(SvSteps)sprintf(cad2,"Part%s  %12.6f  %9.2f  %14s",cad,TimeStep,tseg,fun::GetDateTimeAfter(int(tleft)).c_str());
    else sprintf(cad2,"Part%s  %12.6f  %12d  %7d  %9.2f  %14s",cad,TimeStep,(Ndiv+1),Ndiv-PartNdiv,tseg,fun::GetDateTimeAfter(int(tleft)).c_str());
    Log->Print(cad2);
  }
  else{
    sprintf(cad2,"Part%s        %u particles successfully stored",cad,Np);
    Log->Print(cad2);
  }
  if(Np-NpOk!=PartOut){
    sprintf(cad2,"Particles out: %d  (total: %d)",(Np-NpOk)-PartOut,(Np-NpOk));
    Log->Print(cad2);
    PartOut=(Np-NpOk);
  }
  if(SvDt)SaveDt();
}

//==============================================================================
/// Prints out the allocated memory.
//==============================================================================
void JSph::PrintMemoryAlloc(){
  char cad[128];
  Log->Print("");
  sprintf(cad,"Allocated memory in CPU: %u (%.2f Mb)",GetMemoryCpu(),float(GetMemoryCpu())/(1024*1024));
  Log->Print(cad);
  if(!Cpu){
    sprintf(cad,"Allocated memory in GPU: %u (%.2f Mb)",GetMemoryGpu(),float(GetMemoryGpu())/(1024*1024));
    Log->Print(cad);
  }
  Log->Print("");
}

//==============================================================================
/// Prints out headers of PARTs.
//==============================================================================
void JSph::PrintHeadPart(){
  if(SvSteps){
    Log->Print("PART-Step           PartTime      Time/Seg   Finish time   ");
    Log->Print("==================  ============  =========  ==============");
  }
  else{
    Log->Print("PART       PartTime      TotalSteps    Steps    Time/Seg   Finish time   ");
    Log->Print("=========  ============  ============  =======  =========  ==============");
  }
}

//==============================================================================
/// Adjusts information of InfoDt for new PART or execution beginning.
//==============================================================================
void JSph::InfoDtReset(){
  InfoDt.dt.count=0;  InfoDt.dt.mean=0;  InfoDt.dt.meancount=0;
  InfoDt.dt1.count=0; InfoDt.dt1.mean=0; InfoDt.dt1.meancount=0;
  InfoDt.dt2.count=0; InfoDt.dt2.mean=0; InfoDt.dt2.meancount=0;
  InfoDt.dtmin=-1; InfoDt.dtmax=-1;
  InfoDt.outcount=(Np-NpOk);
  InfoDt.rhopoutcount=RhopOutCount;
}

//==============================================================================
/// Calculates average of values of DT.
//==============================================================================
void JSph::InfoDtMean(StInfoDtMean &sdt){
  if(sdt.count){
    float v=0; 
    for(int c=0;c<sdt.count;c++)v+=sdt.vdt[c];
    sdt.mean=((sdt.mean*sdt.meancount)+v)/(sdt.meancount+sdt.count);
    sdt.meancount+=sdt.count;
    sdt.count=0;
  }
}

//==============================================================================
/// Adds a new value of DT to StInfoDtMean.
//==============================================================================
void JSph::InfoDtAdd(float dt,StInfoDtMean &sdt){
  if(sdt.count>=INFODT_MAX)InfoDtMean(sdt);
  sdt.vdt[sdt.count]=dt;
  sdt.count++;
}

//==============================================================================
/// Adds new values of dt1 and dt2 to InfoDt.
//==============================================================================
void JSph::InfoDtsAdd(float dt1,float dt2){
  InfoDtAdd(dt1,InfoDt.dt1);
  InfoDtAdd(dt2,InfoDt.dt2);
}

//==============================================================================
/// Adds new value of dtstep to InfoDt.
//==============================================================================
void JSph::InfoDtStepAdd(float dtstep){
  if(InfoDt.dt.count||InfoDt.dt.meancount){
    if(InfoDt.dtmin>dtstep)InfoDt.dtmin=dtstep;
    if(InfoDt.dtmax<dtstep)InfoDt.dtmax=dtstep;
  }
  else{
    InfoDt.dtmin=dtstep;
    InfoDt.dtmax=dtstep;
  }
  InfoDtAdd(dtstep,InfoDt.dt);
}

//==============================================================================
/// Stores data of DT.
//==============================================================================
void JSph::SaveDt(){
  const char* met="SaveDt";
  if(!DtPf){
    DtPf=new ofstream;
    string fname=DirOut+"dt.csv";
    DtPf->open(fname.c_str());
    if(!(*DtPf))RunException(met,"Cannot open the file.",fname);
    (*DtPf) << "Part;Steps;TimeStep;RunTime;DtStep;DtStepMin;DtStepMax;Dt1;Dt2;DtModif;PartOut;PartRhopOut" << endl;
  }
  if(DtPf&&Part){
    InfoDtMean(InfoDt.dt);
    InfoDtMean(InfoDt.dt1);
    InfoDtMean(InfoDt.dt2);
    char cad[256];
    TimerTot.Stop();
    float runtime=TimerTot.GetElapsedTimeF()/1000.f;
    sprintf(cad,"%d;%d;%G;%G;%G;%G;%G;%G;%G;",Part,Ndiv,TimeStep,runtime,InfoDt.dt.mean,InfoDt.dtmin,InfoDt.dtmax,InfoDt.dt1.mean,InfoDt.dt2.mean);
    (*DtPf) << cad;
    unsigned nout=(Np-NpOk)-InfoDt.outcount;
    unsigned nrhopout=RhopOutCount-InfoDt.rhopoutcount;
    sprintf(cad,"%d;%u;%u",DtModif,nout,nrhopout);
    (*DtPf) << cad << endl;
  }
  InfoDtReset();
}

//==============================================================================
/// Generates file Run.csv with summary execution
//==============================================================================
void JSph::SaveRes(float tsim,float tseg,float ttot,const string &headplus,const string &detplus){
  const char* met="SaveRes";
  string fname=DirOut+"Run.csv";
  ofstream pf;
  pf.open(fname.c_str());
  if(pf){
    string timershead,timersvars;
    if(SvTimers){
      for(unsigned c=0;c<TimerGetCount();c++){
        timershead=timershead+";"+TimerGetName(c);
        timersvars=timersvars+";"+fun::FloatStr(TimerGetValue(c)/1000);
      }
    }
    pf << "#RunName;Np;TSimul;TSeg;TTotal;MemCpu;MemGpu;Steps;PartFiles;PartsOut;Hw;StepAlgo;Kernel;Viscosity;ViscoValue;KGC;DBC;Shepard;TMax;Nbound;Nfixed;H;RhopOut;PartRhopOut;CellMode" << headplus << timershead << endl;
    pf << RunName << ';' << Np << ';';
    pf << tsim << ';' << tseg << ';' << ttot << ';';
    pf << GetMemoryCpu() << ';' << GetMemoryGpu() << ';';
    pf << Ndiv << ';' << Part << ';' << (Np-NpOk) << ';';
    pf << Hardware << ';' << GetStepName(TStep) << ';' << GetKernelName(TKernel) << ';' << GetViscoName(TVisco) << ';'<< Visco << ';';
    pf << (Kgc? "true": "false") << ';' << DBCSteps << ';' << ShepardSteps << ';' << TimeMax << ';' << Nbound << ';' << Nfixed << ';' << H << ';';
    char rhopcad[256];
    if(RhopOut)sprintf(rhopcad,"(%G-%G)",RhopOutMin,RhopOutMax); else sprintf(rhopcad,"None");
    pf << rhopcad << ';' << RhopOutCount << ';' << GetNameCellMode(CellMode) << detplus << timersvars << endl;
    if(pf.fail())RunException(met,"Failed writing to file.",fname);
    pf.close();
  }
  else RunException(met,"File could not be opened.",fname);
}

//==============================================================================
/// Starts simulation.
//==============================================================================
void JSph::Run(JCfgRun *cfg,JLog *log){     
  const char* met="Run";
  Reset();
  TimerTot.Start();
  DirOut=fun::GetDirWithSlash(cfg->DirOut);   
  CaseName=cfg->CaseName; 
  DirCase=fun::GetDirWithSlash(fun::GetDirParent(CaseName));
  CaseName=CaseName.substr(DirCase.length());
  if(!CaseName.length())RunException(met,"Name of the case for execution was not indicated.");
  RunName=(cfg->RunName.length()? cfg->RunName: CaseName);
  PartBeginDir=cfg->PartBeginDir; PartBegin=cfg->PartBegin; PartBeginFirst=cfg->PartBeginFirst;
  Log=log;

  //-Output format options.
  SvData=byte(SDAT_None);
  if(cfg->Sv_Ascii)SvData|=byte(SDAT_Sphysics);
  if(cfg->Sv_Binx2)SvData|=byte(SDAT_Binx2);
  if(cfg->Sv_Csv)SvData|=byte(SDAT_Csv);
  if(cfg->Sv_Vtk)SvData|=byte(SDAT_Vtk);
  if(cfg->Sv_Flw)SvData|=byte(SDAT_Flw);

  SvDt=cfg->SvDt;
  SvRes=cfg->SvRes;
  SvTimers=cfg->SvTimers;
  SvSteps=false;        //-If "true", stores all the steps.

  printf("\n");
  char cad[256];
  sprintf(cad,"[Initialising %s v%s  %s]",ClassName.c_str(),GetVersionStr().c_str(),fun::GetDateTime().c_str());
  Log->Print(cad);

#ifdef USE_OPENMP
  //-Determines number of threads per host using OpenMP.
  if(Cpu&&cfg->OmpThreads!=1&&cfg->OmpMode!=OMPM_Single){
    OmpThreads=cfg->OmpThreads;
    if(OmpThreads<=0)OmpThreads=omp_get_num_procs();
    if(OmpThreads>MAXTHREADS_OMP)OmpThreads=MAXTHREADS_OMP;
    omp_set_num_threads(OmpThreads);
    Log->Print(string("Threads per host for parallel execution: ")+fun::IntStr(omp_get_max_threads()));
  }
  else{
    OmpThreads=1;
    omp_set_num_threads(OmpThreads);
  }
  OmpMode=cfg->OmpMode;
#else
  OmpThreads=1;
  OmpMode=OMPM_Single;
#endif

  string tx=fun::VarStr("CaseName",CaseName);
  tx=tx+";\n"+fun::VarStr("DirCase",DirCase)+";\n"+fun::VarStr("RunName",RunName)+";\n"+fun::VarStr("DirOut",DirOut)+";";
  if(!PartBeginDir.empty()){
    Log->Print(fun::VarStr("PartBegin",PartBegin));
    Log->Print(fun::VarStr("PartBeginDir",PartBeginDir));
    Log->Print(fun::VarStr("PartBeginFirst",PartBeginFirst));
  }
  LoadCase();
  //-Applies configuration using command line.
  if(cfg->TStep)TStep=cfg->TStep;
  if(cfg->VerletSteps>=0)VerletSteps=cfg->VerletSteps;
  if(cfg->TKernel)TKernel=cfg->TKernel;
  if(cfg->Kgc)Kgc=(cfg->Kgc==1);
  if(cfg->TVisco){ TVisco=cfg->TVisco; Visco=cfg->Visco; }
  if(cfg->ShepardSteps>=0)ShepardSteps=cfg->ShepardSteps;
  if(cfg->DBCSteps>=0)DBCSteps=cfg->DBCSteps;
  if(cfg->Incz>=0)IncZ=cfg->Incz;     //-Allowed increase in Z+ ( ncz =:= ncz*(1+IncZ) ).
  if(cfg->TimeMax>0)TimeMax=cfg->TimeMax;
  if(cfg->TimePart>=0)TimePart=cfg->TimePart;
  CellMode=cfg->CellMode;

  //-Corrects invalid values.
  if(DBCSteps<1)DBCSteps=1;

  RhopOut=cfg->RhopOut; RhopOutMin=cfg->RhopOutMin; RhopOutMax=cfg->RhopOutMax;

  //-Constants for computation.
  Cs0=sqrt(int(Gamma)*CteB/RHOPZERO);
  Dosh=H*2; H2=H*H; Fourh2=H2*4; Eta2=(H*0.1f)*(H*0.1f);
  if(Simulate2D){
    if(TKernel==KERNEL_Cubic){
      CubicCte.a1=float(10.0f/(PI*7.f));
      CubicCte.a2=CubicCte.a1/H2;
      CubicCte.aa=CubicCte.a1/(H*H*H);
      CubicCte.a24=0.25f*CubicCte.a2;
      CubicCte.c1=-3.0f*CubicCte.aa;
      CubicCte.d1=9.0f*CubicCte.aa/4.0f;
      CubicCte.c2=-3.0f*CubicCte.aa/4.0f;
      float deltap=1.f/1.5f;
      float wdeltap=CubicCte.a2*(1.f-1.5f*deltap*deltap+0.75f*deltap*deltap*deltap);
      CubicCte.od_wdeltap=1.f/wdeltap;
      CteShepard=CubicCte.a2;
    }
    if(TKernel==KERNEL_Wendland){
      WendlandCte.awen=0.557f/(H*H);
      WendlandCte.bwen=-2.7852f/(H*H*H);
      CteShepard=WendlandCte.awen;
    }
  }
  else{
    if(TKernel==KERNEL_Cubic){
      CubicCte.a1=float(1.0f/PI);
      CubicCte.a2=CubicCte.a1/(H*H*H);
      CubicCte.aa=CubicCte.a1/(H*H*H*H);
      CubicCte.a24=0.25f*CubicCte.a2;
      CubicCte.c1=-3.0f*CubicCte.aa;
      CubicCte.d1=9.0f*CubicCte.aa/4.0f;
      CubicCte.c2=-3.0f*CubicCte.aa/4.0f;
      float deltap=1.f/1.5f;
      float wdeltap=CubicCte.a2*(1.f-1.5f*deltap*deltap+0.75f*deltap*deltap*deltap);
      CubicCte.od_wdeltap=1.f/wdeltap;
      CteShepard=CubicCte.a2;     
    }
    if(TKernel==KERNEL_Wendland){
      WendlandCte.awen=0.41778f/(H*H*H);
      WendlandCte.bwen=-2.08891f/(H*H*H*H);
      CteShepard=WendlandCte.awen;   
    }
  }
  
  if(TVisco==VISCO_LaminarSPS){       
    float dp_sps=(Simulate2D? sqrt(2.0f*Dp*Dp)/2.0f: sqrt(3.0f*Dp*Dp)/3.0f);
    Smag=pow((0.12f*dp_sps),2);
    Blin=(2.0f/3.0f)*0.0066f*dp_sps*dp_sps; 
  }
  
  VisuConfig();
  if(Nfloat){
    if(Simulate2D)RunException(met,"The 2D version is not compatible with floating bodies.");
    if(Kgc)RunException(met,"The KGC option is not compatible with floating bodies.");
    if(TVisco==VISCO_LaminarSPS)RunException(met,"The viscosity Laminar+SPS is not compatible with floating bodies.");
  }
  if(TStep==STEP_Verlet&&Kgc)Log->Print("\n***** Using KGC with Verlet algorithm is not recommended *****\n");
  if(SvDt)InfoDtReset();
}

//==============================================================================
/// Prints out configuration of the case.
//==============================================================================
void JSph::VisuConfig(){
  const char* met="VisuConfig";
  Log->Print(Simulate2D? "**2D-Simulation parameters:": "**3D-Simulation parameters:");
  Log->Print(fun::VarStr("CaseName",CaseName));
  Log->Print(fun::VarStr("RunName",RunName));
  Log->Print(fun::VarStr("StepAlgorithm",GetStepName(TStep)));
  if(TStep==STEP_None)RunException(met,"StepAlgorithm value is invalid.");
  if(TStep==STEP_Verlet)Log->Print(fun::VarStr("VerletSteps",VerletSteps));
  Log->Print(fun::VarStr("Kernel",GetKernelName(TKernel)));
  Log->Print(fun::VarStr("KernelGradientCorr",Kgc));
  Log->Print(fun::VarStr("Viscosity",GetViscoName(TVisco)));
  Log->Print(fun::VarStr("Visco",Visco));
  Log->Print(fun::VarStr("ShepardSteps",ShepardSteps));
  Log->Print(fun::VarStr("DBCSteps",DBCSteps));
  Log->Print(fun::VarStr("Np",Np));
  Log->Print(fun::VarStr("Nbound",Nbound));
  Log->Print(fun::VarStr("Nfixed",Nfixed));
  Log->Print(fun::VarStr("Nmoving",Nmoving));
  Log->Print(fun::VarStr("Nfloat",Nfloat));
  Log->Print(fun::VarStr("Dx",Dp));
  Log->Print(fun::VarStr("H",H));
  Log->Print(fun::VarStr("CteB",CteB));
  Log->Print(fun::VarStr("Gamma",Gamma));
  Log->Print(fun::VarStr("Rhop0",RHOPZERO));
  Log->Print(fun::VarStr("Eps",Eps));
  Log->Print(fun::VarStr("Cs0",Cs0));
  Log->Print(fun::VarStr("CFLnumber",CFLnumber));
  Log->Print(fun::VarStr("DtIni",DtIni));
  Log->Print(fun::VarStr("DtMin",DtMin));
  Log->Print(fun::VarStr("MassFluid",MassFluid));
  Log->Print(fun::VarStr("MassBound",MassBound));
  if(TKernel==KERNEL_Cubic){                                      
    Log->Print(fun::VarStr("CubicCte.a1",CubicCte.a1));   
    Log->Print(fun::VarStr("CubicCte.a2",CubicCte.a2));           
    Log->Print(fun::VarStr("CubicCte.aa",CubicCte.aa));
    Log->Print(fun::VarStr("CubicCte.a24",CubicCte.a24));
    Log->Print(fun::VarStr("CubicCte.c1",CubicCte.c1));
    Log->Print(fun::VarStr("CubicCte.c2",CubicCte.c2));
    Log->Print(fun::VarStr("CubicCte.d1",CubicCte.d1));
    Log->Print(fun::VarStr("CubicCte.od_wdeltap",CubicCte.od_wdeltap));
  }
  if(TKernel==KERNEL_Wendland){                                     
    Log->Print(fun::VarStr("WendlandCte.awen",WendlandCte.awen));
    Log->Print(fun::VarStr("WendlandCte.bwen",WendlandCte.bwen));
  }
  if(TVisco==VISCO_LaminarSPS){       
    Log->Print(fun::VarStr("Smag",Smag));
    Log->Print(fun::VarStr("Blin",Blin));
  }
  if(ShepardSteps)Log->Print(fun::VarStr("CteShepard",CteShepard));                       
  Log->Print(fun::VarStr("TimeMax",TimeMax));
  Log->Print(fun::VarStr("TimePart",TimePart));
  Log->Print(fun::VarStr("Gravity",Gravity));
  Log->Print(fun::VarStr("IncZ",IncZ));
  Log->Print(fun::VarStr("PartOutMax",(int)PartOutMax));
  Log->Print(fun::VarStr("RhopOut",RhopOut));
  if(RhopOut){
    Log->Print(fun::VarStr("RhopOutMin",RhopOutMin));
    Log->Print(fun::VarStr("RhopOutMax",RhopOutMax));
  }
  if(CteB==0)RunException(met,"Constant \'B\' can not be zero. (B=coefsound2*rho0*g*(h_SWL)/gamma is zero when h_SWL=0 (difference of heights in Z-direction)");  
}

//==============================================================================
/// Returns the name of the time algorithm in text.
//==============================================================================
string JSph::GetStepName(TpStep tstep){
  string tx;
  if(tstep==STEP_Verlet)tx="Verlet";
  else if(tstep==STEP_Symplectic)tx="Symplectic";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns the name of the kenel function in text.
//==============================================================================
string JSph::GetKernelName(TpKernel tkernel){
  string tx;
  if(tkernel==KERNEL_Cubic)tx="Cubic";
  else if(tkernel==KERNEL_Wendland)tx="Wendland";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns the name of the viscosity treatment in text.
//==============================================================================
string JSph::GetViscoName(TpVisco tvisco){
  string tx;
  if(tvisco==VISCO_Artificial)tx="Artificial";
  else if(tvisco==VISCO_LaminarSPS)tx="LaminarSPS";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns the output directory.
//==============================================================================
string JSph::GetDirOut(){ return(DirOut); }
int JSph::GetPart(){ return(Part); }
int JSph::GetNdiv(){ return(Ndiv); }

//==============================================================================
/// Returns string with the name of the timer and its value.
//==============================================================================
string JSph::TimerToText(const string &name,float value)const{
  string ret=name;
  while(ret.length()<33)ret+=".";
  return(ret+": "+fun::FloatStr(value/1000)+" sec.");
}
//==============================================================================
/// Returns the name of an activated timer.
//==============================================================================
string JSph::TimerToText(unsigned ct)const{
  return(TimerIsActive(ct)? TimerToText(TimerGetName(ct),TimerGetValue(ct)): "");
}








