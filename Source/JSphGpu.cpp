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

/// \file JSphGpu.cpp \brief Implements the class \ref JSphGpu

#include "JSphGpu.h"
#include "Functions.h"
#include "JSphMotion.h"
#include "JPtxasInfo.h"

//==============================================================================
/// Constructor.
//==============================================================================
JSphGpu::JSphGpu(){
  ClassName="JSphGpu";
  Cpu=false;
  dev=0;
  CsInit(&Dc);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphGpu::~JSphGpu(){
  Reset();
}

//==============================================================================
// Initialisation of variables.
//==============================================================================
void JSphGpu::Reset(){
  JSph::Reset();
  CsReset(&Dc);
  CsResetCuda(dev);
  BlockSizes="";
}

//==============================================================================
/// Computes limits of the domain to be used in the Neighbour List.
//==============================================================================
void JSphGpu::DivideConfig(tfloat3 &posmin,tfloat3 &difmax){
  Log->Print("");
  tfloat3 pmin,pmax;
  pmin.x=Pos[0].x; pmin.y=Pos[0].y; pmin.z=Pos[0].z;
  pmax=pmin;
  for(unsigned c=1;c<Np;c++){
    tfloat3 *pos=(Pos+c);
    if(pmin.x>pos->x)pmin.x=pos->x;
    if(pmin.y>pos->y)pmin.y=pos->y;
    if(pmin.z>pos->z)pmin.z=pos->z;
    if(pmax.x<pos->x)pmax.x=pos->x;
    if(pmax.y<pos->y)pmax.y=pos->y;
    if(pmax.z<pos->z)pmax.z=pos->z;
  }
  float mod=(Dosh/float(Hdiv))*0.05f;
  pmin.x-=mod; pmin.y-=mod; pmin.z-=mod;
  pmax.x+=mod; pmax.y+=mod; pmax.z+=mod;
  float scell=Dosh/Hdiv;
  Ncells.x=int(ceil((pmax.x-pmin.x)/scell));
  Ncells.y=int(ceil((pmax.y-pmin.y)/scell));
  Ncells.z=int(ceil((pmax.z-pmin.z)/scell));
  char tx[256];
  sprintf(tx,"Cells of the initial domain: %d x %d x %d  (%d)",Ncells.x,Ncells.y,Ncells.z,(Ncells.x*Ncells.y*Ncells.z));
  Log->Print(tx);
  pmax.z=(pmax.z-pmin.z)*(IncZ+1.f)+pmin.z;
  Ncells.x=int(ceil((pmax.x-pmin.x)/scell));
  Ncells.y=int(ceil((pmax.y-pmin.y)/scell));
  Ncells.z=int(ceil((pmax.z-pmin.z+mod)/scell)); //-"mod" is used to get a margin due to precision error with real values.
  Nsheet=Ncells.x*Ncells.y;
  Nct=Nsheet*Ncells.z;
  sprintf(tx,"Domain limits: (%.4f,%.4f,%.4f)-(%.4f,%.4f,%.4f)",pmin.x,pmin.y,pmin.z,pmax.x,pmax.y,pmax.z);
  Log->Print(tx);
  posmin=pmin;
  difmax=TFloat3(pmax.x-pmin.x,pmax.y-pmin.y,pmax.z-pmin.z);
  sprintf(tx,"Dimensions: %f x %f x %f",difmax.x,difmax.y,difmax.z);
  Log->Print(tx);
  sprintf(tx,"Cells of the domain: %d x %d x %d",Ncells.x,Ncells.y,Ncells.z);
  Log->Print(tx);
  Log->Print(fun::VarStr("Nsheet",Nsheet));
  Log->Print(fun::VarStr("Nct",Nct));
}

//==============================================================================
/// Allocates memory of main data.
//==============================================================================
void JSphGpu::AllocMemory(){
  const char met[]="AllocMemory";
  char cad[1024];
  bool err=true;
  if(!CellModeIsGPU(CellMode))CellMode=(Simulate2D? CELLMODE_2H: CELLMODE_None);
  if(CellMode==CELLMODE_None){
    CellMode=CELLMODE_Hneigs;
    ConfigDevice();
    CsAllocMemoryBasic(&Dc,&CteVars,SvTimers);
    err=CsAllocMemoryCellMode(&Dc,&CteVars);
    if(err){ 
      sprintf(cad,"Error: GPU memory is not enough. %.2f Mb (%u + %u)",float(GetMemoryGpu())/(1024*1024),Dc.memgpu,Dc.memgpu2);
      Log->Print(cad);
      CellMode=CELLMODE_H; ConfigDevice(); err=CsAllocMemoryCellMode(&Dc,&CteVars);
    }
    if(err){ 
      sprintf(cad,"Error: GPU memory is not enough. %.2f Mb (%u + %u)",float(GetMemoryGpu())/(1024*1024),Dc.memgpu,Dc.memgpu2);
      Log->Print(cad);
      CellMode=CELLMODE_2H; ConfigDevice(); err=CsAllocMemoryCellMode(&Dc,&CteVars);
    }
  }
  else{
    ConfigDevice();
    CsAllocMemoryBasic(&Dc,&CteVars,SvTimers);
    err=CsAllocMemoryCellMode(&Dc,&CteVars);
  }
  if(err){ 
    sprintf(cad,"Error: GPU memory is not enough. %.2f Mb (%u + %u)",float(GetMemoryGpu())/(1024*1024),Dc.memgpu,Dc.memgpu2);
    Log->Print(cad);
    RunException(met,"GPU memory is not enough.");
  }
}

//==============================================================================
/// Processes movement of moving boundary particles.
//==============================================================================
void JSphGpu::RunMotion(float stepdt){
  TmgStart(Dc.timers,TMG_SuMotion);
  if(Motion->ProcesTime(TimeStep+MotionTimeMod,stepdt)){
    unsigned nmove=Motion->GetMovCount();
    if(TStep==STEP_Symplectic)CsSymplecticMotionPre(&Dc,!nmove); //Vel[]=0 and velPre[]=0 (when necessary).
    if(nmove){
      CsCalcRidpmv(&Dc);
      //-Movevement of boundary particles.
      for(unsigned c=0;c<nmove;c++){
        unsigned ref;
        tfloat3 mvsimple;
        tmatrix4f mvmatrix;
        if(Motion->GetMov(c,ref,mvsimple,mvmatrix)){//-Simple movement.
          unsigned pini=MotionObjBegin[ref];
          CsMoveLinBound(&Dc,pini,MotionObjBegin[ref+1]-pini,Float3(mvsimple),Float3(mvsimple.x/stepdt,mvsimple.y/stepdt,mvsimple.z/stepdt));
        }
        else{//-Movement with matrix.
          unsigned pini=MotionObjBegin[ref];
          CsMoveMatBound(&Dc,pini,MotionObjBegin[ref+1]-pini,mvmatrix,stepdt);
        } 
      }
    }
    Dc.bounddatver++;//-NL of boundaries must be updated when there is movement or when it is finished.
                     // to update boundary values of vel_2[] (being velocity zero).
    if(TStep==STEP_Verlet)Dc.verletresetvel=2;//velnew[]=0 for boundaries.
    if(TStep==STEP_Symplectic)CsSymplecticMotionPost(&Dc);//-Copy boundary values of pos[] and vel[] in pospre[] and velpre[].
  }
  TmgStop(Dc.timers,TMG_SuMotion);
}

//==============================================================================
/// Returns the optimum size according to registers of the kernel and compute capability.
//==============================================================================
unsigned JSphGpu::OptimizeBlockSize(unsigned compute,unsigned nreg){
  if(compute>=20){
    if(nreg<=20)return(256);       // 1-20 -> 192:100%  256:100%  384:100%
    else if(nreg<=22)return(224);  //21-22 -> 192:88%  224:88%  256:83%  288:94%  352:92%  480:94%
    else if(nreg<=24)return(224);  //23-24 -> 192:88%  224:88%  256:83%  448:88%
    else if(nreg<=26)return(192);  //25-26 -> 192:75%  256:67%  288:75%  416:81%
    else if(nreg<=28)return(192);  //27-28 -> 192:75%  256:67%  288:75%
    else if(nreg<=30)return(256);  //29-30 -> 256:67%  352:69%
    else if(nreg<=32)return(256);  //31-32 -> 256:67%
    else if(nreg<=34)return(192);  //33-34 -> 192:63%  256:50%
    else if(nreg<=36)return(224);  //35-36 -> 224:58%  256:50%
    else if(nreg<=38)return(160);  //37-38 -> 160:52%  256:50%  416:54%
    else if(nreg<=40)return(160);  //39-40 -> 160:52%  256:50%
    else if(nreg<=42)return(256);  //41-42 -> 256:50%
    else if(nreg<=46)return(224);  //43-46 -> 224:44%  256:33%  352:46%
    else if(nreg<=48)return(224);  //47-48 -> 224:44%  256:33%
    else if(nreg<=50)return(160);  //49-50 -> 160:42%  256:33%  320:42%
    else if(nreg<=56)return(192);  //51-56 -> 192:38%  256:33%  288:38%
    else if(nreg<=64)return(256);  //57-64 -> 128:33%  256:33%
    else if(nreg<=68)return(160);  //65-68 -> 64:29%  96:31%  128:25%  160:31%  192:25%  224:29%  256:17%
    else if(nreg<=72)return(224);  //69-72 -> 64:29%  128:25%  192:25%  224:29%
    else return(128);              //73-84 -> 64:29%  128:25%  192:25%  256:17%
  }
  else if(compute>=12){
    if(nreg<=16)return(256);       // 1-16 -> 128:100%  256:100%
    else if(nreg<=18)return(448);  //17-18 -> 128:75%  192:75%  256:75%  448:88%
    else if(nreg<=20)return(256);  //19-20 -> 192:75%  256:75%  384:75%
    else if(nreg<=21)return(192);  //21    -> 192:75%  256:50%  288:56%  320:63%  352:69%  384:75%
    else if(nreg<=24)return(128);  //22-24 -> 128:63%  192:56%  288:56%  256:50%  320:63%
    else if(nreg<=25)return(320);  //25    -> 192:56%  288:56%  256:50%  320:63%
    else if(nreg<=26)return(192);  //26    -> 192:56%  256:50%
    else if(nreg<=32)return(256);  //27-32 -> 256:50%
    else if(nreg<=36)return(448);  //33-36 -> 192:38%  256:25%  416:41%  448:44%
    else if(nreg<=42)return(192);  //37-42 -> 192:38%  256:25%
    else if(nreg<=51)return(320);  //43-51 -> 256:25%  288:28%  320:31%
    else if(nreg<=64)return(256);  //52-64 -> 128:25%  256:25%
    else return(192);              //65-85 -> 128:13%  192:19%
  }
  else if(compute>=10){
    if(nreg<=10)return(256);       // 1-10 -> 128:100%  192:100%  256:100%  384:100%
    else if(nreg<=12)return(128);  //11-12 -> 128:83%  192:75%  256:67%  320:83%
    else if(nreg<=13)return(192);  //13    -> 128:67%  192:75%  256:67%
    else if(nreg<=16)return(256);  //14-16 -> 128:67%  192:50%  256:67%
    else if(nreg<=18)return(448);  //17-18 -> 128:50%  192:50%  256:33%  384:50%  448:58%
    else if(nreg<=20)return(128);  //19-20 -> 128:50%  192:50%  256:33%  384:50%
    else if(nreg<=21)return(192);  //21    -> 128:33%  192:50%  256:33%  384:50%
    else if(nreg<=24)return(320);  //22-24 -> 64:42%  128:33%  256:33%  320:42%
    else if(nreg<=25)return(320);  //25    -> 64:33%  128:33%  256:33%  320:42%
    else if(nreg<=32)return(256);  //26-32 -> 64:33%  128:33%  256:33%
    else if(nreg<=40)return(192);  //33-40 -> 64:25%  128:17%  192:25%
    else if(nreg<=42)return(192);  //41-42 -> 64:17%  128:17%  192:25%
    else if(nreg<=64)return(128);  //43-64 -> 64:17%  128:17%
    else return(64);               //65-128-> 64:8%
  }
  return(256);
}

//==============================================================================
/// Returns BlockSize according to registers of the kernel.
//==============================================================================
unsigned JSphGpu::BlockSizeConfig(const string& opname,unsigned compute,unsigned regs){
  char cad[1024];
  unsigned bsize=256;
  if(regs){
    bsize=OptimizeBlockSize(compute,regs);
    sprintf(cad,"%s=%u (%u regs)",opname.c_str(),bsize,regs);
  }
  else sprintf(cad,"%s=%u (? regs)",opname.c_str(),bsize);
  Log->Print(cad);
  if(!BlockSizes.empty())BlockSizes=BlockSizes+", ";
  BlockSizes=BlockSizes+cad;
  return(bsize);
}

//==============================================================================
/// Configures data of DeviceContext and DeviceCtes. Returns true in case of error.
//==============================================================================
void JSphGpu::ConfigDevice(){
  const char met[]="ConfigDevice";
  //-Obtains configuration according to CellMode.
  //--------------------------------------
  Log->Print(" ");
  Log->Print("**CUDA kernel configuration:");
  if(!CellModeIsGPU(CellMode))CellMode=CELLMODE_Hneigs;
  Hdiv=(CellMode==CELLMODE_Hneigs||CellMode==CELLMODE_H? 2: 1);
  unsigned ncellnv=(Hdiv*2+1)*(Hdiv*2+1);
  Log->Print(fun::VarStr("CellMode",string(GetNameCellMode(CellMode))));
  if(Hdiv!=1)Log->Print(fun::VarStr("Hdiv",Hdiv));
  Log->Print(fun::VarStr("PtxasFile",PtxasFile));
  Dc.cellmode=CellMode;
  const unsigned sm=(Dc.compute<20? (Dc.compute<12? 10: 12): 20);
  JPtxasInfo pt;
  if(fun::FileExists(PtxasFile))pt.LoadFile(PtxasFile);
  else Log->Print("Without optimization of registers.");
  //if(fun::FileExists(PtxasFile))pt.SaveCsv(DirOut+"PtxasFile.csv");
  BlockSizes="";
  switch(CellMode){
    case CELLMODE_Hneigs:{
      bool floating=(Nfloat>0);
      if(Kgc)Dc.bsmatrixkgc=BlockSizeConfig("BsMatrixKgc",sm,pt.GetRegs("_KerComputeMatrixKgcNeigs",sm,TKernel,ncellnv));
      if(Kgc||TVisco==VISCO_LaminarSPS){
        Dc.bsforcesfluid=BlockSizeConfig("BsInteractionFluid",sm,pt.GetRegs("_KerComputeForcesFullFluidNeigs",sm,TKernel,TVisco,Kgc,true,ncellnv));
        Dc.bsforcesfluidcorr=BlockSizeConfig("BsInteractionFluidCorr",sm,pt.GetRegs("_KerComputeForcesFullFluidNeigs",sm,TKernel,TVisco,Kgc,false,ncellnv));
      }
      else{
        Dc.bsforcesfluid=BlockSizeConfig("BsInteractionFluid",sm,pt.GetRegs("_KerComputeForcesFluidNeigs",sm,TKernel,floating,true,ncellnv));
        Dc.bsforcesfluidcorr=BlockSizeConfig("BsInteractionFluidCorr",sm,pt.GetRegs("_KerComputeForcesFluidNeigs",sm,TKernel,floating,false,ncellnv));
      }
      Dc.bsforcesbound=BlockSizeConfig("BsInteractionBound",sm,pt.GetRegs("_KerComputeForcesBoundNeigs",sm,TKernel,floating,ncellnv));
      Dc.bsshepard=BlockSizeConfig("BsShepard",sm,pt.GetRegs("_KerComputeForcesShepardNeigs",sm,TKernel,ncellnv));
    }break;
    case CELLMODE_2H:
    case CELLMODE_H:{
      bool floating=(Nfloat>0);
      if(Kgc)Dc.bsmatrixkgc=BlockSizeConfig("BsMatrixKgc",sm,pt.GetRegs("_KerComputeMatrixKgc",sm,TKernel,Hdiv));
      if(Kgc||TVisco==VISCO_LaminarSPS){
        Dc.bsforcesfluid=BlockSizeConfig("BsInteractionFluid",sm,pt.GetRegs("_KerComputeForcesFullFluid",sm,TKernel,TVisco,Kgc,true,Hdiv));
        Dc.bsforcesfluidcorr=BlockSizeConfig("BsInteractionFluidCorr",sm,pt.GetRegs("_KerComputeForcesFullFluid",sm,TKernel,TVisco,Kgc,false,Hdiv));
      }
      else{
        Dc.bsforcesfluid=BlockSizeConfig("BsInteractionFluid",sm,pt.GetRegs("_KerComputeForcesFluid",sm,TKernel,floating,true,Hdiv));
        Dc.bsforcesfluidcorr=BlockSizeConfig("BsInteractionFluidCorr",sm,pt.GetRegs("_KerComputeForcesFluid",sm,TKernel,floating,false,Hdiv));
      }
      Dc.bsforcesbound=BlockSizeConfig("BsInteractionBound",sm,pt.GetRegs("_KerComputeForcesBound",sm,TKernel,floating,Hdiv));
      Dc.bsshepard=BlockSizeConfig("BsShepard",sm,pt.GetRegs("_KerComputeForcesShepard",sm,TKernel,Hdiv));
    }break;
    default: RunException(met,"CellMode unrecognised.");
  }

  //-Obtains configuration of NL according to initial data.
  //--------------------------------------------------------
  tfloat3 posmin,difmax;
  DivideConfig(posmin,difmax); //-Perfomed in cpu due to simplicity and efficiency.

  //-Loads data in CteVars (DeviceCtes).
  //--------------------------------------
  CteVars.fourh2=Fourh2; CteVars.h=H;
  CteVars.cubic_a24=CubicCte.a24;                 
  CteVars.cubic_a2=CubicCte.a2;                                                                              
  CteVars.cubic_c1=CubicCte.c1; CteVars.cubic_c2=CubicCte.c2; CteVars.cubic_d1=CubicCte.d1; CteVars.cubic_odwdeltap=CubicCte.od_wdeltap;
  CteVars.wendland_awen=WendlandCte.awen; CteVars.wendland_bwen=WendlandCte.bwen;
  CteVars.cs0=Cs0; CteVars.visco=Visco; CteVars.eta2=Eta2; CteVars.eps=Eps; 
  CteVars.cteshepard=CubicCte.a2;          
  CteVars.massb=MassBound; CteVars.massf=MassFluid;
  CteVars.nbound=Nbound;

  //-Loads data in Dc (DeviceContext).
  //------------------------------------
  Dc.cflnumber=CFLnumber;
  Dc.hdiv=Hdiv;
  Dc.ncellnv=ncellnv;
  Dc.ovscell=1.0f/(Dosh/Hdiv);
  Dc.posmin=Float3(posmin); Dc.difmax=Float3(difmax);
  Dc.ncx=Ncells.x;  Dc.ncy=Ncells.y;  
  Dc.gamma=Gamma; Dc.b=CteB;
  Dc.gravity=Float3(Gravity);

  Dc.nctotmax=Nct+1;
  Dc.simulate2d=(Simulate2D? 1: 0);
  Dc.np=Np; Dc.nbound=Nbound; Dc.nfluid=Nfluid; Dc.nfixed=Nfixed; Dc.nmoving=Nmoving; Dc.nfloat=Nfloat;
  Dc.npok=NpOk; Dc.npb=Npb; Dc.npf=NpOk-Npb; 
  Dc.massbound=MassBound; Dc.massfluid=MassFluid;
  Dc.bounddatver=BoundDatVer; Dc.bounddivver=0;
  Dc.ndiv=0;
  Dc.tkernel=TKernel;
  Dc.tvisco=TVisco;
  Dc.kgc=Kgc;
  Dc.tstep=TStep;
  if(TStep==STEP_Verlet){
    Dc.verletsteps=VerletSteps;
    Dc.verletresetvel=0;
  }
  if(TVisco==VISCO_LaminarSPS){
    Dc.smag=Smag;
    Dc.blin=Blin;
  }
  Dc.shepardsteps=ShepardSteps;
  Dc.dbcsteps=DBCSteps;
  Dc.dtmin=DtMin; Dc.dtmodif=DtModif;
  Dc.infodt=(SvDt? &InfoDt: NULL);
  Dc.rhopout=RhopOut;
  Dc.rhopoutmin=RhopOutMin;
  Dc.rhopoutmax=RhopOutMax;
  Dc.ftobjs=FtObjs;
  Dc.ftcount=FtCount;
}

//==============================================================================
/// Initialisation of arrays and variables for the execution.
//==============================================================================
void JSphGpu::InitVars(){
  CsUpData(&Dc,&CteVars,Idp,(float3*)Pos,(float3*)Vel,Rhop);
  CsUpCte(&CteVars);
  if(TStep==STEP_Verlet)Dc.verletstep=0;
  if(TStep==STEP_Symplectic)Dc.dtpre=DtIni;
  Part=PartIni; Ndiv=0; PartNdiv=0; PartOut=0;
  TimeStep=TimeStepIni; TimeStepM1=TimeStep;
}

//==============================================================================
/// Starts simulation on GPU.
//==============================================================================
void JSphGpu::Run(JCfgRun *cfg,JLog *log){
  const char* met="Run";
  char cad[256],devname[128];
  if(CsInitCuda(cfg->GpuId,&dev,devname,&Dc))RunException(met,"Failure cuda device Initialisation.");
  log->Print("[GPU Hardware]");
  if(cfg->GpuId<0)sprintf(cad,"Device default: %d?  \"%s\"",dev,devname);
  else sprintf(cad,"Device selected: %d \"%s\"",cfg->GpuId,devname);
  log->Print(cad);
  if(cfg->GpuId<0)sprintf(cad,"Gpu_%d?=\"%s\"",dev,devname);
  else sprintf(cad,"Gpu_%d=\"%s\"",cfg->GpuId,devname);
  Hardware=cad;
  sprintf(cad,"Compute capbility: %.1f",float(Dc.compute)/10); log->Print(cad);
  sprintf(cad,"Memory global: %d MB",int(Dc.mglobal/(1024*1024))); log->Print(cad);
  sprintf(cad,"Memory shared: %u Bytes",Dc.mshared); log->Print(cad);

  //-Loads parameters and input data.
  //-------------------------------------------------
  JSph::Run(cfg,log);
  PtxasFile=cfg->PtxasFile;
  
  //-Configures device and allocates memory.
  //-------------------------------------------------
  AllocMemory();  
  //DivideConfig() called from AllocMemory() since it depends on the allocated memory in GPU.
  LoadPartBegin();
  InitVars();

  //-First Neighbour List computation and PART_0000.
  //-------------------------------------------------
  PrintMemoryAlloc();
  CsCallDivide(DV_All,&Dc,&CteVars);
  BoundDatVer=0;
  SaveData(); PartNdiv=-1;
  Part++;

  //-MAIN LOOP.
  //-------------------------------------------------
  bool partoutstop=false;
  TimerSim.Start();
  TimerPart.Start();   //CkCalc=clock();  CkStartPart=CkCalc;
  sprintf(cad,"\n[Initialising simulation  %s]",fun::GetDateTime().c_str()); Log->Print(cad);
  PrintHeadPart();
  while(TimeStep<TimeMax){
    Dc.ndiv=Ndiv;
    float stepdt=CsCallComputeStep(&Dc,&CteVars,(Ndiv%DBCSteps==0));
    if(Nmoving)RunMotion(stepdt);
    CsCallComputeStepDivide(&Dc,&CteVars);
    if(ShepardSteps&&Ndiv&&( !((Ndiv+1)%ShepardSteps) || (TStep==STEP_Verlet && !(Ndiv%ShepardSteps)) ))CsCallRunShepard(&Dc,&CteVars); //-Shepard density filter.
    //-Using Verlet+Shepard, the density filter is applied at the given step and the following one to smooth Rhop and RhopM1.
    TimeStep+=stepdt;
    partoutstop=(Dc.np-Dc.npok)>=PartOutMax;
    if((TimeStep-TimeStepIni)-TimePart*((Part-PartIni)-1)>=TimePart||SvSteps||partoutstop){
      if(partoutstop){
        Log->Print("\n**** Particles OUT limit reached...\n");
        TimeMax=TimeStep;
      }
      SaveData();
      Part++;
      PartNdiv=Ndiv;
      TimeStepM1=TimeStep;
      TimerPart.Start();
    }
    Ndiv++;
  }
  TimerSim.Stop(); TimerTot.Stop();

  //-End of simulation.
  //---------------------------------------------------------------------------
  GetVarsDevice();
  sprintf(cad,"\n[Simulation %s  %s]",(partoutstop? "INTERRUPTED": "finished"),fun::GetDateTime().c_str());  Log->Print(cad);
  sprintf(cad,"DTs adjusted to DtMin............: %d",DtModif);  Log->Print(cad);
  sprintf(cad,"Excluded particles...............: %d",Np-NpOk);  Log->Print(cad);
  if(RhopOut){ sprintf(cad,"Excluded particles due to RhopOut: %u",RhopOutCount);  Log->Print(cad); }
  float tsim=TimerSim.GetElapsedTimeF()/1000.f;
  float ttot=TimerTot.GetElapsedTimeF()/1000.f;
  float tseg=tsim/TimeStep;
  float ndivseg=float(Ndiv)/tsim;
  sprintf(cad,"Total Runtime....................: %f sec.",ttot);    Log->Print(cad);
  sprintf(cad,"Simulation Runtime...............: %f sec.",tsim);    Log->Print(cad);
  sprintf(cad,"Time per second of simulation....: %f sec.",tseg);    Log->Print(cad);
  sprintf(cad,"Steps per second.................: %f",ndivseg);      Log->Print(cad);
  sprintf(cad,"Steps of simulation..............: %d",Ndiv);         Log->Print(cad);
  sprintf(cad,"PART files.......................: %d",Part-PartIni); Log->Print(cad);
  if(SvTimers)ShowTimers();
  if(SvRes)SaveRes(tsim,tseg,ttot,";BlockSizes",string(";")+BlockSizes);
  Log->Print(" ");
}

//==============================================================================
/// Shows active timers.
//==============================================================================
void JSphGpu::ShowTimers(){
  Log->Print("\n[GPU Timers]");
  for(unsigned c=0;c<TimerGetCount();c++){
    string tx=TimerToText(c); if(!tx.empty())Log->Print(tx);
  }
}
//==============================================================================
/// Recovers values of DeviceContext.
//==============================================================================
void JSphGpu::GetVarsDevice(){
  DtModif=Dc.dtmodif;
  NpOk=Dc.npok; Nct=Dc.nct; Ncells.z=Dc.ncz;
  RhopOutCount=Dc.rhopoutcount;
}
//==============================================================================
/// Creates files with output data.
//==============================================================================
void JSphGpu::SaveData(){
  const char met[]="SaveData";
  GetVarsDevice();
  unsigned pini=(BoundDatVer<Dc.bounddatver? 0: Npb);
  CsDownData(&Dc,pini,Idp,(float3*)Pos,(float3*)Vel,Rhop);  
  BoundDatVer=Dc.bounddatver;
  //-Updates Pdata.
  TmgStart(Dc.timers,TMG_SuSavePart);
  unsigned nout=Dc.outcount;
  if(nout)CsOutGetData(&Dc,Idp+NpOk,(float3*)Pos+NpOk,(float3*)Vel+NpOk,Rhop+NpOk);
  if(Pdata.SetDataUnsorted(Part,TimeStep,false,NpOk,nout,Idp,Pos,Vel,Rhop,ProbeVel,ProbeRhop))RunException(met,"Some excluded particles appear again in the simulation.");
  JSph::SaveData();
  TmgStop(Dc.timers,TMG_SuSavePart);
}








