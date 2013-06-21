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

/// \file JCfgRun.cpp \brief Implements the class \ref JCfgRun.

#include "JCfgRun.h"
#include "JSpaceEParms.h"

using namespace std;
using namespace fun;

//==============================================================================
// Constructor.
//==============================================================================
JCfgRun::JCfgRun(){
  ClassName="JCfgRun";
  Reset();
}

//==============================================================================
// Initialisation of variables.
//==============================================================================
void JCfgRun::Reset(){
  PrintInfo=false; SvDef=false; DirsDef=0;
  Cpu=false;
  Gpu=false; GpuId=-1;
  OmpThreads=0;
  OmpMode=OMPM_Single;
  #ifdef USE_OPENMP
    OmpMode=OMPM_Dynamic;
  #endif
  TStep=STEP_None; VerletSteps=-1;
  TKernel=KERNEL_None;
  Kgc=0;                
  TVisco=VISCO_None; Visco=0;
  ShepardSteps=-1;
  DBCSteps=-1;
  SvDt=false; SvRes=true; SvTimers=true;
  Sv_Binx2=false; Sv_Csv=false; Sv_Ascii=false;
  Sv_Vtk=false; Sv_Flw=false;
  CaseName=""; DirOut=""; RunName=""; 
  PartBegin=0; PartBeginFirst=0; PartBeginDir="";
  Incz=-1;
  TimeMax=-1; TimePart=-1;
  RhopOut=true; RhopOutMin=700; RhopOutMax=1300;
  CellMode=CELLMODE_None;
  PtxasFile="";
}

//==============================================================================
/// Shows information about execution parameters.
//==============================================================================
void JCfgRun::VisuInfo()const{
  printf("Information about execution parameters:\n\n");
  printf("  DualSPHysics [name_case [dir_out]] [options]\n\n");
  printf("  Options:\n");
  printf("    -h          Shows information about parameters\n\n");
  printf("    -opt <file> Load a file configuration\n\n");
  printf("    -cpu        Execution on Cpu (option by default)\n");
  printf("    -gpu[:id]   Execucion on Gpu and id of the device\n\n");
#ifdef USE_OPENMP
  printf("    -ompthreads:<int>  Only for Cpu execution, indicates the number of threads\n");
  printf("                   by host for parallel execution, it takes the number of \n");
  printf("                   cores of the device by default (or using zero value)\n");
  printf("    -ompdynamic    Parallel execution with symmetry in interaction \n");
  printf("                   and dynamic load balancing\n");
  printf("    -ompstatic     Parallel execution with symmetry in interaction \n");
  printf("                   and static load balancing\n\n");
#endif
  printf("    -symplectic      Symplectic algorithm as time step algorithm\n");
  printf("    -verlet[:steps]  Verlet algorithm as time step algorithm and number of\n");
  printf("                     time steps to switch equations\n\n");
  printf("    -cubic           Cubic spline kernel\n");
  printf("    -wendland        Wendland kernel\n");
  printf("    -kgc:<0/1>       Kernel Gradient Correction\n\n");            
  printf("    -viscoart:<float>      Artifitical viscosity [0-1]\n");
  printf("    -viscolamsps:<float>   Laminar+SPS viscosity [order of 1E-6]\n");
  printf("\n");
  printf("    -shepard:steps   Shepard filter and number of steps to be applied\n");
  printf("    -dbc:steps       Hughes and Graham correction and number of steps \n");
  printf("                     to update the density of the boundaries\n\n");
  printf("    -sv:[formas,...] Specify the output formats.\n");
  printf("        none    No particles files are generated\n");
  printf("        binx2   Bynary files (option by default)\n");
  printf("        ascii   ASCII files (PART_xxxx of SPHysics)\n");
  printf("        vtk     VTK files\n");
  printf("        csv     CSV files\n");
  printf("        flw     FLW files\n");
  printf("    -svdt:<0/1>      Generate file with information about the time step dt\n");
  printf("    -svres:<0/1>     Generate file that summarizes the execution process\n");
  printf("    -svtimers:<0/1>  Obtain timing for each individual process\n");
  printf("    -name <string>      Specify path and name of the case \n");
  printf("    -runname <string>   Specify name for case execution\n");
  printf("    -dirout <dir>       Specify the out directory \n\n");
  printf("    -partbegin:begin[:first] dir \n");
  printf("     Specify the beginning of the simulation starting from a given PART\n");
  printf("     (begin) and located in the directory (dir), (first) indicates the\n");
  printf("     number of the first PART to be generated\n\n");
  printf("    -incz:<float>    Allowed increase in Z+ direction \n");
  printf("    -rhopout:min:max Exclude fluid particles out of these density limits\n\n");
  printf("    -tmax:<float>   Maximum time of simulation\n");
  printf("    -tout:<float>   Time between output files\n\n");
  printf("    -cellmode:<mode>  Specify the cell division mode. By default, the fastest\n");
  printf("                      approach is chosen \n");
  printf("        hneigs    fastest and the most expensive in memory (only for gpu)\n");
  printf("        h         intermediate\n");
  printf("        2h        slowest and the least expensive in memory \n\n");
  printf("    -ptxasfile <file> Indicate the file with information about the compilation\n");
  printf("     of kernels in CUDA to adjust the size of the blocks depending on the \n");
  printf("     needed registers for each kernel (only for gpu)\n\n");
  printf("  Examples:\n");
  printf("    DualSPHysics case out_case -sv:binx2,csv \n");
  printf("    DualSPHysics case -gpu -svdt:1 \n\n");
}

//==============================================================================
/// Shows current configuration.
//==============================================================================
void JCfgRun::VisuConfig()const{
  printf("\nConfiguration of execution:\n");
  string ln="\n";
  PrintVar("  CaseName",CaseName,ln);
  PrintVar("  RunName",RunName,ln);
  PrintVar("  DirOut",DirOut,ln);
  PrintVar("  PartBegin",PartBegin,ln);
  PrintVar("  PartBeginFirst",PartBeginFirst,ln);
  PrintVar("  PartBeginDir",PartBeginDir,ln);
  PrintVar("  Cpu",Cpu,ln);
  printf("  %s  %s\n",VarStr("Gpu",Gpu).c_str(),VarStr("GpuId",GpuId).c_str());
  PrintVar("  OmpThreads",OmpThreads,ln);
  PrintVar("  OmpMode",GetNameOmpMode(OmpMode),ln);
  PrintVar("  TStep",TStep,ln);
  PrintVar("  VerletSteps",VerletSteps,ln);
  PrintVar("  TKernel",TKernel,ln);
  PrintVar("  Kgc",Kgc,ln);      
  PrintVar("  TVisco",TVisco,ln);
  PrintVar("  Visco",Visco,ln);
  PrintVar("  ShepardSteps",ShepardSteps,ln);
  PrintVar("  DBCSteps",DBCSteps,ln);
  PrintVar("  SvDt",SvDt,ln);
  PrintVar("  SvRes",SvRes,ln);
  PrintVar("  SvTimers",SvTimers,ln);
  PrintVar("  Sv_Binx2",Sv_Binx2,ln);
  PrintVar("  Sv_Csv",Sv_Csv,ln);
  PrintVar("  Sv_Ascii",Sv_Ascii,ln);
  PrintVar("  Sv_Vtk",Sv_Vtk,ln);
  PrintVar("  Sv_Flw",Sv_Flw,ln);
  PrintVar("  Incz",Incz,ln);
  PrintVar("  RhopOut",RhopOut,ln);
  PrintVar("  RhopOutMin",RhopOutMin,ln);
  PrintVar("  RhopOutMax",RhopOutMax,ln);
  PrintVar("  TimeMax",TimeMax,ln);
  PrintVar("  TimePart",TimePart,ln);
  PrintVar("  CellMode",GetNameCellMode(CellMode),ln);
  PrintVar("  PtxasFile",PtxasFile,ln);
}

//==============================================================================
/// Loads execution parameters from the command line.
//==============================================================================
void JCfgRun::LoadArgv(int argc,char** argv){
  const char met[]="LoadArgv";
  Reset();
  const int MAXOPTS=100;
  string *optlis=new string[MAXOPTS];
  int optn=0;
  for(int c=0;c<argc-1;c++){
    string tex=StrTrim(argv[c+1]);
    int pos=int(tex.find(" "));
    if(pos>0){
      while(pos>0){
        if(optn>=MAXOPTS)RunException(met,"Has exceeded the maximum configuration options.");
        optlis[optn]=tex.substr(0,pos); optn++;
        tex=tex.substr(pos+1);
        pos=int(tex.find(" "));
      }
    }
    if(optn>=MAXOPTS)RunException(met,"Has exceeded the maximum configuration options.");
    optlis[optn]=tex; optn++;
  }
  if(optn)LoadOpts(optlis,optn,0,"");
  delete[] optlis;
  if(!optn)PrintInfo=true;
  if(!PrintInfo){ //-Configuration by default.
    if(!Cpu&&!Gpu)Cpu=true;
    if(!SvDef)Sv_Binx2=true;
  }
  else VisuInfo();
  if(PtxasFile.empty())PtxasFile=GetWithoutExtension(argv[0])+"_ptxasinfo";
}

//==============================================================================
/// Loads execution parameters from a text file.
//==============================================================================
void JCfgRun::LoadFile(string fname,int lv){
  const char met[]="LoadFile";
  //printf("\nFile:[%s] lv:%d\n",fname.c_str(),lv);
  const int MAXOPTS=50;
  int optn=0;
  string *optlis=new string[MAXOPTS];
  ifstream pf;
  pf.open(fname.c_str());
  if(pf){
    while(!pf.eof()&&optn<MAXOPTS){
      string tex;  pf >> tex;
      if(tex!=""){
        if(optn<MAXOPTS)optlis[optn]=tex;
        optn++;
      }
    } 
    if(!pf.eof()&&pf.fail())RunException(met,"Error reading data from the file.",fname);
    pf.close();
  }
  else RunException(met,"The file can not be opened.",fname);
  if(optn>=MAXOPTS){
    char cad[128];
    sprintf(cad,"File with too many lines (Maximum=%d)",MAXOPTS);
    RunException(met,cad,fname);
  }
  if(optn>0)LoadOpts(optlis,optn,lv,fname);
  delete[] optlis;
}

//==============================================================================
/// Generates error of unknown parameter.
//==============================================================================
void JCfgRun::ErrorParm(const std::string &opt,int optc,int lv,const std::string &file)const{
  const char met[]="ErrorParm";
  char cad[256];
  std::string tx;
  sprintf(cad,"Parameter \"%s\" unrecognised or invalid. ",opt.c_str());  tx=cad;
  sprintf(cad,"(Level cfg:%d, Parameter:%d)",lv,optc);  tx=tx+cad;
  if(file!="")RunException(met,tx,file); else RunException(met,tx);
}

//==============================================================================
/// Loads execution parameters.
//==============================================================================
void JCfgRun::LoadOpts(string *optlis,int optn,int lv,string file){
  const char met[]="LoadOpts";
  if(lv>=10)RunException(met,"No more than 10 levels of recursive configuration.");
  for(int c=0;c<optn;c++){
    string opt=optlis[c];
    if(opt[0]!='-'&&opt[0]!='#'){
      if(!DirsDef){ CaseName=opt; DirsDef++; }
      else if(DirsDef==1){ DirOut=opt; DirsDef++; }
      else ErrorParm(opt,c,lv,file);
    }
    else if(opt[0]=='-'){
      string tx=opt.substr(1);
      int pos=int(tx.find("#"));
      if(pos>0)tx=tx.substr(0,pos);
      pos=int(tx.find(":"));
      string txopt="",txopt2="",txopt3="",txword=StrUpper(pos>0? tx.substr(0,pos): tx);
      if(pos>=0)txopt=tx.substr(pos+1);
      tx=txopt;
      pos=int(tx.find(":"));
      txopt=(pos>=0? tx.substr(0,pos): tx);
      if(pos>=0)txopt2=tx.substr(pos+1);
      tx=txopt2;
      pos=int(tx.find(":"));
      txopt2=(pos>=0? tx.substr(0,pos): tx);
      if(pos>=0)txopt3=tx.substr(pos+1);
      if(txword=="CPU"){ Cpu=true; Gpu=false; }
      else if(txword=="GPU"){ Gpu=true; Cpu=false;
        if(txopt!="")GpuId=atoi(txopt.c_str()); 
      }
#ifdef USE_OPENMP
      else if(txword=="OMPTHREADS"){ 
        OmpThreads=atoi(txopt.c_str()); if(OmpThreads<0)OmpThreads=0;
        if(OmpThreads==1)OmpMode=OMPM_Single;
        else if(OmpMode==OMPM_Single)OmpMode=OMPM_Dynamic;
      } 
      else if(txword=="OMPDYNAMIC")OmpMode=OMPM_Dynamic;
      else if(txword=="OMPSTATIC")OmpMode=OMPM_Static;
#endif
      else if(txword=="SYMPLECTIC")TStep=STEP_Symplectic;
      else if(txword=="VERLET"){ TStep=STEP_Verlet; 
        if(txopt!="")VerletSteps=atoi(txopt.c_str()); 
      }
      else if(txword=="CUBIC")TKernel=KERNEL_Cubic;
      else if(txword=="WENDLAND")TKernel=KERNEL_Wendland;
      else if(txword=="KGC"){                 
        bool on=(txopt!=""? atoi(txopt.c_str()): 1)!=0;
        Kgc=(on? 1: 2);
      }
      else if(txword=="VISCOART"){ 
        Visco=float(atof(txopt.c_str())); 
        if(Visco>1)ErrorParm(opt,c,lv,file);
        TVisco=VISCO_Artificial;
      }
      else if(txword=="VISCOLAMSPS"){ 
        Visco=float(atof(txopt.c_str())); 
        if(Visco>0.001)ErrorParm(opt,c,lv,file);
        TVisco=VISCO_LaminarSPS;
      }
      else if(txword=="SHEPARD"){
        if(txopt!="")ShepardSteps=atoi(txopt.c_str()); 
        if(ShepardSteps<0)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="DBC"){
        if(txopt!="")DBCSteps=atoi(txopt.c_str()); 
        if(DBCSteps<1)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="SVDT")SvDt=(txopt!=""? atoi(txopt.c_str()): 1)!=0;
      else if(txword=="SVRES")SvRes=(txopt!=""? atoi(txopt.c_str()): 1)!=0;
      else if(txword=="SVTIMERS")SvTimers=(txopt!=""? atoi(txopt.c_str()): 1)!=0;
      else if(txword=="SV"){
        string txop=StrUpper(txopt);
        while(txop.length()>0){
          pos=int(txop.find(","));
          string op=(pos>=0? txop.substr(0,pos): txop);
          txop=(pos>=0? txop.substr(pos+1): "");
          if(op=="NONE"){ 
            SvDef=true; Sv_Ascii=false; Sv_Binx2=false; 
            Sv_Csv=false; Sv_Vtk=false; Sv_Flw=false;
          }
          else if(op=="BINX2"){   SvDef=true; Sv_Binx2=true; }
          else if(op=="ASCII"){   SvDef=true; Sv_Ascii=true; }
          else if(op=="CSV"){     SvDef=true; Sv_Csv=true; }
          else if(op=="VTK"){     SvDef=true; Sv_Vtk=true; }
          else if(op=="FLW"){     SvDef=true; Sv_Flw=true; }
          else ErrorParm(opt,c,lv,file);
        }
      }
      else if(txword=="NAME"&&c+1<optn){ CaseName=optlis[c+1]; c++; }
      else if(txword=="RUNNAME"&&c+1<optn){ RunName=optlis[c+1]; c++; }
      else if(txword=="DIROUT"&&c+1<optn){ DirOut=optlis[c+1]; c++; }
      else if(txword=="PARTBEGIN"&&c+1<optn){ 
        int v1=atoi(txopt.c_str());
        int v2=atoi(txopt2.c_str());
        if(v1<0||v2<0)ErrorParm(opt,c,lv,file);
        else{
          PartBegin=unsigned(v1);
          PartBeginFirst=(txopt2.empty()? PartBegin: unsigned(v2));
        }
        PartBeginDir=optlis[c+1]; c++; 
      }
      else if(txword=="INCZ"){ 
        Incz=float(atof(txopt.c_str())); 
        if(Incz>25||Incz<0)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="RHOPOUT"){ 
        RhopOutMin=float(atof(txopt.c_str())); 
        RhopOutMax=float(atof(txopt2.c_str())); 
        RhopOut=(RhopOutMin<RhopOutMax);
      }
      else if(txword=="TMAX"){ 
        TimeMax=float(atof(txopt.c_str())); 
        if(TimeMax<0)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="TOUT"){ 
        TimePart=float(atof(txopt.c_str())); 
        if(TimePart<0)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="CELLMODE"){
        bool ok=true;
        if(!txopt.empty()){
          txopt=StrUpper(txopt);
          if(txopt=="HNEIGS")CellMode=CELLMODE_Hneigs;
          else if(txopt=="H")CellMode=CELLMODE_H;
          else if(txopt=="2H")CellMode=CELLMODE_2H;
          else ok=false;
        }
        else ok=false;
        if(!ok)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="PTXASFILE"&&c+1<optn){ PtxasFile=optlis[c+1]; c++; }
      else if(txword=="OPT"&&c+1<optn){ LoadFile(optlis[c+1],lv+1); c++; }
      else if(txword=="H"||txword=="HELP"||txword=="?")PrintInfo=true;
      else ErrorParm(opt,c,lv,file);
    }
  }
}







