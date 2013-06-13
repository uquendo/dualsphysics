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

/** \mainpage DualSPHysics Documentation
\section main-des Description
DualSPHysics is based on the Smoothed Particle Hydrodynamics <a href="http://www.sphysics.org">SPHysics code.</a> \n
The package is a set of C++ and CUDA codes. \n
DualSPHysics is developed to deal with real-life engineering problems <a href="http://www.youtube.com/user/DualSPHysics">DualSPHysics animations.</a> \n

EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain. \n
School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.
\section compile_sec Project files
Please download source files and documentation from <a href="http://dual.sphysics.org">DualSPHysics website.</a> \n
\author Jose M. Dominguez
\author Dr. Alejandro Crespo 
\author Prof. M. Gomez Gesteira
\author Anxo Barreiro
\author Dr. Benedict Rogers
\version 2.00
\date 15-03-2012
\copyright GNU Public License
*/

/// \file main.cpp \brief Main file of the project that executes the code on CPU or GPU

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;
#include "JLog.h"
#include "JCfgRun.h"
#include "JSphCpu.h"
#include "JException.h"
#ifndef _ONLY_CPU
  #include "JSphGpu.h"
#endif

//==============================================================================
//==============================================================================
std::string getlicense(const char* name){
  std::string tx="";
  tx=tx+"\n\n <"+name+">  Copyright (C) 2012 by Jose M. Dominguez, Dr. Alejandro ";
  tx=tx+"\n Crespo, Prof. M. Gomez Gesteira, Anxo Barreiro, Dr. Benedict Rogers \n";
  tx=tx+"\n DualSPHysics code belongs to SPHysics project, www.sphysics.org,  ";
  tx=tx+"\n an international collaboration between University of Vigo (Spain), ";
  tx=tx+"\n University of Manchester (UK) and Johns Hopkins University (USA) \n";
  tx=tx+"\n This file is part of DualSPHysics project. \n";
  tx=tx+"\n DualSPHysics is free software: you can redistribute it and/or modify";
  tx=tx+"\n it under the terms of the GNU General Public License as published by";
  tx=tx+"\n the Free Software Foundation, either version 3 of the License, or";
  tx=tx+"\n (at your option) any later version. \n";
  tx=tx+"\n DualSPHysics is distributed in the hope that it will be useful,";
  tx=tx+"\n but WITHOUT ANY WARRANTY; without even the implied warranty of";
  tx=tx+"\n MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the";
  tx=tx+"\n GNU General Public License for more details. \n";
  tx=tx+"\n You should have received a copy of the GNU General Public License,";
  tx=tx+"\n along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. \n\n";
  return(tx);
}

//==============================================================================
//==============================================================================
int main(int argc, char** argv)
{
  int errcode=1;
  std::string license=getlicense("DUALSPHYSICS");
  printf("%s",license.c_str());
  char appname[256],appnamesub[256];
  sprintf(appname,"DualSPHysics v%s (08-03-2012)",JSph::GetVersionStr().c_str());
  for(unsigned c=0;c<=strlen(appname);c++)appnamesub[c]='='; appnamesub[strlen(appname)+1]='\0';
  printf("\n%s\n%s\n",appname,appnamesub);
  JCfgRun cfg;
  JLog log;
  try{
    cfg.LoadArgv(argc,argv);
    if(!cfg.PrintInfo){
      log.Init(cfg.DirOut+"/Run.out");
      log.Print(license,JLog::Out_File);
      log.Print(appname,JLog::Out_File);
      log.Print(appnamesub,JLog::Out_File);
      //- SPH Execution
#ifdef _ONLY_CPU
      cfg.Cpu=true;
#endif
      if(cfg.Cpu){ 
        JSphCpu sph;
        sph.Run(&cfg,&log);
      }
#ifndef _ONLY_CPU
      else{
        JSphGpu sph;
        sph.Run(&cfg,&log);
      }
#endif
    }
    errcode=0;
  }
  catch(char *cad){
    string tx=string("\n*** Excepcion: ")+cad+"\n";
    if(log.IsOk())log.Print(tx); else printf("%s",tx.c_str());
  }
  catch(const string &e){
    string tx=string("\n*** Excepcion: ")+e+"\n";
    if(log.IsOk())log.Print(tx); else printf("%s",tx.c_str());
  }
  catch (const exception &e){
    string tx=string("\n*** ")+e.what()+"\n";
    if(log.IsOk())log.Print(tx); else printf("%s",tx.c_str());
  }
  catch(...){
    printf("\n*** Atention: Unknown exception...\n");
  }
  return(errcode);
}








