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

/// \file JPtxasInfo.cpp \brief Implements the class \ref JPtxasInfo

#include "JPtxasInfo.h"
#include <algorithm>
#include <cstring>

using std::string;
using std::ifstream;
using std::ofstream;
using std::endl;

//==============================================================================
/// Constructor.
//==============================================================================
JPtxasInfoKer::JPtxasInfoKer(){
  ClassName="JPtxasInfoKer";
  Reset();
}
//==============================================================================
JPtxasInfoKer::JPtxasInfoKer(const JPtxasInfoKer &ker){
  ClassName="JPtxasInfoKer";
  Reset();
  *this=ker;
}

//==============================================================================
/// Overloading operator for a correct allocation.
//==============================================================================
JPtxasInfoKer &JPtxasInfoKer::operator=(const JPtxasInfoKer &ker){
  SetName(ker.GetName());
  SetNameSpace(ker.GetNameSpace());
  SetArgs(ker.GetArgs());
  SetRegs(10,ker.GetRegs(10));
  SetRegs(12,ker.GetRegs(12));
  SetRegs(20,ker.GetRegs(20));
  for(unsigned c=0;c<ker.CountTemplateArgs();c++)AddTemplateArg(ker.GetTemplateArgsType(c),ker.GetTemplateArgsValue(c));
  return(*this);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPtxasInfoKer::Reset(){
  Code=""; CodeMin="";
  NameSpace=""; Name=""; Args="";
  Regs_sm10=0; Regs_sm12=0; Regs_sm20=0;
  for(unsigned c=0;c<TemplateArgs.size();c++)delete TemplateArgs[c];
  TemplateArgs.clear();
}

//==============================================================================
/// Adds an argument to the template.
//==============================================================================
void JPtxasInfoKer::UpdateCode(){
  CodeMin=NameSpace+"_"+Name;
  for(unsigned c=0;c<TemplateArgs.size();c++)CodeMin=CodeMin+"_"+TemplateArgs[c]->value;
  Code=CodeMin+"_"+Args;
}

//==============================================================================
/// Adds an argument to the template.
//==============================================================================
void JPtxasInfoKer::AddTemplateArg(const std::string &type,const std::string &value){
  StTemplateArg* arg=new StTemplateArg;
  arg->type=type; arg->value=value;
  TemplateArgs.push_back(arg);
  UpdateCode();
}

//==============================================================================
/// Modifies the number of registers according to sm.
//==============================================================================
void JPtxasInfoKer::SetRegs(unsigned sm,unsigned regs){
  if(sm==10)Regs_sm10=regs;
  else if(sm==12)Regs_sm12=regs;
  else if(sm==20)Regs_sm20=regs;
}

//==============================================================================
/// Returns the number of registers according to sm.
//==============================================================================
unsigned JPtxasInfoKer::GetRegs(unsigned sm)const{
  if(sm==10)return(Regs_sm10);
  else if(sm==12)return(Regs_sm12);
  else if(sm==20)return(Regs_sm20);
  return(0);
}

//==============================================================================
/// Shows data for debug.
//==============================================================================
void JPtxasInfoKer::Print()const{
  printf("JPtxasInfoKer{\n");
  printf("  Code=[%s]\n",Code.c_str());
  printf("  Name=[%s]\n",Name.c_str());
  printf("  Args=[%s]\n",Args.c_str());
  printf("  Regs_sm10=%u\n",Regs_sm10);
  printf("  Regs_sm12=%u\n",Regs_sm12);
  printf("  Regs_sm20=%u\n",Regs_sm20);
  for(unsigned c=0;c<TemplateArgs.size();c++)printf("  TemplateArgs[%u]=[%s]:[%s]\n",c,TemplateArgs[c]->type.c_str(),TemplateArgs[c]->value.c_str());
}

//==============================================================================
//##############################################################################
//==============================================================================

//==============================================================================
/// Constructor.
//==============================================================================
JPtxasInfo::JPtxasInfo(){
  ClassName="JPtxasInfo";
  SmValues[0]=10;
  SmValues[1]=12;
  SmValues[2]=20;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPtxasInfo::Reset(){
  for(unsigned c=0;c<Kernels.size();c++)delete Kernels[c];
  Kernels.clear();
}

//==============================================================================
/// Returns the index of the requested kernel (-1: It does not exist)
//==============================================================================
int JPtxasInfo::GetIndexKernel(const std::string &code)const{
  int pos=-1;
  for(unsigned c=0;c<Kernels.size()&&pos<0;c++)if(Kernels[c]->GetCode()==code)pos=c;
  return(pos);
}

//==============================================================================
/// Returns the pointer to the requested kernel.
//==============================================================================
const JPtxasInfoKer* JPtxasInfo::GetKernel(unsigned pos)const{ 
  const char met[]="GetKernel";
  if(pos>=Count())RunException(met,"Kernel number requested is not valid.");
  return(Kernels[pos]);
}

//==============================================================================
/// Adds a new kernel.
//==============================================================================
void JPtxasInfo::AddKernel(const JPtxasInfoKer &ker){
  int pos=GetIndexKernel(ker.GetCode());
  if(pos<0){
    JPtxasInfoKer* k=new JPtxasInfoKer(ker);
    Kernels.push_back(k);
    pos=int(Kernels.size())-1;
  }
  for(unsigned c=0;c<SM_COUNT;c++){
    unsigned sm=SmValues[c];
    if(ker.GetRegs(sm))Kernels[pos]->SetRegs(sm,ker.GetRegs(sm));
  }
}

//==============================================================================
/// Loads ptxas info from the file obtained after compilation.
//==============================================================================
void JPtxasInfo::LoadFile(const std::string &file){
  const char met[]="LoadFile";
  ifstream pf;
  pf.open(file.c_str());
  if(pf){
    JPtxasInfoKer ker;
    unsigned kersm;
    bool fin=false;
    while(!pf.eof()&&!fin){
      char buff[1024];
      pf.getline(buff,1024);
      string line=buff;
      int pos=int(line.find("ptxas info"));
      if(pos>=0&&int(line.find("Function properties"))<0){
        pos=int(line.find("Compiling entry function '_Z")); 
        int pos2=int(line.find("' for 'sm_")); 
        if(pos>=0&&pos2>=0){//-Name of the kernel.
          ker.Reset(); kersm=0;
          //-Obtains the compute capability.
          string tx=line.substr(pos2+string("' for 'sm_").length());
          int len=int(strspn(tx.c_str(),"0123456789"));
          kersm=atoi((tx.substr(0,len)).c_str());
          //-Obtains the name of namespace.
          tx=line.substr(0,pos2);
          tx=tx.substr(pos+string("Compiling entry function '_Z").length());
          if(!tx.empty()&&tx[0]=='N'){
            tx=tx.substr(1);
            len=int(strspn(tx.c_str(),"0123456789"));
            int n=atoi((tx.substr(0,len)).c_str());
            ker.SetNameSpace(tx.substr(len,n));
            tx=tx.substr(len+n);
          }
          //-Obtains name of the function.
          len=int(strspn(tx.c_str(),"0123456789"));
          int n=atoi((tx.substr(0,len)).c_str());
          ker.SetName(tx.substr(len,n));
          tx=tx.substr(len+n);
          //-Obtains parameters of the template.
          if(!tx.empty()&&tx[0]=='I'){
            tx=tx.substr(1);
            while(!tx.empty()&&tx[0]=='L'){
              tx=tx.substr(1);
              int len=int(strspn(tx.c_str(),"0123456789"));
              //-Obtains type of argument.
              string typearg;
              if(len){//-Type with name.
                int n=atoi((tx.substr(0,len)).c_str());
                typearg=tx.substr(len,n);  tx=tx.substr(len+n);
              }
              else{//-Basic type (b:bool, j:unsigned, i:int,...)
                typearg=tx[0];  tx=tx.substr(1);
              }
              //-Obtains value of argument.
              pos=int(tx.find("E")); 
              if(pos<0)RunException(met,"Error in interpreting the template arguments.",file);
              ker.AddTemplateArg(typearg,tx.substr(0,pos));
              tx=tx.substr(pos+1);
            }
            tx=tx.substr(1);
          }
          ker.SetArgs(tx);
          AddKernel(ker);
        }
        else if(kersm){
          pos=int(line.find("Used ")); 
          int pos2=int(line.find(" registers")); 
          if(pos>=0&&pos2>=0){//-Registers of the kernel.
            pos+=int(string("Used ").length());
            string tx=line.substr(pos,pos2-pos);
            ker.SetRegs(kersm,atoi(tx.c_str()));
            AddKernel(ker);
          }
          ker.Reset(); kersm=0;
        }
      }
    } 
    if(!pf.eof()&&pf.fail())RunException(met,"Failed reading to file.",file);
    pf.close();
  }
  else RunException(met,"File could not be opened.",file);
  Sort();
}

//==============================================================================
/// Function to order kernels.
//==============================================================================
bool JPtxasInfoKerSort(JPtxasInfoKer* i,JPtxasInfoKer* j){ 
  return (i->GetCode()<j->GetCode());
}

//==============================================================================
/// Sorts kernels by code.
//==============================================================================
void JPtxasInfo::Sort(){
  sort(Kernels.begin(),Kernels.end(),JPtxasInfoKerSort);
}

//==============================================================================
/// Shows data for debug.
//==============================================================================
void JPtxasInfo::Print()const{
  for(unsigned c=0;c<Count();c++){
    printf("\nJPtxasInfoKer[%u]\n",c);
    GetKernel(c)->Print();
  }
}

//==============================================================================
/// Stores data in csv file.
//==============================================================================
void JPtxasInfo::SaveCsv(const std::string &file)const{
  const char met[]="SaveCsv";
  ofstream pf;
  pf.open(file.c_str());
  if(pf){
    pf << "Sp;Kernel;TemplateVars;SM_10;SM_12;SM_20;Args" << endl;
    for(unsigned c=0;c<Count();c++){
      JPtxasInfoKer* ker=Kernels[c];
      string tmpvar;
      for(unsigned t=0;t<ker->CountTemplateArgs();t++){
        if(t)tmpvar=tmpvar+", ";
        tmpvar=tmpvar+ker->GetTemplateArgsType(t)+"="+ker->GetTemplateArgsValue(t);
      }
      pf << ker->GetNameSpace() << ";" << ker->GetName() << ";" << tmpvar << ";";
      pf << ker->GetRegs(10) << ";" << ker->GetRegs(12) << ";" << ker->GetRegs(20) << ";";
      pf << ker->GetArgs() << endl;
    }   
    if(pf.fail())RunException(met,"Failed writing to file.",file);
    pf.close();
  }
  else RunException(met,"File could not be opened.",file);
}

//==============================================================================
/// Returns the number of registers of the requested kernel.
//==============================================================================
unsigned JPtxasInfo::GetRegs(const std::string &kername,unsigned sm)const{
  unsigned reg=0,n=0;
  for(unsigned c=0;c<Count();c++){
    JPtxasInfoKer* ker=Kernels[c];
    if(ker->GetCodeMin()==kername){ reg=ker->GetRegs(sm); n++; }
  }
  return(n!=1? 0: reg);
}
//==============================================================================
unsigned JPtxasInfo::GetRegs(const std::string &kername,unsigned sm,unsigned v1)const{
  char cad[1024]; sprintf(cad,"_%u",v1);
  return(GetRegs(kername+cad,sm));
}
//==============================================================================
unsigned JPtxasInfo::GetRegs(const std::string &kername,unsigned sm,unsigned v1,unsigned v2)const{
  char cad[1024]; sprintf(cad,"_%u_%u",v1,v2);
  return(GetRegs(kername+cad,sm));
}
//==============================================================================
unsigned JPtxasInfo::GetRegs(const std::string &kername,unsigned sm,unsigned v1,unsigned v2,unsigned v3)const{
  char cad[1024]; sprintf(cad,"_%u_%u_%u",v1,v2,v3);
  return(GetRegs(kername+cad,sm));
}
//==============================================================================
unsigned JPtxasInfo::GetRegs(const std::string &kername,unsigned sm,unsigned v1,unsigned v2,unsigned v3,unsigned v4)const{
  char cad[1024]; sprintf(cad,"_%u_%u_%u_%u",v1,v2,v3,v4);
  return(GetRegs(kername+cad,sm));
}
//==============================================================================
unsigned JPtxasInfo::GetRegs(const std::string &kername,unsigned sm,unsigned v1,unsigned v2,unsigned v3,unsigned v4,unsigned v5)const{
  char cad[1024]; sprintf(cad,"_%u_%u_%u_%u_%u",v1,v2,v3,v4,v5);
  return(GetRegs(kername+cad,sm));
}








