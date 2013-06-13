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

/// \file JPtxasInfo.h \brief Declares the class \ref JPtxasInfo

#ifndef _JPtxasInfo_
#define _JPtxasInfo_

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

#include "JObject.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>


class JPtxasInfoKer : protected JObject
{
public:
  typedef struct{
    std::string type;
    std::string value;
  }StTemplateArg; 

protected:
  std::string Code,CodeMin;
  std::string NameSpace;
  std::string Name;
  std::string Args;

  std::vector<StTemplateArg*> TemplateArgs;  

  unsigned Regs_sm10; ///<Number of registers according to compute capability 1.0 (0:No info).
  unsigned Regs_sm12; ///<Number of registers according to compute capability 1.2 (0:No info).
  unsigned Regs_sm20; ///<Number of registers according to compute capability 2.0 (0:No info).
  
  void UpdateCode();

public:

  JPtxasInfoKer();
  JPtxasInfoKer(const JPtxasInfoKer &ker);
  ~JPtxasInfoKer(){ Reset(); }
  JPtxasInfoKer & operator=(const JPtxasInfoKer &ker);

  void Reset();
  void SetName(const std::string &name){ Name=name; UpdateCode(); }
  void SetNameSpace(const std::string &namesp){ NameSpace=namesp; UpdateCode(); }
  void SetArgs(const std::string &args){ Args=args; UpdateCode(); }
  void AddTemplateArg(const std::string &type,const std::string &value);
  void SetRegs(unsigned sm,unsigned regs);
  
  std::string GetCode()const{ return(Code); };
  std::string GetCodeMin()const{ return(CodeMin); };
  std::string GetName()const{ return(Name); };
  std::string GetNameSpace()const{ return(NameSpace); };
  std::string GetArgs()const{ return(Args); };
  unsigned GetRegs(unsigned sm)const;
  unsigned CountTemplateArgs()const{ return(unsigned(TemplateArgs.size())); }
  std::string GetTemplateArgsType(unsigned pos)const{ return(pos<CountTemplateArgs()? TemplateArgs[pos]->type: std::string("")); }
  std::string GetTemplateArgsValue(unsigned pos)const{ return(pos<CountTemplateArgs()? TemplateArgs[pos]->value: std::string("")); }

  void Print()const;
};

//##############################################################################
//# JPtxasInfo
//##############################################################################
/// \brief Returns the number of registers of each CUDA kernel. 

class JPtxasInfo : protected JObject
{
protected:
  static const unsigned SM_COUNT=3;
  unsigned SmValues[SM_COUNT];
  std::vector<JPtxasInfoKer*> Kernels;  

public:

  JPtxasInfo();
  void Reset();
  void LoadFile(const std::string &file);
  void AddKernel(const JPtxasInfoKer &ker);
  int GetIndexKernel(const std::string &code)const;

  unsigned Count()const{ return(unsigned(Kernels.size())); }
  const JPtxasInfoKer* GetKernel(unsigned pos)const;
  void SaveCsv(const std::string &file)const;

  unsigned GetRegs(const std::string &kername,unsigned sm)const;
  unsigned GetRegs(const std::string &kername,unsigned sm,unsigned v1)const;
  unsigned GetRegs(const std::string &kername,unsigned sm,unsigned v1,unsigned v2)const;
  unsigned GetRegs(const std::string &kername,unsigned sm,unsigned v1,unsigned v2,unsigned v3)const;
  unsigned GetRegs(const std::string &kername,unsigned sm,unsigned v1,unsigned v2,unsigned v3,unsigned v4)const;
  unsigned GetRegs(const std::string &kername,unsigned sm,unsigned v1,unsigned v2,unsigned v3,unsigned v4,unsigned v5)const;

  void Sort();

  void Print()const;
};

#endif







