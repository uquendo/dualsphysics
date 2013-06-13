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

/// \file JVarsAscii.cpp \brief Implements the class \ref JVarsAscii.

#include "JVarsAscii.h"

using std::string;
using std::ifstream;
using std::ofstream;
using std::endl;

//==============================================================================
/// Constructor.
//==============================================================================
JVarsAscii::JVarsAscii(){
  ClassName="JVarsAscii";
  Names=NULL; Values=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JVarsAscii::~JVarsAscii(){
  delete[] Names;  Names=NULL;
  delete[] Values; Values=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JVarsAscii::Reset(){
  delete[] Names;  Names=new string[MAXVARS];
  delete[] Values; Values=new string[MAXVARS];
  Count=0;
}

//==============================================================================
/// Gets a string without spaces.
//==============================================================================
std::string JVarsAscii::StrTrim(const std::string &cad){  //static
  string ret;
  int lsp=0,rsp=0;
  for(int c=0;c<int(cad.length())&&cad[c]==' ';c++)lsp++;
  for(int c=int(cad.length())-1;c<int(cad.length())&&cad[c]==' ';c--)rsp++;
  int size=int(cad.length())-(lsp+rsp);
  return(size>0? cad.substr(lsp,size): "");
}

//==============================================================================
/// Gets string in uppercase.
//==============================================================================
std::string JVarsAscii::StrUpper(const std::string &cad){  //static
  string ret;
  for(unsigned c=0;c<cad.length();c++)ret=ret+char(toupper(cad[c]));
  return(ret);
}

//==============================================================================
/// Stores data in a file.
//==============================================================================
void JVarsAscii::SaveFile(const std::string &file)const{
  const char met[]="SaveFile";
  ofstream pf;
  pf.open(file.c_str());
  if(pf){
    for(unsigned c=0;c<GetCount();c++)pf << GetName(c) << "=" << GetValue(c) << endl;
    if(pf.fail())RunException(met,"File writing failure.",file);
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Loads data from a file.
//==============================================================================
void JVarsAscii::LoadFile(const std::string &file){
  const char met[]="LoadFile";
  Reset();
  char caderr[512];
  ifstream pf;
  pf.open(file.c_str());
  if(pf){
    while(!pf.eof()){
      char buff[1024];
      pf.getline(buff,1024);
      string tx=buff;
      if(tx!=""){
        int pos=int(tx.find("#"));
        if(pos>0)tx=tx.substr(0,pos);
        tx=StrTrim(tx);
        if(tx!=""){
          pos=int(tx.find("="));
          if(pos<0){ sprintf(caderr,"invalid line:[%s].",tx.c_str()); RunException(met,caderr); }
          AddValueString(StrTrim(tx.substr(0,pos)),StrTrim(tx.substr(pos+1)));
        }
      }
    } 
    if(!pf.eof()&&pf.fail())RunException(met,"Failure reading data from file.",file);
    pf.close();
  }
  else RunException(met,"Cannot open the file.",file);
}

//==============================================================================
/// Indicates whether the value of a requested variable exists or not.
//==============================================================================
int JVarsAscii::FindName(const std::string &name)const{
  int pos=-1;
  string name2=StrUpper(name);
  for(unsigned c=0;c<Count&&pos<0;c++)if(name2==StrUpper(Names[c]))pos=int(c);
  return(pos);
}

//==============================================================================
/// Returns the number of values of a variable.
//==============================================================================
unsigned JVarsAscii::ValuesCount(unsigned pos)const{
  unsigned vals=0;
  if(pos<Count){
    string tx=Values[pos];
    while(tx.length()>0){
      int pos=int(tx.find(","));
      tx=(pos>=0? tx.substr(pos+1): "");
      vals++;
    }
  }
  return(vals);
}

//==============================================================================
/// Returns the requested value as string.
//==============================================================================
std::string JVarsAscii::ValueString(const std::string &name,unsigned value)const{
  const char met[]="GetValueString";
  string texvar="The variable \'";
  texvar=texvar+name+"\' ";
  int pos=FindName(name);
  if(pos<0)RunException(met,texvar+"does not exist.");
  string val;
  unsigned vals=0;
  string tx=Values[pos];
  while(tx.length()>0){
    int pos=int(tx.find(","));
    if(vals==value)val=(pos>=0? tx.substr(0,pos): tx);
    tx=(pos>=0? tx.substr(pos+1): "");
    vals++;
  }
  if(value>=vals)RunException(met,texvar+"there is not the requested value.");
  return(val);
}

//==============================================================================
/// Returns the value of a requested variable.
//==============================================================================
int JVarsAscii::ValueInt(const std::string &name,unsigned value)const{
  string v=ValueString(name,value);
  return(atoi(v.c_str()));
}
//==============================================================================
double JVarsAscii::ValueDouble(const std::string &name,unsigned value)const{
  string v=ValueString(name,value);
  return(atof(v.c_str()));
}

//==============================================================================
/// Adds a new value and returns its position.
//==============================================================================
unsigned JVarsAscii::AddValue(const std::string &name){
  const char met[]="AddValue";
  int pos=FindName(name);
  if(pos<0){
    if(Count>=MAXVARS)RunException(met,"Error can not create more variables.");
    Names[Count]=StrTrim(name);
    Values[Count]="";
    pos=Count;
    Count++;
  }
  return(unsigned(pos));
}

//==============================================================================
/// Adds new variables.
//==============================================================================
unsigned JVarsAscii::AddValueString(const std::string &name,const std::string &v1,const std::string &v2,const std::string &v3){
  unsigned pos=AddValue(name);
  Values[pos]=v1;
  if(!v2.empty())Values[pos]=Values[pos]+","+v2;
  if(!v3.empty())Values[pos]=Values[pos]+","+v3;
  return(pos);
}
//==============================================================================
unsigned JVarsAscii::AddValueInt(const std::string &name,int v1){
  char buff[1024]; sprintf(buff,"%d",v1);
  return(AddValueString(name,buff));
}
//==============================================================================
unsigned JVarsAscii::AddValueInt(const std::string &name,int v1,int v2){
  char buff[1024]; sprintf(buff,"%d,%d",v1,v2);
  return(AddValueString(name,buff));
}
//==============================================================================
unsigned JVarsAscii::AddValueInt(const std::string &name,int v1,int v2,int v3){
  char buff[1024]; sprintf(buff,"%d,%d,%d",v1,v2,v3);
  return(AddValueString(name,buff));
}
//==============================================================================
unsigned JVarsAscii::AddValueDouble(const std::string &name,double v1,const char fmt[]){
  char buff[1024]; sprintf(buff,fmt,v1);
  return(AddValueString(name,buff));
}
//==============================================================================
unsigned JVarsAscii::AddValueDouble(const std::string &name,double v1,double v2,const char fmt[]){
  string v; char buff[1024]; 
  sprintf(buff,fmt,v1); v=buff;
  sprintf(buff,fmt,v2); v=v+","+buff;
  return(AddValueString(name,v));
}
//==============================================================================
unsigned JVarsAscii::AddValueDouble(const std::string &name,double v1,double v2,double v3,const char fmt[]){
  string v; char buff[1024]; 
  sprintf(buff,fmt,v1); v=buff;
  sprintf(buff,fmt,v2); v=v+","+buff;
  sprintf(buff,fmt,v3); v=v+","+buff;
  return(AddValueString(name,v));
}




