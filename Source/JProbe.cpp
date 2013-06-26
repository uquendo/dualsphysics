#include "JProbe.h"
#include <fstream>
#include <iostream>
#include <cstring>

JProbe::JProbe(unsigned np, tfloat3* Pos, tfloat3* Vel, float* Rhop) : Nprobes(np), ProbePos(Pos), ProbeVel(Vel), ProbeRhop(Rhop)
{
}

JProbe::~JProbe() 
{
}

void JProbe::LoadFileG3D(std::string filename)
{
        const char met[]="LoadFileG3D";
        std::ifstream pf;
        pf.open(filename.c_str());
        if(pf) {
                pf.ignore(256,'\n');
                float fpx,fpy,fpz; unsigned n;
                for(unsigned i=0;i<Nprobes;i++){
                        pf >> fpx >> fpy >> fpz >> n;
                        ProbePos[i]=TFloat3(fpx,fpy,fpz);
                }
                if(pf.fail()) RunException(met,"Error reading data file.",filename);
                pf.close();
        } else RunException(met,"Cannot open the file.",filename);
        memset(ProbeVel, 0, Nprobes*sizeof(tfloat3));
        memset(ProbeRhop, 0, Nprobes*sizeof(float));
}


void JProbe::LoadCastNodesCountG3D(std::string filename){
        const char met[]="LoadCastNodesCountG3D";
        std::ifstream pf;
        pf.open(filename.c_str());
        if(pf) {
                pf >> Nprobes;
                if(pf.fail()) RunException(met,"Error reading data file.",filename);
                pf.close();
        } else RunException(met,"Cannot open the file.",filename);
}

