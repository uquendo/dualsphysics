#pragma once

#include "Types.h"
#include "TypesDef.h"
#include "JObject.h"

class JProbe : protected JObject
{
        friend class JSph;
private:
        unsigned Nprobes;      ///<Number of probe particles.
        tfloat3 *ProbePos;         ///<Position of the probe particles (X,Y,Z).
        tfloat3 *ProbeVel;         ///<Velocity of the probe particles (X,Y,Z).
        float *ProbeRhop;          ///<Density of the probe particles.
public:
        JProbe(unsigned np=0, tfloat3* Pos = NULL, tfloat3* Vel = NULL, float* Rhop = NULL);
        ~JProbe();

        void LoadCastNodesCountG3D(std::string filename);
        void LoadFileG3D(std::string filename);
        unsigned GetNProbe(){ return Nprobes; };
        void SetProbeArrays( tfloat3* Pos = NULL, tfloat3* Vel = NULL, float* Rhop = NULL){ ProbePos=Pos; ProbeVel=Vel; ProbeRhop=Rhop; };
};

