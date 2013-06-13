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

/// \file CudaSphFC.cu \brief Implements all the functions to compute forces on GPU.


//##############################################################################
//# Kernels for force computation (including Floating objects).
//##############################################################################
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction of a particle with a set of particles (Fluid-Fluid & Fluid-Bound).  
//------------------------------------------------------------------------------
template<TpKernel tkernel,bool floating,bool xsph> __device__ void KerComputeForcesFluidBox
  (unsigned p1,const unsigned &pini,const unsigned &pfin,float4 *pospres,float4 *velrhop,unsigned *idp,float massp2,float massftp1,float3 posp1,float3 velp1,float3 devsp1,float rhopp1,float3 &acep1,float3 &vcor,float &arp1,float &visc)
{
  for(int p2=pini;p2<pfin;p2++)if(p1!=p2){
    float4 pospres2=pospres[p2];
    float drx=posp1.x-pospres2.x;
    float dry=posp1.y-pospres2.y;
    float drz=posp1.z-pospres2.z;
    float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.fourh2&&rr2>=1e-18f){
      const float4 velrhop2=velrhop[p2];
      const float prrhop2=pospres2.w/(velrhop2.w*velrhop2.w);
      float prs=devsp1.x+prrhop2;
      float wab,frx,fry,frz;
      {//===== Kernel =====
        const float rad=sqrt(rr2);
        const float qq=rad/CTE.h;
        float fac;
        if(tkernel==KERNEL_Cubic){//-Cubic kernel.
          const bool radgt=qq>1;
          const float wqq1=(radgt? 2.0f-qq: qq);
          const float wqq2=wqq1*wqq1;
          const float wqq3=wqq2*wqq1;
          wab=(radgt? CTE.cubic_a24*wqq3: CTE.cubic_a2*(1.0f-1.5f*wqq2+0.75f*wqq3));  
          fac=(radgt? CTE.cubic_c2*wqq2: (CTE.cubic_c1*qq+CTE.cubic_d1*wqq2))/rad;
          //-Tensile correction.
          float fab=wab*CTE.cubic_odwdeltap;
          fab*=fab; fab*=fab; //fab=fab^4
          prs+=fab*(devsp1.y+ prrhop2*(pospres2.w>0? 0.01f: -0.2f) );
        }
        if(tkernel==KERNEL_Wendland){//-Wendland kernel.
          const float wqq=2.f*qq+1.f;
          const float wqq1=1.f-0.5f*qq;
          const float wqq2=wqq1*wqq1;
          wab=CTE.wendland_awen*wqq*wqq2*wqq2;
          fac=CTE.wendland_bwen*qq*wqq2*wqq1/rad;
        } 
        frx=fac*drx; fry=fac*dry; frz=fac*drz;
      }
      if(floating)massp2=(p2<CTE.nbound||idp[p2]<CTE.nbound? CTE.massb: CTE.massf);
      {//===== Aceleration ===== 
        const float p_vpm=-prs*(floating? massp2*massftp1: massp2);
        acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
      }
      float dvx=velp1.x-velrhop2.x, dvy=velp1.y-velrhop2.y, dvz=velp1.z-velrhop2.z;
      float robar=(rhopp1+velrhop2.w)*0.5f;
      {//===== Viscosity ===== 
        const float dot=drx*dvx + dry*dvy + drz*dvz;
        const float dot_rr2=dot/(rr2+CTE.eta2);
        {//-Artificial viscosity 
          if(dot<0){ 
            const float csound=velrhop2.w*OVERRHOPZERO;                        //const float csound=CTE.cs0*powf(rrhop*OVERRHOPZERO,3);          
            const float cbar=(devsp1.z+ CTE.cs0*(csound*csound*csound) )*0.5f;
            const float amubar=CTE.h*dot_rr2;                                  //float amubar=CTE.h*dot/(rr2+CTE.eta2);
            const float pi_visc=(-CTE.visco*cbar*amubar/robar)*massp2;
            acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
          }
        }
        visc=max(dot_rr2,visc);  //reduction of viscdt will be performed on GPU.
      }
      //-Density derivative.
      arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz);
      //-XSPH correction.
      if(xsph){
        const float wab_rhobar=massp2*(wab/robar);
        vcor.x-=wab_rhobar * dvx;
        vcor.y-=wab_rhobar * dvy;
        vcor.z-=wab_rhobar * dvz;
      }
    }
  }
}

//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction of a particle with a set of particles (Bound-Fluid).
//------------------------------------------------------------------------------
template<TpKernel tkernel,bool floating> __device__ void KerComputeForcesBoundBox
  (unsigned p1,const unsigned &pini,const unsigned &pfin,float4 *pospres,float4 *velrhop,unsigned* idp,float massf,float3 posp1,float3 velp1,float &arp1,float &visc)
{
  for(int p2=pini;p2<pfin;p2++)if(p1!=p2){
    float4 pospres2=pospres[p2];
    float drx=posp1.x-pospres2.x;
    float dry=posp1.y-pospres2.y;
    float drz=posp1.z-pospres2.z;
    float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.fourh2&&rr2>=1e-18f){
      const float4 velrhop2=velrhop[p2];
      float frx,fry,frz;
      {//===== Kernel =====
        const float rad=sqrt(rr2);
        const float qq=rad/CTE.h;
        float fac;
        if(tkernel==KERNEL_Cubic){//-Cubic kernel.
          const bool radgt=qq>1;
          float wqq2=(radgt? 2.0f-qq: qq); wqq2*=wqq2;
          fac=(radgt? CTE.cubic_c2*wqq2: (CTE.cubic_c1*qq+CTE.cubic_d1*wqq2))/rad;
        }
        if(tkernel==KERNEL_Wendland){//-Wendland kernel.
          const float wqq1=1.f-0.5f*qq;
          fac=CTE.wendland_bwen*qq*wqq1*wqq1*wqq1/rad;
        } 
        frx=fac*drx; fry=fac*dry; frz=fac*drz;
      }
      float dvx=velp1.x-velrhop2.x, dvy=velp1.y-velrhop2.y, dvz=velp1.z-velrhop2.z;
      {//===== Viscosity ===== 
        const float dot=drx*dvx + dry*dvy + drz*dvz;
        const float dot_rr2=dot/(rr2+CTE.eta2);
        visc=max(dot_rr2,visc);   // reduction of viscdt will be performed on GPU.
      }
      //-Density derivative.
      if(floating)arp1+=(idp[p2]<CTE.nbound? CTE.massb: CTE.massf)*(dvx*frx+dvy*fry+dvz*frz);
      else arp1+=massf*(dvx*frx+dvy*fry+dvz*frz);
    }
  }
}
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction between particles using NeigsCell (Fluid-Fluid & Fluid-Bound).
//------------------------------------------------------------------------------
template<TpKernel tkernel,bool floating,bool xsph,unsigned ncellnv> __global__ void KerComputeForcesFluidNeigs
  (unsigned pini,unsigned n,unsigned ncells,unsigned *cell,uint2 *cellnvb,uint2 *cellnvf,float4 *pospres,float4 *velrhop,unsigned *idp,float3 *ace,float *ar,float3 *velxcor,float *viscdt)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    unsigned p1=p+pini; 
    float arp1=0,visc=0;
    float3 acep1=make_float3(0,0,0);
    float3 vcor;
    if(xsph)vcor=acep1;
    float4 r=velrhop[p1];
    float3 velp1=make_float3(r.x,r.y,r.z);
    float rhopp1=r.w;
    r=pospres[p1];
    float3 posp1=make_float3(r.x,r.y,r.z);
    const float csoun=rhopp1*OVERRHOPZERO;
    float3 devsp1; devsp1.x=r.w/(rhopp1*rhopp1); devsp1.y=devsp1.x*(r.w>0? 0.01f: -0.2f); devsp1.z=CTE.cs0*(csoun*csoun*csoun);
    float massftp1=(floating? (idp[p1]<CTE.nbound? CTE.massb: 1.f): 0);
    const int cel=cell[p];
    //-Interaction with fluid particles.
    uint2 rg;
    for(int r=0;r<ncellnv;r++){
      rg=(!r||rg.y? cellnvf[r*ncells+cel]: make_uint2(0,0));
      if(rg.y)KerComputeForcesFluidBox<tkernel,floating,xsph> (p1,rg.x,rg.y,pospres,velrhop,idp,(floating? 0: CTE.massf),massftp1,posp1,velp1,devsp1,rhopp1,acep1,vcor,arp1,visc);
    }
    //-Interaction with boundaries.
    for(int r=0;r<ncellnv;r++){
      rg=(!r||rg.y? cellnvb[r*ncells+cel]: make_uint2(0,0));
      if(rg.y)KerComputeForcesFluidBox<tkernel,floating,xsph> (p1,rg.x,rg.y,pospres,velrhop,idp,(floating? 0: CTE.massb),massftp1,posp1,velp1,devsp1,rhopp1,acep1,vcor,arp1,visc);
    }
    //-Stores results.
    if(arp1||acep1.x||acep1.y||acep1.z||visc){
      ar[p1]+=arp1;
      float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
      if(xsph){
        r=velxcor[p1]; r.x+=vcor.x; r.y+=vcor.y; r.z+=vcor.z; velxcor[p1]=r;
      }
      if(visc>viscdt[p1])viscdt[p1]=visc;
    }
  }
}

//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction between particles using NeigsCell (Bound-Fluid).
//------------------------------------------------------------------------------
template<TpKernel tkernel,bool floating,unsigned ncellnv> __global__ void KerComputeForcesBoundNeigs
  (unsigned n,unsigned ncells,unsigned *cell,uint2 *cellnvf,float4 *pospres,float4 *velrhop,unsigned *idp,float *ar,float *viscdt)
{
  unsigned p1=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p1<n){
    float arp1=0,visc=0;
    float4 r=velrhop[p1];
    float3 velp1=make_float3(r.x,r.y,r.z);
    r=pospres[p1];
    float3 posp1=make_float3(r.x,r.y,r.z);
    const int cel=cell[p1];
    //-Interaction of boundary with fluid particles.
    uint2 rg;
    for(int r=0;r<ncellnv;r++){
      rg=(!r||rg.y? cellnvf[r*ncells+cel]: make_uint2(0,0));
      if(rg.y)KerComputeForcesBoundBox<tkernel,floating>(p1,rg.x,rg.y,pospres,velrhop,idp,(floating? 0: CTE.massf),posp1,velp1,arp1,visc);
    }
    //-Stores results.
    if(arp1||visc){
      ar[p1]+=arp1;
      if(visc>viscdt[p1])viscdt[p1]=visc;
    }
  }
}

//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction between particles (Fluid-Fluid & Fluid-Bound).
//------------------------------------------------------------------------------
template<TpKernel tkernel,bool floating,bool xsph,unsigned hdiv> __global__ void KerComputeForcesFluid
  (unsigned pini,unsigned n,unsigned ncx,unsigned ncy,unsigned ncz,unsigned *cell,int2 *cellbegb,int2 *cellbegf,float4 *pospres,float4 *velrhop,unsigned *idp,float3 *ace,float *ar,float3 *velxcor,float *viscdt)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    unsigned p1=p+pini; 
    float arp1=0,visc=0;
    float3 acep1=make_float3(0,0,0);
    float3 vcor;
    if(xsph)vcor=acep1;
    float4 r=velrhop[p1];
    float3 velp1=make_float3(r.x,r.y,r.z);
    float rhopp1=r.w;
    r=pospres[p1];
    float3 posp1=make_float3(r.x,r.y,r.z);
    const float csoun=rhopp1*OVERRHOPZERO;
    float3 devsp1; devsp1.x=r.w/(rhopp1*rhopp1); devsp1.y=devsp1.x*(r.w>0? 0.01f: -0.2f); devsp1.z=CTE.cs0*(csoun*csoun*csoun);
    float massftp1=(floating? (idp[p1]<CTE.nbound? CTE.massb: 1.f): 0);
    const int cel=cell[p];
    //-Obtains limits of interaction.
    int cx=cel%ncx;
    int cz=int(cel/(ncx*ncy));
    int cy=int((cel%(ncx*ncy))/ncx);
    int cxini=cx-min(cx,hdiv);
    int cxfin=cx+min(ncx-cx-1,hdiv)+1;
    int yini=cy-min(cy,hdiv);
    int yfin=cy+min(ncy-cy-1,hdiv)+1;
    int zini=cz-min(cz,hdiv);
    int zfin=cz+min(ncz-cz-1,hdiv)+1;
    //-Interaction with fluid particles.
    for(int z=zini;z<zfin;z++){
      int zmod=(ncx*ncy)*z;
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+ncx*y;
        unsigned pini,pfin=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=cellbegf[x+ymod];
          if(cbeg.y){
            if(!pfin)pini=cbeg.x;
            pfin=cbeg.y;
          }
        }
        if(pfin)KerComputeForcesFluidBox<tkernel,floating,xsph> (p1,pini,pfin,pospres,velrhop,idp,(floating? 0: CTE.massf),massftp1,posp1,velp1,devsp1,rhopp1,acep1,vcor,arp1,visc);
      }
    }
    //-Interaction with boundaries.
    for(int z=zini;z<zfin;z++){
      int zmod=(ncx*ncy)*z;
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+ncx*y;
        unsigned pini,pfin=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=cellbegb[x+ymod];
          if(cbeg.y){
            if(!pfin)pini=cbeg.x;
            pfin=cbeg.y;
          }
        }
        if(pfin)KerComputeForcesFluidBox<tkernel,floating,xsph> (p1,pini,pfin,pospres,velrhop,idp,(floating? 0: CTE.massb),massftp1,posp1,velp1,devsp1,rhopp1,acep1,vcor,arp1,visc);
      }
    }
    //-Stores results.
    if(arp1||acep1.x||acep1.y||acep1.z||visc){
      ar[p1]+=arp1;
      float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
      if(xsph){
        r=velxcor[p1]; r.x+=vcor.x; r.y+=vcor.y; r.z+=vcor.z; velxcor[p1]=r;
      }
      if(visc>viscdt[p1])viscdt[p1]=visc;
    }
  }
}

//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction between particles (Bound-Fluid).
//------------------------------------------------------------------------------
template<TpKernel tkernel,bool floating,unsigned hdiv> __global__ void KerComputeForcesBound(unsigned n,unsigned ncx,unsigned ncy,unsigned ncz,unsigned *cell,int2 *cellbegf,float4 *pospres,float4 *velrhop,unsigned *idp,float *ar,float *viscdt)
{
  unsigned p1=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p1<n){
    float arp1=0,visc=0;
    float4 r=velrhop[p1];
    float3 velp1=make_float3(r.x,r.y,r.z);
    r=pospres[p1];
    float3 posp1=make_float3(r.x,r.y,r.z);
    const int cel=cell[p1];
    //-Obtains limits of interaction.
    int cx=cel%ncx;
    int cz=int(cel/(ncx*ncy));
    int cy=int((cel%(ncx*ncy))/ncx);
    int cxini=cx-min(cx,hdiv);
    int cxfin=cx+min(ncx-cx-1,hdiv)+1;
    int yini=cy-min(cy,hdiv);
    int yfin=cy+min(ncy-cy-1,hdiv)+1;
    int zini=cz-min(cz,hdiv);
    int zfin=cz+min(ncz-cz-1,hdiv)+1;
    //-Interaction of boundary with fluid particles.
    for(int z=zini;z<zfin;z++){
      int zmod=(ncx*ncy)*z;
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+ncx*y;
        unsigned pini,pfin=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=cellbegf[x+ymod];
          if(cbeg.y){
            if(!pfin)pini=cbeg.x;
            pfin=cbeg.y;
          }
        }
        if(pfin)KerComputeForcesBoundBox<tkernel,floating>(p1,pini,pfin,pospres,velrhop,idp,(floating? 0: CTE.massf),posp1,velp1,arp1,visc);
      }
    }
    //-Stores results.
    if(arp1||visc){
      ar[p1]+=arp1;
      if(visc>viscdt[p1])viscdt[p1]=visc;
    }
  }
}

//##############################################################################
//# Kernels for force computation when using Laminar+SPS viscosity and KGC. 
//##############################################################################
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction of a particle with a set of particles (Fluid-Fluid & Fluid-Bound). 
//------------------------------------------------------------------------------
template<TpKernel tkernel,TpVisco tvisco,bool kgc,bool xsph> __device__ void KerComputeForcesFullFluidBox
  (unsigned p1,const unsigned &pini,const unsigned &pfin,float4 *pospres,float4 *velrhop,unsigned *idp,const tsymatrix3f *tau,float massp2,float3 posp1,float3 velp1,float3 devsp1,float rhopp1,const tsymatrix3f &matkgcp1,const tsymatrix3f &taup1,tsymatrix3f &csphp1,float3 &acep1,float3 &vcor,float &arp1,float &visc)
{
  for(int p2=pini;p2<pfin;p2++)if(p1!=p2){
    float4 pospres2=pospres[p2];
    float drx=posp1.x-pospres2.x;
    float dry=posp1.y-pospres2.y;
    float drz=posp1.z-pospres2.z;
    float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.fourh2&&rr2>=1e-18f){
      const float4 velrhop2=velrhop[p2];
      const float prrhop2=pospres2.w/(velrhop2.w*velrhop2.w);
      float prs=devsp1.x+prrhop2;
      float wab,frx,fry,frz;
      {//===== Kernel =====
        const float rad=sqrt(rr2);
        const float qq=rad/CTE.h;
        float fac;
        if(tkernel==KERNEL_Cubic){//-Cubic kernel.
          const bool radgt=qq>1;
          const float wqq1=(radgt? 2.0f-qq: qq);
          const float wqq2=wqq1*wqq1;
          const float wqq3=wqq2*wqq1;
          wab=(radgt? CTE.cubic_a24*wqq3: CTE.cubic_a2*(1.0f-1.5f*wqq2+0.75f*wqq3)); 
          fac=(radgt? CTE.cubic_c2*wqq2: (CTE.cubic_c1*qq+CTE.cubic_d1*wqq2))/rad;
          //-Tensile correction.
          float fab=wab*CTE.cubic_odwdeltap;
          fab*=fab; fab*=fab; //fab=fab^4
          prs+=fab*(devsp1.y+ prrhop2*(pospres2.w>0? 0.01f: -0.2f) );
        }
        if(tkernel==KERNEL_Wendland){//-Wendland kernel.
          const float wqq=2.f*qq+1.f;
          const float wqq1=1.f-0.5f*qq;
          const float wqq2=wqq1*wqq1;
          wab=CTE.wendland_awen*wqq*wqq2*wqq2;
          fac=CTE.wendland_bwen*qq*wqq2*wqq1/rad;
        } 
        frx=fac*drx; fry=fac*dry; frz=fac*drz;
      }
      //===== Kernel Gradient Correction ===== 
      if(kgc==true){
        float frx2=matkgcp1.xx*frx + matkgcp1.xy*fry + matkgcp1.xz*frz;      
        float fry2=matkgcp1.xy*frx + matkgcp1.yy*fry + matkgcp1.yz*frz;
        float frz2=matkgcp1.xz*frx + matkgcp1.yz*fry + matkgcp1.zz*frz;
        frx=frx2; fry=fry2; frz=frz2;
      }
      {//===== Aceleration ===== 
        const float p_vpm=-prs*massp2;
        acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
      }
      float dvx=velp1.x-velrhop2.x, dvy=velp1.y-velrhop2.y, dvz=velp1.z-velrhop2.z;
      float robar=(rhopp1+velrhop2.w)*0.5f;
      {//===== Viscosity ===== 
        const float dot=drx*dvx + dry*dvy + drz*dvz;
        const float dot_rr2=dot/(rr2+CTE.eta2);
        if(tvisco==VISCO_Artificial){//-Artificial viscosity. 
          if(dot<0){
            const float csound=velrhop2.w*OVERRHOPZERO;                        //const float csound=CTE.cs0*powf(rrhop*OVERRHOPZERO,3); 
            const float cbar=(devsp1.z+ CTE.cs0*(csound*csound*csound) )*0.5f;
            const float amubar=CTE.h*dot_rr2;                                  //float amubar=CTE.h*dot/(rr2+CTE.eta2);
            const float pi_visc=(-CTE.visco*cbar*amubar/robar)*massp2;
            acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
          }
        }
        if(tvisco==VISCO_LaminarSPS){//-LaminarSPS viscosity. 
          const float temp=2.0f*CTE.visco/((rr2+CTE.eta2)*robar);
          const float vtemp=massp2*temp*(drx*frx+dry*fry+drz*frz);  
          acep1.x+=vtemp*dvx; acep1.y+=vtemp*dvy; acep1.z+=vtemp*dvz;
          // SPS turbulence model.
          tsymatrix3f tausum=taup1;
          if(p2>CTE.nbound){
            tausum=tau[p2-CTE.nbound];
            tausum.xx+=taup1.xx;
            tausum.xy+=taup1.xy;
            tausum.xz+=taup1.xz;
            tausum.yy+=taup1.yy;
            tausum.yz+=taup1.yz;
            tausum.zz+=taup1.zz;
          }
          acep1.x+=massp2*(tausum.xx*frx+tausum.xy*fry+tausum.xz*frz);
          acep1.y+=massp2*(tausum.xy*frx+tausum.yy*fry+tausum.yz*frz);
          acep1.z+=massp2*(tausum.xz*frx+tausum.yz*fry+tausum.zz*frz);
          // CSPH terms.
          const float volp2=-massp2/velrhop2.w;
          float dv=dvx*volp2; csphp1.xx+=dv*frx; csphp1.xy+=dv*fry; csphp1.xz+=dv*frz;
                dv=dvy*volp2; csphp1.xy+=dv*frx; csphp1.yy+=dv*fry; csphp1.yz+=dv*frz;
                dv=dvz*volp2; csphp1.xz+=dv*frx; csphp1.yz+=dv*fry; csphp1.zz+=dv*frz;
        }
        visc=max(dot_rr2,visc);   // reduction of viscdt will be performed on GPU.
      }
      //-Density derivative.
      arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz);
      //-XSPH correction.
      if(xsph){
        const float wab_rhobar=massp2*(wab/robar);
        vcor.x-=wab_rhobar * dvx;
        vcor.y-=wab_rhobar * dvy;
        vcor.z-=wab_rhobar * dvz;
      }
    }
  }
}
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction between particles using NeigsCell (Fluid-Fluid & Fluid-Bound).
//------------------------------------------------------------------------------
template<TpKernel tkernel,TpVisco tvisco,bool kgc,bool xsph,unsigned ncellnv> __global__ void KerComputeForcesFullFluidNeigs
  (unsigned pini,unsigned n,unsigned ncells,unsigned *cell,uint2 *cellnvb,uint2 *cellnvf,float4 *pospres,float4 *velrhop,unsigned *idp,float3 *ace,float *ar,float3 *velxcor,float *viscdt,const tsymatrix3f *matkgc,const tsymatrix3f *tau,tsymatrix3f *csph)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    unsigned p1=p+pini;      
    float arp1=0,visc=0;
    float3 acep1=make_float3(0,0,0);
    float3 vcor;
    if(xsph)vcor=acep1;
    float4 r=velrhop[p1];
    float3 velp1=make_float3(r.x,r.y,r.z);
    float rhopp1=r.w;
    r=pospres[p1];
    float3 posp1=make_float3(r.x,r.y,r.z);
    const float csoun=rhopp1*OVERRHOPZERO;
    float3 devsp1; devsp1.x=r.w/(rhopp1*rhopp1); devsp1.y=devsp1.x*(r.w>0? 0.01f: -0.2f); devsp1.z=CTE.cs0*(csoun*csoun*csoun);
    tsymatrix3f matkgcp1=matkgc[p];
    tsymatrix3f taup1=tau[p];
    tsymatrix3f csphp1={0,0,0,0,0,0};
    const int cel=cell[p];
    //-Interaction with fluid particles.
    uint2 rg;
    for(int r=0;r<ncellnv;r++){
      rg=(!r||rg.y? cellnvf[r*ncells+cel]: make_uint2(0,0));
      if(rg.y)KerComputeForcesFullFluidBox<tkernel,tvisco,kgc,xsph> (p1,rg.x,rg.y,pospres,velrhop,idp,tau,CTE.massf,posp1,velp1,devsp1,rhopp1,matkgcp1,taup1,csphp1,acep1,vcor,arp1,visc);
    }
    //-Interaction with boundaries.
    for(int r=0;r<ncellnv;r++){
      rg=(!r||rg.y? cellnvb[r*ncells+cel]: make_uint2(0,0));
      if(rg.y)KerComputeForcesFullFluidBox<tkernel,tvisco,kgc,xsph> (p1,rg.x,rg.y,pospres,velrhop,idp,tau,CTE.massb,posp1,velp1,devsp1,rhopp1,matkgcp1,taup1,csphp1,acep1,vcor,arp1,visc);
    }
    //-Stores results.
    if(arp1||acep1.x||acep1.y||acep1.z||visc){
      ar[p1]+=arp1;
      float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
      if(xsph){ r=velxcor[p1]; r.x+=vcor.x; r.y+=vcor.y; r.z+=vcor.z; velxcor[p1]=r; }
      if(visc>viscdt[p1])viscdt[p1]=visc;
      if(tvisco==VISCO_LaminarSPS)csph[p]=csphp1;
    }
  }
}
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction between particles (Fluid-Fluid & Fluid-Bound).
//------------------------------------------------------------------------------
template<TpKernel tkernel,TpVisco tvisco,bool kgc,bool xsph,unsigned hdiv> __global__ void KerComputeForcesFullFluid
  (unsigned pini,unsigned n,unsigned ncx,unsigned ncy,unsigned ncz,unsigned *cell,int2 *cellbegb,int2 *cellbegf,float4 *pospres,float4 *velrhop,unsigned *idp,float3 *ace,float *ar,float3 *velxcor,float *viscdt,const tsymatrix3f *matkgc,const tsymatrix3f *tau,tsymatrix3f *csph)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    unsigned p1=p+pini; 
    float arp1=0,visc=0;
    float3 acep1=make_float3(0,0,0);
    float3 vcor;
    if(xsph)vcor=acep1;
    float4 r=velrhop[p1];
    float3 velp1=make_float3(r.x,r.y,r.z);
    float rhopp1=r.w;
    r=pospres[p1];
    float3 posp1=make_float3(r.x,r.y,r.z);
    const float csoun=rhopp1*OVERRHOPZERO;
    float3 devsp1; devsp1.x=r.w/(rhopp1*rhopp1); devsp1.y=devsp1.x*(r.w>0? 0.01f: -0.2f); devsp1.z=CTE.cs0*(csoun*csoun*csoun);
    tsymatrix3f matkgcp1=matkgc[p];
    tsymatrix3f taup1=tau[p];
    tsymatrix3f csphp1={0,0,0,0,0,0};
    const int cel=cell[p];
    //-Obtains limits of interaction.
    int cx=cel%ncx;
    int cz=int(cel/(ncx*ncy));
    int cy=int((cel%(ncx*ncy))/ncx);
    int cxini=cx-min(cx,hdiv);
    int cxfin=cx+min(ncx-cx-1,hdiv)+1;
    int yini=cy-min(cy,hdiv);
    int yfin=cy+min(ncy-cy-1,hdiv)+1;
    int zini=cz-min(cz,hdiv);
    int zfin=cz+min(ncz-cz-1,hdiv)+1;
    //-Interaction with fluid particles.
    for(int z=zini;z<zfin;z++){
      int zmod=(ncx*ncy)*z;
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+ncx*y;
        unsigned pini,pfin=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=cellbegf[x+ymod];
          if(cbeg.y){
            if(!pfin)pini=cbeg.x;
            pfin=cbeg.y;
          }
        }
        if(pfin)KerComputeForcesFullFluidBox<tkernel,tvisco,kgc,xsph> (p1,pini,pfin,pospres,velrhop,idp,tau,CTE.massf,posp1,velp1,devsp1,rhopp1,matkgcp1,taup1,csphp1,acep1,vcor,arp1,visc);
      }
    }
    //-Interaction with boundaries.
    for(int z=zini;z<zfin;z++){
      int zmod=(ncx*ncy)*z;
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+ncx*y;
        unsigned pini,pfin=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=cellbegb[x+ymod];
          if(cbeg.y){
            if(!pfin)pini=cbeg.x;
            pfin=cbeg.y;
          }
        }
        if(pfin)KerComputeForcesFullFluidBox<tkernel,tvisco,kgc,xsph> (p1,pini,pfin,pospres,velrhop,idp,tau,CTE.massb,posp1,velp1,devsp1,rhopp1,matkgcp1,taup1,csphp1,acep1,vcor,arp1,visc);
      }
    }
    //-Stores results.
    if(arp1||acep1.x||acep1.y||acep1.z||visc){
      ar[p1]+=arp1;
      float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
      if(xsph){ r=velxcor[p1]; r.x+=vcor.x; r.y+=vcor.y; r.z+=vcor.z; velxcor[p1]=r; }
      if(visc>viscdt[p1])viscdt[p1]=visc;
      if(tvisco==VISCO_LaminarSPS)csph[p]=csphp1;
    }
  }
}


//##############################################################################
//# Kernel: CELLMODE_Hneigs for Shepard density filter.
//##############################################################################
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction of a particle with a set of particles.
//------------------------------------------------------------------------------
template<TpKernel tkernel> __device__ void KerComputeForcesShepardBox(unsigned p1,const unsigned &pini,const unsigned &pfin,float massfluid,float4 *posrhop,float4 posrhop1,float &fdwabp1,float &fdrhopp1)
{
  for(int p2=pini;p2<pfin;p2++)if(p1!=p2){
    float4 posrhop2=posrhop[p2];
    float drx=posrhop1.x-posrhop2.x;
    float dry=posrhop1.y-posrhop2.y;
    float drz=posrhop1.z-posrhop2.z;
    float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.fourh2&&rr2>=1e-18f){
      float wab;
      {//===== Kernel =====
        const float rad=sqrt(rr2);
        const float qq=rad/CTE.h;
        if(tkernel==KERNEL_Cubic){//-Cubic kernel.
          const bool radgt=qq>1;
          const float wqq1=(radgt? 2.0f-qq: qq);
          const float wqq2=wqq1*wqq1;
          const float wqq3=wqq2*wqq1;
          wab=(radgt? CTE.cubic_a24*wqq3: CTE.cubic_a2*(1.0f-1.5f*wqq2+0.75f*wqq3));  
        }
        if(tkernel==KERNEL_Wendland){//-Wendland kernel.
          float wqq2=1.f-0.5f*qq; wqq2*=wqq2;
          wab=CTE.wendland_awen*(2*qq+1)*wqq2*wqq2;
        }
      }
      fdwabp1+=wab*(massfluid/posrhop2.w);
      fdrhopp1+=wab;
    }
  }
}
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction Shepard between fluid particles.
//------------------------------------------------------------------------------
template<TpKernel tkernel,unsigned ncellnv> __global__ void KerComputeForcesShepardNeigs(unsigned pini,unsigned n,unsigned ncells,unsigned *cell,uint2 *cellnvf,float massfluid,float4 *posrhop,float *fdrhop)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    unsigned p1=p+pini; 
    float4 posrhop1=posrhop[p1];
    float fdrhopp1=CTE.cteshepard;    
    float fdwabp1=fdrhopp1*(massfluid/posrhop1.w);
    const int cel=cell[p];
    //-Interaction with fluid particles.
    uint2 rg;
    for(int r=0;r<ncellnv;r++){
      rg=(!r||rg.y? cellnvf[r*ncells+cel]: make_uint2(0,0));
      if(rg.y)KerComputeForcesShepardBox <tkernel> (p1,rg.x,rg.y,massfluid,posrhop,posrhop1,fdwabp1,fdrhopp1);
    }
    //-Stores results.
    fdrhop[p]=(fdrhopp1*massfluid)/fdwabp1;
  }
}
//##############################################################################
//# Kernel: CELLMODE_2H & CELLMODE_H for Shepard density filter
//##############################################################################
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction between particles 
//------------------------------------------------------------------------------
template<TpKernel tkernel,unsigned hdiv> __global__ void KerComputeForcesShepard(unsigned pini,unsigned n,unsigned ncx,unsigned ncy,unsigned ncz,unsigned *cell,int2 *cellbegf,float massfluid,float4 *posrhop,float *fdrhop)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    unsigned p1=p+pini; 
    float4 posrhop1=posrhop[p1];
    float fdrhopp1=CTE.cteshepard;     
    float fdwabp1=fdrhopp1*(massfluid/posrhop1.w);
    const int cel=cell[p];
    //-Obtains limits of interaction.
    int cx=cel%ncx;
    int cz=int(cel/(ncx*ncy));
    int cy=int((cel%(ncx*ncy))/ncx);
    int cxini=cx-min(cx,hdiv);
    int cxfin=cx+min(ncx-cx-1,hdiv)+1;
    int yini=cy-min(cy,hdiv);
    int yfin=cy+min(ncy-cy-1,hdiv)+1;
    int zini=cz-min(cz,hdiv);
    int zfin=cz+min(ncz-cz-1,hdiv)+1;
    //-Interaction with fluid particles.
    for(int z=zini;z<zfin;z++){
      int zmod=(ncx*ncy)*z;
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+ncx*y;
        unsigned pini,pfin=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=cellbegf[x+ymod];
          if(cbeg.y){
            if(!pfin)pini=cbeg.x;
            pfin=cbeg.y;
          }
        }
        if(pfin)KerComputeForcesShepardBox <tkernel> (p1,pini,pfin,massfluid,posrhop,posrhop1,fdwabp1,fdrhopp1);
      }
    }
    //-Stores results.
    fdrhop[p]=(fdrhopp1*massfluid)/fdwabp1;
  }
}

//------------------------------------------------------------------------------
/// CUDA KERNEL that makes ace[].y equal to zero.
//------------------------------------------------------------------------------
__global__ void KerResetAcey(unsigned n,unsigned npb,float3 *ace)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n)ace[p+npb].y=0;
}



//##############################################################################
//# Kernel: CELLMODE_Hneigs for MatrixKGC when kernel gradient correction.
//##############################################################################
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction of a particle with a set of particles 
//------------------------------------------------------------------------------
template<TpKernel tkernel> __device__ void KerComputeMatrixKgcBox(unsigned p1,const unsigned &pini,const unsigned &pfin,float massp2,const float4 *posrhop,const float3 &posp1,tsymatrix3f &matkgcp1)
{
  for(int p2=pini;p2<pfin;p2++)if(p1!=p2){
    float4 posrhopp2=posrhop[p2];
    float drx=posp1.x-posrhopp2.x;
    float dry=posp1.y-posrhopp2.y;
    float drz=posp1.z-posrhopp2.z;
    float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.fourh2&&rr2>=1e-18f){
      float frx,fry,frz;
      {//===== Kernel =====
        const float rad=sqrt(rr2);
        const float qq=rad/CTE.h;
        float fac;
        if(tkernel==KERNEL_Cubic){//-Cubic kernel.
          const float wqq1=2.0f-qq;
          fac=(qq>1? CTE.cubic_c2*(wqq1*wqq1): (CTE.cubic_c1+CTE.cubic_d1*qq)*qq) /rad;
        }
        if(tkernel==KERNEL_Wendland){//-Wendland kernel.
          const float wqq1=1.f-0.5f*qq;
          fac=CTE.wendland_bwen*qq*wqq1*wqq1*wqq1/rad;
        } 
        frx=fac*drx; fry=fac*dry; frz=fac*drz;
      }
      //===== KGC Matrix  =====
      const float volp2=-massp2/posrhopp2.w;
      float fr=frx*volp2; matkgcp1.xx+=fr*drx; matkgcp1.xy+=fr*dry; matkgcp1.xz+=fr*drz; 
            fr=fry*volp2; matkgcp1.xy+=fr*drx; matkgcp1.yy+=fr*dry; matkgcp1.yz+=fr*drz; 
            fr=frz*volp2; matkgcp1.xz+=fr*drx; matkgcp1.yz+=fr*dry; matkgcp1.zz+=fr*drz; 
    }
  }
}
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction between particles using NeigsCell.
//------------------------------------------------------------------------------
template<TpKernel tkernel,unsigned ncellnv> __global__ void KerComputeMatrixKgcNeigs
  (unsigned pini,unsigned n,unsigned ncells,const unsigned *cell,const uint2 *cellnvb,const uint2 *cellnvf,const float4 *posrhop,float massb,float massf,tsymatrix3f *matkgc)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    unsigned p1=p+pini; 
    tsymatrix3f matkgcp1={0,0,0,0,0,0};
    float4 r=posrhop[p1];
    float3 posp1=make_float3(r.x,r.y,r.z);
    const int cel=cell[p];
    //-Interaction with fluid particles.
    uint2 rg;
    for(int r=0;r<ncellnv;r++){
      rg=(!r||rg.y? cellnvf[r*ncells+cel]: make_uint2(0,0));
      if(rg.y)KerComputeMatrixKgcBox<tkernel> (p1,rg.x,rg.y,massf,posrhop,posp1,matkgcp1);
    }
    //-Interaction with boundaries.
    for(int r=0;r<ncellnv;r++){
      rg=(!r||rg.y? cellnvb[r*ncells+cel]: make_uint2(0,0));
      if(rg.y)KerComputeMatrixKgcBox<tkernel> (p1,rg.x,rg.y,massb,posrhop,posp1,matkgcp1);
    }
    //-Stores results.
    matkgc[p]=matkgcp1;
  }
}
//------------------------------------------------------------------------------
/// CUDA KERNEL that computes interaction between particles (Fluid-Fluid & Fluid-Bound).
//------------------------------------------------------------------------------
template<TpKernel tkernel,unsigned hdiv> __global__ void KerComputeMatrixKgc
  (unsigned pini,unsigned n,unsigned ncx,unsigned ncy,unsigned ncz,const unsigned *cell,const int2 *cellbegb,const int2 *cellbegf,const float4 *posrhop,float massb,float massf,tsymatrix3f *matkgc)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    unsigned p1=p+pini;  
    tsymatrix3f matkgcp1={0,0,0,0,0,0};
    float4 r=posrhop[p1];
    float3 posp1=make_float3(r.x,r.y,r.z);
    const int cel=cell[p];
    //-Obtains limits of interaction.
    int cx=cel%ncx;
    int cz=int(cel/(ncx*ncy));
    int cy=int((cel%(ncx*ncy))/ncx);
    int cxini=cx-min(cx,hdiv);
    int cxfin=cx+min(ncx-cx-1,hdiv)+1;
    int yini=cy-min(cy,hdiv);
    int yfin=cy+min(ncy-cy-1,hdiv)+1;
    int zini=cz-min(cz,hdiv);
    int zfin=cz+min(ncz-cz-1,hdiv)+1;
    //-Interaction with fluid particles.
    for(int z=zini;z<zfin;z++){
      int zmod=(ncx*ncy)*z;
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+ncx*y;
        unsigned pini,pfin=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=cellbegf[x+ymod];
          if(cbeg.y){
            if(!pfin)pini=cbeg.x;
            pfin=cbeg.y;
          }
        }
        if(pfin)KerComputeMatrixKgcBox<tkernel> (p1,pini,pfin,massf,posrhop,posp1,matkgcp1);
      }
    }
    //-Interaction with boundaries.
    for(int z=zini;z<zfin;z++){
      int zmod=(ncx*ncy)*z;
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+ncx*y;
        unsigned pini,pfin=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=cellbegb[x+ymod];
          if(cbeg.y){
            if(!pfin)pini=cbeg.x;
            pfin=cbeg.y;
          }
        }
        if(pfin)KerComputeMatrixKgcBox<tkernel> (p1,pini,pfin,massb,posrhop,posp1,matkgcp1);
      }
    }
    //-Stores results.
    matkgc[p]=matkgcp1;
  }
}

//------------------------------------------------------------------------------
/// CUDA KERNEL that prepares variables for interaction MatrixKgc.
//------------------------------------------------------------------------------
__global__ void KerPreInteraction_MatrixKgc(unsigned np,const float3 *pos,const float *rhop,float4 *posrhop)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if(p<np){
    float3 rpos=pos[p];
    posrhop[p]=make_float4(rpos.x,rpos.y,rpos.z,rhop[p]);
  }
}

//------------------------------------------------------------------------------
/// CUDA KERNEL that inverts MatrixKgc.
//------------------------------------------------------------------------------
template<bool simula2d> __global__ void KerInvertMatrixKgc(unsigned npf,tsymatrix3f *matkgc)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if(p<npf){
    tsymatrix3f mat=matkgc[p];
    //-Matrix in case of 2D simulation.
    if(simula2d)mat.yy=1.0f;
    //-Inversion of a symmetric matrix.
    mat.xy*=0.5f;
    mat.xz*=0.5f;
    mat.yz*=0.5f;
    float det=mat.xx*mat.yy*mat.zz + 2.f*mat.xy*mat.yz*mat.xz - mat.xz*mat.yy*mat.xz - mat.xx*mat.yz*mat.yz - mat.xy*mat.xy*mat.zz; //-Determinant of the matrix.
    //-Matrix inversion.
    tsymatrix3f invmat={1,0,0 ,1,0 ,1}; //-Matrix when no correction will be applied.
    if(fabs(det)>0.01 && fabs(mat.xx)>MATRIXKGC_CONTROL && fabs(mat.yy)>MATRIXKGC_CONTROL && fabs(mat.zz)>MATRIXKGC_CONTROL){
      invmat.xx=(mat.yy*mat.zz-mat.yz*mat.yz)/det;
      invmat.xy=(mat.xz*mat.yz-mat.xy*mat.zz)/det;
      invmat.xz=(mat.xy*mat.yz-mat.yy*mat.xz)/det;
      invmat.yy=(mat.xx*mat.zz-mat.xz*mat.xz)/det;
      invmat.yz=(mat.xy*mat.xz-mat.xx*mat.yz)/det;
      invmat.zz=(mat.xx*mat.yy-mat.xy*mat.xy)/det;
    }
    matkgc[p]=invmat;
  }
}

//==============================================================================
/// Computes matrix of kernel gradient correction (KGC).
//==============================================================================
template<TpKernel tkernel,unsigned ncellnv> void CsInteraction_MatrixKgc(StDeviceContext *dc,StDeviceCte *cte){
  TmgStart(dc->timers,TMG_CfMatrixKgc);
  cudaMemset(dc->matkgc,0,sizeof(tsymatrix3f)*dc->npf);  //matkgc[]=0
  dim3 sgrid=CsGetGridSize(dc->npok,BLOCKSIZE);
  KerPreInteraction_MatrixKgc <<<sgrid,BLOCKSIZE>>>(dc->npok,dc->pos,dc->rhop,dc->pospres);
  switch(dc->cellmode){
    case CELLMODE_Hneigs:{
      //-Interaction Fluid-Fluid & Fluid-Boundary.
      unsigned bsize=dc->bsmatrixkgc;
      dim3 sgridf=CsGetGridSize(dc->npf,bsize);
      KerComputeMatrixKgcNeigs<tkernel,ncellnv> <<<sgridf,bsize>>> (dc->npb,dc->npf,dc->nctotmax-1,dc->cellf,dc->cellnvb,dc->cellnvf,dc->pospres,cte->massb,cte->massf,dc->matkgc);
    }break;
    case CELLMODE_2H:
    case CELLMODE_H:{
      //-Interaction Fluid-Fluid & Fluid-Boundary.
      unsigned bsize=dc->bsmatrixkgc;
      dim3 sgridf=CsGetGridSize(dc->npf,bsize);
      KerComputeMatrixKgc<tkernel,(ncellnv==9? 1: 2)> <<<sgridf,bsize>>> (dc->npb,dc->npf,dc->ncx,dc->ncy,dc->ncz,dc->cellf,dc->cellbegb,dc->cellbegf,dc->pospres,cte->massb,cte->massf,dc->matkgc);
    }break;
    default: throw std::string("CsInteraction_MatrixKgc> CellMode unrecognised.");
  }
  //-Inversion of a matrix matkgc.
  if(dc->simulate2d)KerInvertMatrixKgc<true> <<<sgrid,BLOCKSIZE>>>(dc->npf,dc->matkgc);
  else KerInvertMatrixKgc<false> <<<sgrid,BLOCKSIZE>>>(dc->npf,dc->matkgc);  
  CheckErrorCuda("Failed while executing kernels of MatrixKgc interaction.");
  TmgStop(dc->timers,TMG_CfMatrixKgc);
}



//------------------------------------------------------------------------------
/// CUDA KERNEL that prepares variables for force computation.
//------------------------------------------------------------------------------
__global__ void KerPreInteraction_Forces(unsigned n,unsigned npb,float gamma,float b,float3 gravity,const float3 *pos,const float3 *vel,const float *rhop,float *csound,float3 *ace,float4 *pospres,float4 *velrhop)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    float rrhop=rhop[p];
    float press=b*(powf(rrhop*OVERRHOPZERO,gamma)-1.0f);
    float csoun=rrhop*OVERRHOPZERO;
    csound[p]=CTE.cs0*(csoun*csoun*csoun);
    float3 r=pos[p]; pospres[p]=make_float4(r.x,r.y,r.z,press);
    r=vel[p]; velrhop[p]=make_float4(r.x,r.y,r.z,rrhop);
    ace[p]=(p<npb? make_float3(0,0,0): gravity);
  }
}

//==============================================================================
/// Prepares variables for force computation.
//==============================================================================
void CsPreInteraction_Forces(StDeviceContext *dc,bool xsph){
  TmgStart(dc->timers,TMG_CfPreForces);
  cudaMemset(dc->ar,0,sizeof(float)*dc->npok);  //ar[]=0
  if(xsph)cudaMemset(dc->velxcor,0,sizeof(float3)*dc->npok);  //velxcor[]=0
  if(dc->tvisco==VISCO_LaminarSPS)cudaMemset(dc->csph,0,sizeof(tsymatrix3f)*(dc->npf));  //csph[]={0,..,0}
  cudaMemset(dc->viscdt,0,sizeof(float)*dc->npok);  
  dim3 sgrid=CsGetGridSize(dc->npok,BLOCKSIZE);
  KerPreInteraction_Forces <<<sgrid,BLOCKSIZE>>> (dc->npok,dc->npb,dc->gamma,dc->b,dc->gravity,dc->pos,dc->vel,dc->rhop,dc->csound,dc->ace,dc->pospres,dc->velrhop);
  CheckErrorCuda("Failed preparing vars for interaction forces.");
  TmgStop(dc->timers,TMG_CfPreForces);
}

//------------------------------------------------------------------------------
/// CUDA KERNEL that computes sub-particle stress tensor (Tau) for SPS turbulence model.
//------------------------------------------------------------------------------
__global__ void KerSPSCalcTau(unsigned n,unsigned npb,float smag,float blin,const float *rhop,const tsymatrix3f *csph,tsymatrix3f *tau)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    tsymatrix3f rcsph=csph[p];
    const float pow1=rcsph.xx*rcsph.xx + rcsph.yy*rcsph.yy + rcsph.zz*rcsph.zz;
    const float prr=pow1+pow1 + rcsph.xy*rcsph.xy + rcsph.xz*rcsph.xz + rcsph.yz*rcsph.yz;
    const float visc_SPS=smag*sqrt(prr);
    const float div_u=rcsph.xx+rcsph.yy+rcsph.zz;
    const float sps_k=(2.0f/3.0f)*visc_SPS*div_u;
    const float sps_Blin=blin*prr;
    const float sumsps=-(sps_k+sps_Blin);
    const float twovisc_SPS=(visc_SPS+visc_SPS);
    const float one_rho2 = 1.0f/rhop[p+npb];
    tsymatrix3f rtau;
    rtau.xx=one_rho2*(twovisc_SPS*rcsph.xx +sumsps);
    rtau.xy=one_rho2*(visc_SPS*rcsph.xy);
    rtau.xz=one_rho2*(visc_SPS*rcsph.xz);
    rtau.yy=one_rho2*(twovisc_SPS*rcsph.yy +sumsps);
    rtau.yz=one_rho2*(visc_SPS*rcsph.yz);
    rtau.zz=one_rho2*(twovisc_SPS*rcsph.zz +sumsps);
    tau[p]=rtau;
  }
}

//==============================================================================
/// Interaction between particles.
//==============================================================================
template<bool xsph,TpKernel tkernel,TpVisco tvisco,bool kgc,bool floating,unsigned ncellnv> void CsInteraction_Forces(StDeviceContext *dc,StDeviceCte *cte){
  //-Interaction to compute matrix for Kernel Gradient Correction.
  if(dc->kgc)CsInteraction_MatrixKgc<tkernel,ncellnv>(dc,cte);
  //-Interaction Forces
  CsPreInteraction_Forces(dc,xsph);
  switch(dc->cellmode){
    case CELLMODE_Hneigs:{
      //-Interaction Fluid-Fluid & Fluid-Boundary.
      unsigned bsize=(xsph? dc->bsforcesfluid: dc->bsforcesfluidcorr);
      dim3 sgridf=CsGetGridSize(dc->npf,bsize);
      TmgStart(dc->timers,TMG_CfForcesFluid);
      if(kgc||tvisco==VISCO_LaminarSPS)KerComputeForcesFullFluidNeigs<tkernel,tvisco,kgc,xsph,ncellnv> <<<sgridf,bsize>>> (dc->npb,dc->npf,dc->nctotmax-1,dc->cellf,dc->cellnvb,dc->cellnvf,dc->pospres,dc->velrhop,dc->idp,dc->ace,dc->ar,dc->velxcor,dc->viscdt,dc->matkgc,dc->tau,dc->csph);
      else KerComputeForcesFluidNeigs<tkernel,floating,xsph,ncellnv> <<<sgridf,bsize>>> (dc->npb,dc->npf,dc->nctotmax-1,dc->cellf,dc->cellnvb,dc->cellnvf,dc->pospres,dc->velrhop,dc->idp,dc->ace,dc->ar,dc->velxcor,dc->viscdt);
      TmgStop(dc->timers,TMG_CfForcesFluid);
      //-Interaction Boundary-Fluid.
      bsize=dc->bsforcesbound;
      dim3 sgridb=CsGetGridSize(dc->npb,bsize);
      TmgStart(dc->timers,TMG_CfForcesBound);
      KerComputeForcesBoundNeigs<tkernel,floating,ncellnv> <<<sgridb,bsize>>> (dc->npb,dc->nctotmax-1,dc->cellb,dc->cellnvf,dc->pospres,dc->velrhop,dc->idp,dc->ar,dc->viscdt);
      TmgStop(dc->timers,TMG_CfForcesBound);
    }break;
    case CELLMODE_2H:
    case CELLMODE_H:{
      //-Interaction Fluid-Fluid & Fluid-Boundary.
      unsigned bsize=(xsph? dc->bsforcesfluid: dc->bsforcesfluidcorr);
      dim3 sgridf=CsGetGridSize(dc->npf,bsize);
      TmgStart(dc->timers,TMG_CfForcesFluid);
      if(kgc||tvisco==VISCO_LaminarSPS)KerComputeForcesFullFluid<tkernel,tvisco,kgc,xsph,(ncellnv==9? 1: 2)> <<<sgridf,bsize>>> (dc->npb,dc->npf,dc->ncx,dc->ncy,dc->ncz,dc->cellf,dc->cellbegb,dc->cellbegf,dc->pospres,dc->velrhop,dc->idp,dc->ace,dc->ar,dc->velxcor,dc->viscdt,dc->matkgc,dc->tau,dc->csph);
      else KerComputeForcesFluid<tkernel,floating,xsph,(ncellnv==9? 1: 2)> <<<sgridf,bsize>>> (dc->npb,dc->npf,dc->ncx,dc->ncy,dc->ncz,dc->cellf,dc->cellbegb,dc->cellbegf,dc->pospres,dc->velrhop,dc->idp,dc->ace,dc->ar,dc->velxcor,dc->viscdt);
      TmgStop(dc->timers,TMG_CfForcesFluid);
      //-Interaction Boundary-Fluid.
      bsize=dc->bsforcesbound;
      dim3 sgridb=CsGetGridSize(dc->npb,bsize);
      TmgStart(dc->timers,TMG_CfForcesBound);
      KerComputeForcesBound<tkernel,floating,(ncellnv==9? 1: 2)> <<<sgridb,bsize>>> (dc->npb,dc->ncx,dc->ncy,dc->ncz,dc->cellb,dc->cellbegf,dc->pospres,dc->velrhop,dc->idp,dc->ar,dc->viscdt);
      TmgStop(dc->timers,TMG_CfForcesBound);
    }break;
    default: throw std::string("CsInteraction_Forces> CellMode unrecognised.");
  }
  //-Computes values of tau[] starting from csph[].
  if(tvisco==VISCO_LaminarSPS){
    dim3 sgridf=CsGetGridSize(dc->npf,BLOCKSIZE);
    KerSPSCalcTau <<<sgridf,BLOCKSIZE>>> (dc->npf,dc->npb,dc->smag,dc->blin,dc->rhop,dc->csph,dc->tau);
  }
  //-Cancels forces in Y direction when 2D case.
  if(dc->simulate2d){
    dim3 sgrid=CsGetGridSize(dc->npf,BLOCKSIZE);
    KerResetAcey <<<sgrid,BLOCKSIZE>>> (dc->npf,dc->npb,dc->ace);
  }
  CheckErrorCuda("Failed while executing kernels of interaction forces.");
}
//==============================================================================
/// Calls CsInteraction_Forces().
//==============================================================================
template<bool xsph,TpKernel tkernel,TpVisco tvisco,bool kgc,bool floating>void CsCallInteraction_Forces3(StDeviceContext *dc,StDeviceCte *cte){
  if(dc->ncellnv==9)CsInteraction_Forces<xsph,tkernel,tvisco,kgc,floating,9>(dc,cte);
  else if(dc->ncellnv==25)CsInteraction_Forces<xsph,tkernel,tvisco,kgc,floating,25>(dc,cte);
}
//==============================================================================
template<bool xsph,TpKernel tkernel,TpVisco tvisco,bool floating>void CsCallInteraction_Forces2(StDeviceContext *dc,StDeviceCte *cte){
  if(dc->kgc)CsCallInteraction_Forces3<xsph,tkernel,tvisco,true,floating>(dc,cte);
  else CsCallInteraction_Forces3<xsph,tkernel,tvisco,false,floating>(dc,cte);
}
//==============================================================================
template<bool xsph,TpKernel tkernel,bool floating>void CsCallInteraction_Forces1(StDeviceContext *dc,StDeviceCte *cte){
  if(dc->tvisco==VISCO_Artificial)CsCallInteraction_Forces2<xsph,tkernel,VISCO_Artificial,floating>(dc,cte);
  else if(dc->tvisco==VISCO_LaminarSPS)CsCallInteraction_Forces2<xsph,tkernel,VISCO_LaminarSPS,floating>(dc,cte);
}
//==============================================================================
template<bool xsph,bool floating>void CsCallInteraction_Forces(StDeviceContext *dc,StDeviceCte *cte){
  if(dc->tkernel==KERNEL_Cubic)CsCallInteraction_Forces1<xsph,KERNEL_Cubic,floating>(dc,cte);
  else if(dc->tkernel==KERNEL_Wendland)CsCallInteraction_Forces1<xsph,KERNEL_Wendland,floating>(dc,cte);
}

//------------------------------------------------------------------------------
/// CUDA KERNEL that prepares variables for Shepard interaction.
//------------------------------------------------------------------------------
template<bool floating> __global__ void KerPreShepard(unsigned pini,unsigned npf,float3 *pos,float *rhop,unsigned *idp,unsigned nbound,float ftposx,float4 *posrhop)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;  
  if(p<npf){
    p+=pini;
    float3 rpos=(floating&&idp[p]<nbound? make_float3(ftposx,0,0): pos[p]); 
    posrhop[p]=make_float4(rpos.x,rpos.y,rpos.z,rhop[p]);
  }
}

//==============================================================================
/// Applies Shepard density filter.
//==============================================================================
template<TpKernel tkernel,bool floating,unsigned ncellnv> void CsRunShepard(StDeviceContext *dc,StDeviceCte *cte){
  TmgStart(dc->timers,TMG_CfShepard);
  dim3 sgrid=CsGetGridSize(dc->npf,BLOCKSIZE);
  KerPreShepard<floating> <<<sgrid,BLOCKSIZE>>>(dc->npb,dc->npf,dc->pos,dc->rhop,dc->idp,dc->nbound,dc->posmin.x-(cte->h*2),dc->pospres);
  //-Shepard interaction.
  const unsigned bsize=dc->bsshepard;
  dim3 sgridf=CsGetGridSize(dc->npf,bsize);
  switch(dc->cellmode){
    case CELLMODE_Hneigs:
      KerComputeForcesShepardNeigs<tkernel,ncellnv> <<<sgridf,bsize>>> (dc->npb,dc->npf,dc->nctotmax-1,dc->cellf,dc->cellnvf,cte->massf,dc->pospres,dc->fdrhop);
    break;
    case CELLMODE_2H:
    case CELLMODE_H:
      KerComputeForcesShepard<tkernel,(ncellnv==9? 1: 2)> <<<sgridf,bsize>>> (dc->npb,dc->npf,dc->ncx,dc->ncy,dc->ncz,dc->cellf,dc->cellbegf,cte->massf,dc->pospres,dc->fdrhop);
    break;
    default: throw std::string("CsRunShepard> CellMode unrecognised.");
  }

  CheckErrorCuda("Failed while executing kernel Shepard.");
  //-Applies new density values.
  cudaMemcpy(dc->rhop+dc->npb,dc->fdrhop,sizeof(float)*dc->npf,cudaMemcpyDeviceToDevice);
  TmgStop(dc->timers,TMG_CfShepard);
}

//==============================================================================
/// Calls CsRunShepard().
//==============================================================================
template<TpKernel tkernel,bool floating> void CsCallRunShepard3(StDeviceContext *dc,StDeviceCte *cte){
  if(dc->ncellnv==9)CsRunShepard<tkernel,floating,9>(dc,cte);
  else if(dc->ncellnv==25)CsRunShepard<tkernel,floating,25>(dc,cte);
}
//==============================================================================
template<TpKernel tkernel> void CsCallRunShepard2(StDeviceContext *dc,StDeviceCte *cte){
  if(dc->nfloat)CsCallRunShepard3<tkernel,true>(dc,cte);
  else CsCallRunShepard3<tkernel,false>(dc,cte);
}
//==============================================================================
void CsCallRunShepard(StDeviceContext *dc,StDeviceCte *cte){
  if(dc->tkernel==KERNEL_Cubic)CsCallRunShepard2<KERNEL_Cubic>(dc,cte);
  else if(dc->tkernel==KERNEL_Wendland)CsCallRunShepard2<KERNEL_Wendland>(dc,cte);
}






