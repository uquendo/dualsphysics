#DualSPHysics compilation 

#Linux architecture
#ARCH=32
ARCH=64

#To include OpenMP libraries
#WITHOMP=0
WITHOMP=1

EXECNAME=DualSPHysics_linux$(ARCH)

JLIBS=-L./ -ljxml_$(ARCH) -ljformatfiles2_$(ARCH) -ljsphmotion_$(ARCH)

OBJ_CPU=main.o Functions.o JCfgRun.o JDivideCpu.o JException.o JLog.o JObject.o JPartData.o JSpaceCtes.o JSpaceEParms.o JSpaceParts.o JSph.o JSphCpu.o JVarsAscii.o 
OBJ_CUDA=CudaSphApi.o
OBJ_GPU=JPtxasInfo.o JSphGpu.o

OBJECTS=$(OBJ_CPU) $(OBJ_GPU) $(OBJ_CUDA)

GENCODE:=$(GENCODE) -gencode=arch=compute_10,code=\"sm_10,compute_10\"
GENCODE:=$(GENCODE) -gencode=arch=compute_12,code=\"sm_12,compute_12\"
GENCODE:=$(GENCODE) -gencode=arch=compute_20,code=\"sm_20,compute_20\"

NVCC=nvcc
NVCCFLAGS=-c $(GENCODE) -use_fast_math -O3 -I./
NVCCLINKFLAGS=

ifeq ($(ARCH),64)
NVCCFLAGS := $(NVCCFLAGS) -m64
NVCCLINKFLAGS := $(NVCCLINKFLAGS) -m64
endif

ifeq ($(WITHOMP),1)
NVCCFLAGS := $(NVCCFLAGS) -Xcompiler -fopenmp
NVCCLINKFLAGS := $(NVCCLINKFLAGS) -Xcompiler -fopenmp -lgomp
endif

all:dualsphysics 
	mv dualsphysics $(EXECNAME)
	mv dualsphysics_ptxasinfo $(EXECNAME)_ptxasinfo
	rm -rf *.o

dualsphysics: $(OBJECTS)
	$(NVCC) ${NVCCLINKFLAGS} $(OBJECTS) -o $@ $(JLIBS)

.cpp.o: 
	$(NVCC) $(NVCCFLAGS) $< 

CudaSphApi.o: CudaSphApi.cu
	$(NVCC) $(NVCCFLAGS) --ptxas-options -v CudaSphApi.cu 2>dualsphysics_ptxasinfo

clean:
	rm -rf *.o DualSPHysics*


