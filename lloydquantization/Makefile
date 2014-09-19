CC=gcc
CC_FLAGS= -Wall
OCT=mkoctfile
NVCC=nvcc
NV_FLAGS=-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\" --compiler-options '-fPIC'
#-gencode arch=compute_10,code=sm_10 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\" --compiler-options '-fPIC'

CPATH=src/oct/
CUDAPATH=src/cuda/
OCTNAME=lloydcuda
SOURCE=$(CPATH)$(OCTNAME).cpp
CSOURCES=$(wildcard $(CUDAPATH)*.cpp)
COBJECTS=$(CSOURCES:.cpp=.o)
CUDANAME=lloyd
CUDASOURCES=$(wildcard $(CUDAPATH)*.cu)
CUDAOBJECTS=$(CUDASOURCES:.cu=.o)
CUDALIB=lib$(CUDANAME).so


all:cuda
	$(OCT) -L. -l$(CUDANAME) $(SOURCE)

cuda:$(CUDAOBJECTS) $(COBJECTS)
	$(NVCC) -o $(CUDALIB) --shared $(CUDAOBJECTS) $(COBJECTS)

%.o:%.cpp
	$(NVCC) $(NV_FLAGS) -o $@ -c $<

%.o:%.cu
	$(NVCC) $(NV_FLAGS) -o $@ -c $<

clean:
	rm -f $(CUDAOBJECTS) $(COBJECTS) $(CUDALIB) $(OCTNAME).oct $(OCTNAME).o

