default: haar

all: haar WaveletSampleImage.class test_haar_cpu normalize convert \
  cudahaar.mex test_huffman test_bit_stream

java: WaveletSampleImage.class

oct: cudahaar.mex

# enable this to generate code for multiple card generations
NVCC=nvcc --compiler-options -fPIC -gencode arch=compute_20,code=sm_20 \
  -gencode arch=compute_30,code=sm_30 \
  -gencode arch=compute_35,code=sm_35

# enable one of these to generate code for just one generation of GPU
# (reduces compile time by 30%)
# NVCC=nvcc -arch sm_20
# NVCC=nvcc -arch sm_30

CC = gcc -Wall -g

MKOCT=mkoctfile

WaveletSampleImage.class: WaveletSampleImage.java
	javac $<


IS_CYGWIN=

ifeq "$(shell uname | head -c 9)" "CYGWIN_NT"

IS_CYGWIN=YES
CUDA_OBJS=dwt_cpu.obj dwt_gpu.obj data_io.obj transpose_gpu.obj nixtimer.obj
# NVCC_LIBS=-lws2_32
%.obj: %.cc
	$(NVCC) -c $<
%.obj: %.cu
	$(NVCC) -c $<
CLASSPATH_DIR="$(shell cygpath --windows `pwd`)"

else

CUDA_OBJS=dwt_cpu.o dwt_gpu.o data_io.o transpose_gpu.o nixtimer.o
LIBS=-lrt
%.o: %.cc
	$(NVCC) -c $<
%.o: %.cu
	$(NVCC) -c $<
CLASSPATH_DIR=$(CURDIR)

endif


haar: $(CUDA_OBJS) haar.cu
	$(NVCC) -g $^ -o $@

test_haar_cpu: test_haar_cpu.cc dwt_cpu.cc data_io.cc
	gcc -Wall -g $^ -o $@ -lstdc++ $(LIBS)

normalize: normalize.cc data_io.cc
	gcc -Wall -g $^ -o $@ -lstdc++ $(LIBS)

test_rle: test_rle.cc rle.h data_io.cc data_io.h huffman.h huffman.cc
	$(CC) test_rle.cc huffman.cc data_io.cc -o $@ -lstdc++ $(LIBS)

test_huffman: test_huffman.cc huffman.cc huffman.h
	$(CC) test_huffman.cc huffman.cc -o $@ -lstdc++ $(LIBS)

test_bit_stream: test_bit_stream.cc bit_stream.h nixtimer.h nixtimer.cc
	$(CC) test_bit_stream.cc nixtimer.cc -o $@ -lstdc++ $(LIBS)

libwaveletcuda.so: $(CUDA_OBJS) Octave/octave_wrapper.cu
	$(NVCC) -I. -c Octave/octave_wrapper.cu
	$(NVCC) -o $@ --shared $(CUDA_OBJS) octave_wrapper.o

cudahaar.mex: Octave/cudahaar.cc libwaveletcuda.so
	$(MKOCT) -L. -lwaveletcuda Octave/cudahaar.cc

convert: Makefile WaveletSampleImage.class
	@echo Write $@ wrapper for \"java WaveletSampleImage\"
	@echo '#!/bin/sh' > convert
	@echo java -cp \"$(CLASSPATH_DIR)\" WaveletSampleImage \"\$$@\" >> convert
	chmod 755 convert

send: .senddevel

.senddevel: *.cu *.cc *.h *.java Makefile
	scp $? devel:tmp/Wavelets
	touch .senddevel

sendscu: .sendscu

.sendscu: *.cu *.cc *.h *.java Makefile
	scp $? scu:Wavelets
	touch .sendscu

clean:
	rm -f *.class *.obj *.o *.exp *.lib *.pdb *~ \
	  convert haar test_haar_cpu normalize libwaveletcuda.so cudahaar.oct

