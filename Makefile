# Requirements:
#   nvcc (NVIDIA CUDA compiler)
#     Install the CUDA toolkit from:
#     https://developer.nvidia.com/cuda-downloads
#     Add nvcc to your PATH.
#     
#   javac (Java compiler)
#     Download the Java development kit (aka Java Platform, aka JDK):
#     http://www.oracle.com/technetwork/java/javase/downloads/index.html
#     Add javac to your PATH.
# 
#   protoc (Google Protocol Buffers)
#     On Unix or Unix-like (Cygwin) systems, it's pretty easy to either
#     find it as pre-built package (look for protobuf, libprotobuf-devel,
#     and libprotobuf-java) or to build it from the source and install it:
#       https://developers.google.com/protocol-buffers/docs/downloads
#     On Windows, you'll either need to build it from the source, or download
#     pre-built binaries from our FTP site where we've got large images
#     archived (ftp.dynapsecorp.com). They're in the "Protobuf binaries"
#     subdirectory. Download the ZIP file for your version of Visual Studio
#     and extract the contents into this directory. It will create a
#     directory "protobuf-2.6.0" which contains the binaries for Debug
#     and Release builds.
#
#     If you want to build from the source:
#       1. Download the source and expand the archive in this directory.
#       2. Look in the "vsprojects" subdirectory. Open "protobuf.sln"
#          in Visual Studio. Build the project for Debug and Release.
#          You'll want to disable incremental linking and set the
#          Debug Information Format to Program Database (/Zi).
#       3. In a command prompt, navigate to the "vsprojects" subdirectory
#          and run "extract_includes.bat".
#       4. In Makefile.nmake and Makefile, check that PROTOBUF_DIR,
#          PROTOBUF_LIB, and PROTOC point to the locations matching
#          your build.
#
#   mkoctfile (Octave C interface)
#     Install liboctave-dev package.
# 
# In short, on a fresh Ubuntu 14 install, start by installing these:
#    sudo apt-get install g++ openjdk-7-jdk make git protobuf-compiler libprotobuf-dev libprotobuf-java unzip liboctave-dev
#
#    Then install CUDA
#      http://developer.download.nvidia.com/compute/cuda/6_5/rel/installers/cuda_6.5.14_linux_64.run

# build all the tools that don't require CUDA
default: test_compress_cpu java convert histogram cubelet_convert

EXECS = test_compress_cpu test_haar_cpu haar \
  test_haar_thresh_quantUnif_cpu test_haar_thresh_quantLog_cpu \
  normalize test_rle test_huffman test_bit_stream test_quant_count \
  list_data image_error test_transform test_lloyd \
  histogram test_cubelet_file cubelet_convert test_wavelet

# build the tools that require CUDA
cuda: haar test_compress_gpu

java: makeall.class

makeall.class Cubelet.class CubeletFile.class CubeletViewer.class \
  ImageDiff.class WaveletCompress.class WaveletSampleImage.class \
  : makeall.java Cubelet.java CubeletFile.java CubeletViewer.java \
    ImageDiff.java WaveletCompress.java WaveletSampleImage.java
	$(JAVAC) makeall.java

all: convert $(EXECS) java

# this is the only part that require JavaCV, so if you don't have
# JavaCV installed, this is the only part you can't build
movieconvert: MovieToCubelets.class

octave: libwaveletcuda.so cudahaar.mex

# Set this to YES or NO, to select between a Debug or Release build
# IS_DEBUG=YES
IS_DEBUG=NO


ifeq ($(IS_DEBUG),YES)
BUILD=Debug
CC_OPT_FLAG=-g
CL_OPT_FLAG=/MDd
else
BUILD=Release
CC_OPT_FLAG=-O2 -DNDEBUG
CL_OPT_FLAG=/MD
endif

# Enable OpenMP?
#  -fopenmp

# Either -m32 or -m64. Currently we use 64-bit on Linux and 32-bit on
# Windows, because not everyone in the group has 64-bit support in
# Visual Studio.
NVCC_ARCH_SIZE = -m64

# By default, build CUDA code for many architectures
NVCC_ARCH = \
  -gencode arch=compute_11,code=sm_11 \
  -gencode arch=compute_20,code=sm_20 \
  -gencode arch=compute_30,code=sm_30 \
  -gencode arch=compute_35,code=sm_35

# enable one of these to generate code for just one generation of GPU
# (reduces compile time by 30%)
# NVCC_ARCH=-arch sm_20
NVCC_ARCH=-arch sm_30

# use this to direct NVCC to use a different host compiler, if necessary
# NVCC_COMPILER_BINDIR=--compiler-bindir='C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin'

CC = g++ -std=c++11 -Wall $(CC_OPT_FLAG)

MKOCT=mkoctfile


IS_CYGWIN=

ifeq "$(shell uname | head -c 9)" "CYGWIN_NT"

IS_CYGWIN=YES
OBJ_EXT=obj
# NVCC_LIBS=-lws2_32
LIBS=-lstdc++
PROTOBUF_DIR_VC = protobuf-2.6.0/vsprojects
PROTOBUF_DIR = /usr/local
PROTOBUF_LIB_NVCC = $(PROTOBUF_DIR_VC)/$(BUILD)/libprotobuf.lib
PROTOBUF_LIB = -L$(PROTOBUF_DIR)/lib -lprotobuf
PROTOBUF_JAR = protobuf-2.6.0/protobuf.jar
JAVACV_JAR = c:/Apps/javacv-bin/javacv.jar
PROTOC_VC = $(PROTOBUF_DIR_VC)/$(BUILD)/protoc.exe
PROTOC = protoc
NVCC_ARCH_SIZE = -m32
NVCC_OPT = --compiler-options $(CL_OPT_FLAG) -D_SCL_SECURE_NO_WARNINGS  -I$(PROTOBUF_DIR_VC)/include $(NVCC_COMPILER_BINDIR)
CLASSPATH_DIR="$(shell cygpath --windows `pwd`);protobuf-2.6.0/protobuf.jar;c:/Apps/javacv-bin/javacv.jar"
JAVAC = javac -cp '.;$(PROTOBUF_JAR);$(JAVACV_JAR)'

else

OBJ_EXT=o
LIBS=-lstdc++ -lrt
PROTOBUF_DIR = /usr/local
PROTOBUF_LIB = -L$(PROTOBUF_DIR)/lib -lprotobuf
PROTOBUF_LIB_NVCC = $(PROTOBUF_LIB)
PROTOBUF_JAR = /usr/share/java/protobuf.jar
JAVACV_JAR = /usr/local/javacv/bin/javacv.jar
PROTOC = protoc
NVCC_OPT = -std=c++11
NVCC_SHLIB_OPT=--compiler-options -fPIC
CLASSPATH_DIR=$(CURDIR):$(PROTOBUF_JAR)
JAVAC = javac -cp .:$(PROTOBUF_JAR):$(JAVACV_JAR)

endif

NVCC = nvcc $(NVCC_OPT) $(NVCC_ARCH) $(NVCC_ARCH_SIZE) $(NVCC_COMPILER_BINDIR) $(CC_OPT_FLAG) 

%.$(OBJ_EXT): %.cc
	$(NVCC) -c $<

%.$(OBJ_EXT): %.cu
	$(NVCC) -c $<

CUDA_OBJS=dwt_cpu.$(OBJ_EXT) dwt_gpu.$(OBJ_EXT) data_io.$(OBJ_EXT) \
  transpose_gpu.$(OBJ_EXT) nixtimer.$(OBJ_EXT) octave_wrapper.$(OBJ_EXT) \
  wavelet.$(OBJ_EXT) wavelet_compress.pb.$(OBJ_EXT)

HAAR_OBJS=haar.$(OBJ_EXT) dwt_cpu.$(OBJ_EXT) dwt_gpu.$(OBJ_EXT) \
  data_io.$(OBJ_EXT) transpose_gpu.$(OBJ_EXT) nixtimer.$(OBJ_EXT) \
  wavelet.$(OBJ_EXT) wavelet_compress.pb.$(OBJ_EXT)

LLOYD_INC=-IOctave/LloydsAlgorithm/src/c++

haar.$(OBJ_EXT): wavelet_compress.pb.h haar.cu dwt_cpu.h data_io.h dwt_gpu.h

haar: $(HAAR_OBJS)
	$(NVCC) $(HAAR_OBJS) $(PROTOBUF_LIB_NVCC) -o $@

cubelet: test_cubelet_file MovieToCubelets.class cubelet_convert

MovieToCubelets.class: MovieToCubelets.java
	$(JAVAC) $<

test_haar_cpu.o: test_haar_cpu.cc dwt_cpu.h data_io.h cubelet_file.h \
	  wavelet_compress.pb.h
	$(CC) -c $<

test_haar_cpu: test_haar_cpu.o dwt_cpu.o data_io.o wavelet.o \
	  wavelet_compress.pb.o cubelet_file.o
	$(CC) $^ -o $@ $(LIBS) $(PROTOBUF_LIB)

test_lloyd: test_lloyd.cc Octave/LloydsAlgorithm/src/c++/lloyds.cpp \
	  Octave/LloydsAlgorithm/src/c++/lloyds.h
	$(CC) test_lloyd.cc Octave/LloydsAlgorithm/src/c++/lloyds.cpp \
	  $(LLOYD_INC) -o $@

test_haar_thresh_quantUnif_cpu: test_haar_thresh_quantUnif_cpu.cc \
  dwt_cpu.cc dwt_cpu.h data_io.cc data_io.h nixtimer.cc nixtimer.h \
  thresh_cpu.cc thresh_cpu.h quant_unif_cpu.cc quant_unif_cpu.h \
  dquant_unif_cpu.cc dquant_unif_cpu.h quant.h quant.cc
	$(CC) test_haar_thresh_quantUnif_cpu.cc dwt_cpu.cc data_io.cc \
	  nixtimer.cc thresh_cpu.cc quant_unif_cpu.cc dquant_unif_cpu.cc \
	  quant.cc wavelet.cc wavelet_compress.pb.o \
	  $(PROTOBUF_LIB) $(LIBS) -o $@

test_haar_thresh_quantLog_cpu: test_haar_thresh_quantLog_cpu.cc \
  dwt_cpu.cc dwt_cpu.h data_io.cc data_io.h nixtimer.cc nixtimer.h \
  thresh_cpu.cc thresh_cpu.h quant_log_cpu.cc quant_log_cpu.h \
  dquant_log_cpu.cc dquant_log_cpu.h quant.h quant.cc
	$(CC) test_haar_thresh_quantLog_cpu.cc dwt_cpu.cc data_io.cc \
	  nixtimer.cc thresh_cpu.cc quant_log_cpu.cc dquant_log_cpu.cc \
	  quant.cc wavelet.cc wavelet_compress.pb.o \
	  $(PROTOBUF_LIB) $(LIBS) -o $@

normalize: normalize.cc data_io.cc
	$(CC) $^ -o $@ $(LIBS)

image_error: image_error.cc data_io.cc
	$(CC) $^ -o $@ $(LIBS)

test_rle: test_rle.cc rle.h data_io.cc data_io.h huffman.h huffman.cc
	$(CC) test_rle.cc huffman.cc data_io.cc -o $@ $(LIBS)

test_huffman: test_huffman.cc huffman.cc huffman.h bit_stack.h
	$(CC) test_huffman.cc huffman.cc -o $@ $(LIBS)

test_bit_stream: test_bit_stream.cc bit_stream.h nixtimer.h nixtimer.cc
	$(CC) test_bit_stream.cc nixtimer.cc -o $@ $(LIBS)

test_quant_count: test_quant_count.cc quant_count.h quant_count.cc quant.h \
	  quant.cc thresh_cpu.cc quant_unif_cpu.cc quant_log_cpu.cc \
	  dquant_unif_cpu.cc dquant_log_cpu.cc
	$(CC) test_quant_count.cc quant_count.cc quant.cc \
	  data_io.cc thresh_cpu.cc quant_unif_cpu.cc quant_log_cpu.cc \
	  dquant_unif_cpu.cc dquant_log_cpu.cc \
	  -o $@ $(LIBS)

test_transform: test_transform.cu
	$(NVCC) $< -o $@

test_compress: test_compress_cpu test_compress_gpu

test_compress_cpu.o: test_compress_cpu.cc test_compress_common.h \
	  dwt_cpu.h nixtimer.h thresh_cpu.h quant.h bit_stream.h huffman.h \
	  Octave/LloydsAlgorithm/src/c++/lloyds.h wavelet_compress.pb.h \
	  wavelet.h cubelet_file.h
	$(CC) $(LLOYD_INC) -c $<

wavelet.o: wavelet.cc wavelet.h wavelet_compress.pb.h
	$(CC) -c $<

optimize.o: optimize.cc optimize.h test_compress_cpu.h wavelet.h
	$(CC) -c $<

test_compress_common.o: test_compress_common.cc test_compress_common.h \
	  rle.h bit_stream.h nixtimer.h huffman.h dwt_cpu.h \
	  wavelet_compress.pb.h
	$(CC) -c $<

dwt_cpu.o: dwt_cpu.cc dwt_cpu.h nixtimer.h wavelet.h
	$(CC) -c $<

data_io.o: data_io.cc data_io.h
	$(CC) -c $<

nixtimer.o: nixtimer.cc nixtimer.h
	$(CC) -c $<

thresh_cpu.o: thresh_cpu.cc thresh_cpu.h
	$(CC) -c $<

param_string.o: param_string.cc param_string.h
	$(CC) -c $<

quant.o: quant.cc quant.h nixtimer.h
	$(CC) -c $<

huffman.o: huffman.cc huffman.h bit_stream.h bit_stack.h
	$(CC) -c $<

wavelet_compress.pb.o: wavelet_compress.pb.cc wavelet_compress.pb.h
	$(CC) -c $<

lloyds.o: Octave/LloydsAlgorithm/src/c++/lloyds.cpp \
	  Octave/LloydsAlgorithm/src/c++/lloyds.h
	$(CC) -c -IOctave/LloydsAlgorithm/src/c++ $< -o $@

cubelet_file.o: cubelet_file.cc cubelet_file.h wavelet_compress.pb.h
	$(CC) -c $<

test_cubelet_file: test_cubelet_file.cc cubelet_file.o wavelet_compress.pb.o \
	  wavelet.o
	$(CC) $^ $(PROTOBUF_LIB) -o $@

cubelet_convert: cubelet_convert.cc cubelet_file.o wavelet_compress.pb.o \
	  wavelet.o
	$(CC) $^ $(PROTOBUF_LIB) -o $@

test_wavelet: test_wavelet.cc wavelet.o wavelet_compress.pb.o nixtimer.o
	$(CC) $^ $(PROTOBUF_LIB) -o $@

TEST_COMPRESS_CPU_OBJS=test_compress_cpu.o wavelet.o \
	test_compress_common.o dwt_cpu.o data_io.o \
	nixtimer.o thresh_cpu.o \
	quant.o wavelet_compress.pb.o \
	lloyds.o huffman.o cubelet_file.o optimize.o

test_compress_cpu: $(TEST_COMPRESS_CPU_OBJS)
	$(CC) $(TEST_COMPRESS_CPU_OBJS) -o $@ $(LIBS) $(PROTOBUF_LIB)

test_compress_gpu.$(OBJ_EXT): test_compress_gpu.cu wavelet_compress.pb.h quant.h

test_compress_common.$(OBJ_EXT): test_compress_common.cc test_compress_common.h  rle.h

TEST_COMPRESS_GPU_OBJS=test_compress_gpu.$(OBJ_EXT) \
  test_compress_common.$(OBJ_EXT) \
  dwt_cpu.$(OBJ_EXT) dwt_gpu.$(OBJ_EXT) huffman.$(OBJ_EXT) \
  data_io.$(OBJ_EXT) transpose_gpu.$(OBJ_EXT) nixtimer.$(OBJ_EXT) \
  wavelet_compress.pb.$(OBJ_EXT) quant_count.$(OBJ_EXT)

test_compress_gpu: $(TEST_COMPRESS_GPU_OBJS)
	$(NVCC) $(TEST_COMPRESS_GPU_OBJS) -o $@ $(PROTOBUF_LIB_NVCC)

proto: wavelet_compress.pb.h
wavelet_compress.pb.h wavelet_compress.pb.cc WaveletCompress.java: wavelet_compress.proto
	$(PROTOC) $< --cpp_out=. --java_out=.

wavelet_compress.pb.$(OBJ_EXT): wavelet_compress.pb.cc wavelet_compress.pb.h

list_data: list_data.cc data_io.cc data_io.h
	$(CC) list_data.cc data_io.cc -o $@ $(LIBS)

histogram: histogram.cc data_io.cc data_io.h
	$(CC) histogram.cc data_io.cc -o $@ $(LIBS)

test_array_data: test_array_data.cc array_data.h
	$(CC) test_array_data.cc -o $@

libwaveletcuda.so: dwt_cpu.h dwt_cpu.cc dwt_gpu.h dwt_gpu.cu data_io.h data_io.cc \
  transpose_gpu.h transpose_gpu.cu nixtimer.h nixtimer.cc wavelet.cc wavelet.h \
  Octave/LloydsAlgorithm/octave_wrapper.cu \
  wavelet_compress.pb.cc wavelet_compress.pb.h
	rm -f $(CUDA_OBJS)
	$(NVCC) $(NVCC_SHLIB_OPT) -I. -c dwt_cpu.cc dwt_gpu.cu data_io.cc \
	  transpose_gpu.cu nixtimer.cc wavelet.cc wavelet_compress.pb.cc \
	  Octave/LloydsAlgorithm/octave_wrapper.cu 
	$(NVCC) -o $@ --shared $(CUDA_OBJS) $(PROTOBUF_LIB_NVCC)
	rm -f $(CUDA_OBJS)

cudahaar.mex: Octave/LloydsAlgorithm/cudahaar.cc libwaveletcuda.so
	$(MKOCT) -L. -lwaveletcuda Octave/LloydsAlgorithm/cudahaar.cc

convert: Makefile
	@echo Write $@ wrapper for \"java WaveletSampleImage\"
	@echo '#!/bin/sh' > convert
	@echo 'unset DISPLAY' >> convert
	@echo java -cp \"$(CLASSPATH_DIR)\" WaveletSampleImage \"\$$@\" >> convert
	chmod 755 convert

send: .senddevel

.senddevel: *.cu *.cc *.h *.java *.proto Makefile
	scp $? devel:tmp/Wavelets
	touch .senddevel

sendscu: .sendscu

.sendscu: *.cu *.cc *.h *.java *.proto Makefile
	scp $? scu:Wavelets
	touch .sendscu

clean:
	rm -f *.class *.obj *.o *.exp *.lib *.pdb *.so *.gch *.ilk *~ \
	  *.stackdump $(EXECS) convert libwaveletcuda.so cudahaar.oct \
	  wavelet_compress.pb.{h,cc} WaveletCompress.java
