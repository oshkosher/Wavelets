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
#     find it as pre-built package (look for protobuf and libprotobuf-devel)
#     or to build it from the source and install it:
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

default: convert test_haar_cpu haar test_compress_cpu test_compress_gpu

EXECS = test_haar_cpu haar test_compress \
  test_haar_thresh_quantUnif_cpu test_haar_thresh_quantLog_cpu \
  normalize test_rle test_huffman test_bit_stream test_quant_count \
  test_compress_gpu list_data image_error test_transform test_lloyd

all: convert $(EXECS) libwaveletcuda.so cudahaar.mex

java: WaveletSampleImage.class ImageDiff.class

oct: cudahaar.mex

# Set this to YES or NO, to select between a Debug or Release build
IS_DEBUG=YES
# IS_DEBUG=NO


ifeq ($(IS_DEBUG),YES)
BUILD=Debug
CC_OPT_FLAG=-g
CL_OPT_FLAG=/MDd
else
BUILD=Release
CC_OPT_FLAG=-O2
CL_OPT_FLAG=/MD
endif

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
PROTOC_VC = $(PROTOBUF_DIR_VC)/$(BUILD)/protoc.exe
PROTOC = protoc
NVCC_ARCH_SIZE = -m32
NVCC_OPT = --compiler-options $(CL_OPT_FLAG) -D_SCL_SECURE_NO_WARNINGS  -I$(PROTOBUF_DIR_VC)/include $(NVCC_COMPILER_BINDIR)
CLASSPATH_DIR="$(shell cygpath --windows `pwd`)"

else

OBJ_EXT=o
LIBS=-lstdc++ -lrt
PROTOBUF_DIR = /usr/local
PROTOBUF_LIB = -L$(PROTOBUF_DIR)/lib -lprotobuf
PROTOBUF_LIB_NVCC = $(PROTOBUF_LIB)
PROTOC = protoc
NVCC_OPT=--compiler-options -fPIC
CLASSPATH_DIR=$(CURDIR)

endif

NVCC = nvcc $(NVCC_OPT) $(NVCC_ARCH) $(NVCC_ARCH_SIZE) $(NVCC_COMPILER_BINDIR) $(CC_OPT_FLAG) 

%.$(OBJ_EXT): %.cc
	$(NVCC) -c $<

%.$(OBJ_EXT): %.cu
	$(NVCC) -c $<

CUDA_OBJS=dwt_cpu.$(OBJ_EXT) dwt_gpu.$(OBJ_EXT) data_io.$(OBJ_EXT) \
  transpose_gpu.$(OBJ_EXT) nixtimer.$(OBJ_EXT)

HAAR_OBJS=haar.$(OBJ_EXT) dwt_cpu.$(OBJ_EXT) dwt_gpu.$(OBJ_EXT) \
  data_io.$(OBJ_EXT) transpose_gpu.$(OBJ_EXT) nixtimer.$(OBJ_EXT)

haar: $(HAAR_OBJS)
	$(NVCC) $(HAAR_OBJS) -o $@

WaveletSampleImage.class: WaveletSampleImage.java
	javac $<

ImageDiff.class: ImageDiff.java
	javac $<

test_haar_cpu: test_haar_cpu.cc dwt_cpu.cc data_io.cc
	$(CC) $^ -o $@ $(LIBS)

test_lloyd: test_lloyd.cc Octave/LloydsAlgorithm/src/c++/lloyds.cpp \
	  Octave/LloydsAlgorithm/src/c++/lloyds.h
	$(CC) test_lloyd.cc Octave/LloydsAlgorithm/src/c++/lloyds.cpp \
	  -IOctave/LloydsAlgorithm/src/c++ -o $@

test_haar_thresh_quantUnif_cpu: test_haar_thresh_quantUnif_cpu.cc \
  dwt_cpu.cc dwt_cpu.h data_io.cc data_io.h nixtimer.cc nixtimer.h \
  thresh_cpu.cc thresh_cpu.h quant_unif_cpu.cc quant_unif_cpu.h \
  dquant_unif_cpu.cc dquant_unif_cpu.h quant.h quant.cc
	$(CC) test_haar_thresh_quantUnif_cpu.cc dwt_cpu.cc data_io.cc \
	  nixtimer.cc thresh_cpu.cc quant_unif_cpu.cc dquant_unif_cpu.cc \
	  quant.cc $(LIBS) -o $@

test_haar_thresh_quantLog_cpu: test_haar_thresh_quantLog_cpu.cc \
  dwt_cpu.cc dwt_cpu.h data_io.cc data_io.h nixtimer.cc nixtimer.h \
  thresh_cpu.cc thresh_cpu.h quant_log_cpu.cc quant_log_cpu.h \
  dquant_log_cpu.cc dquant_log_cpu.h quant.h quant.cc
	$(CC) test_haar_thresh_quantLog_cpu.cc dwt_cpu.cc data_io.cc \
	  nixtimer.cc thresh_cpu.cc quant_log_cpu.cc dquant_log_cpu.cc \
	  quant.cc $(LIBS) -o $@

normalize: normalize.cc data_io.cc
	$(CC) $^ -o $@ $(LIBS)

image_error: image_error.cc data_io.cc
	$(CC) $^ -o $@ $(LIBS)

test_rle: test_rle.cc rle.h data_io.cc data_io.h huffman.h huffman.cc
	$(CC) test_rle.cc huffman.cc data_io.cc -o $@ $(LIBS)

test_huffman: test_huffman.cc huffman.cc huffman.h
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

test_compress: test_compress_cpu

test_compress_cpu: test_compress_cpu.cc test_compress_common.cc \
	dwt_cpu.cc data_io.cc \
	quant_unif_cpu.cc quant_log_cpu.cc nixtimer.cc \
	dquant_unif_cpu.cc dquant_log_cpu.cc thresh_cpu.cc \
	dwt_cpu.h data_io.h \
	quant_unif_cpu.h quant_log_cpu.h \
	dquant_unif_cpu.h dquant_log_cpu.h thresh_cpu.h \
	bit_stream.h nixtimer.h rle.h param_string.h param_string.cc \
	quant_count.h quant_count.cc quant.h quant.cc \
	wavelet_compress.pb.h wavelet_compress.pb.cc \
	Octave/LloydsAlgorithm/src/c++/lloyds.cpp \
	Octave/LloydsAlgorithm/src/c++/lloyds.h
	$(CC) test_compress_cpu.cc test_compress_common.cc \
	  dwt_cpu.cc thresh_cpu.cc \
	  quant_unif_cpu.cc quant_log_cpu.cc quant_count.cc quant.cc \
	  dquant_unif_cpu.cc dquant_log_cpu.cc param_string.cc \
	  data_io.cc nixtimer.cc wavelet_compress.pb.cc \
	  Octave/LloydsAlgorithm/src/c++/lloyds.cpp \
	  -IOctave/LloydsAlgorithm/src/c++ \
	  -o $@ $(LIBS) $(PROTOBUF_LIB)

test_compress_gpu.$(OBJ_EXT): test_compress_gpu.cu wavelet_compress.pb.h quant.h

test_compress_common.$(OBJ_EXT): test_compress_common.cc test_compress_common.h  rle.h param_string.h

TEST_COMPRESS_GPU_OBJS=test_compress_gpu.$(OBJ_EXT) \
  test_compress_common.$(OBJ_EXT) \
  dwt_cpu.$(OBJ_EXT) dwt_gpu.$(OBJ_EXT) \
  data_io.$(OBJ_EXT) transpose_gpu.$(OBJ_EXT) nixtimer.$(OBJ_EXT) \
  wavelet_compress.pb.$(OBJ_EXT) quant_count.$(OBJ_EXT) param_string.$(OBJ_EXT)

test_compress_gpu: $(TEST_COMPRESS_GPU_OBJS)
	$(NVCC) $(TEST_COMPRESS_GPU_OBJS) -o $@ $(PROTOBUF_LIB_NVCC)

proto: wavelet_compress.pb.h
wavelet_compress.pb.h wavelet_compress.pb.cc: wavelet_compress.proto
	$(PROTOC) $< --cpp_out=.

wavelet_compress.pb.$(OBJ_EXT): wavelet_compress.pb.cc wavelet_compress.pb.h

list_data: list_data.cc data_io.cc data_io.h
	$(CC) list_data.cc data_io.cc -o $@ $(LIBS)

libwaveletcuda.so: $(CUDA_OBJS) Octave/LloydsAlgorithm/octave_wrapper.cu
	$(NVCC) -I. -c Octave/LloydsAlgorithm/octave_wrapper.cu
	$(NVCC) -o $@ --shared $(CUDA_OBJS) octave_wrapper.$(OBJ_EXT)

cudahaar.mex: Octave/LloydsAlgorithm/cudahaar.cc libwaveletcuda.so
	$(MKOCT) -L. -lwaveletcuda Octave/LloydsAlgorithm/cudahaar.cc

convert: Makefile WaveletSampleImage.class
	@echo Write $@ wrapper for \"java WaveletSampleImage\"
	@echo '#!/bin/sh' > convert
	@echo 'unset DISPLAY' >> convert
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
	rm -f *.class *.obj *.o *.exp *.lib *.pdb *.so *~ $(EXECS) \
	  convert libwaveletcuda.so cudahaar.oct \
	  wavelet_compress.pb.{h,cc}
