default: haar

all: haar WaveletSampleImage.class test_haar_cpu normalize convert \
  cudahaar.mex test_huffman test_haar_thresh_quantUnif_cpu \
  test_haar_thresh_quantLog_cpu test_bit_stream test_compress

java: WaveletSampleImage.class

oct: cudahaar.mex

# enable this to generate code for multiple card generations
NVCC=nvcc \
  -gencode arch=compute_11,code=sm_11 \
  -gencode arch=compute_20,code=sm_20 \
  -gencode arch=compute_30,code=sm_30 \
  -gencode arch=compute_35,code=sm_35

# enable one of these to generate code for just one generation of GPU
# (reduces compile time by 30%)
# NVCC=nvcc -arch sm_20
# NVCC=nvcc -arch sm_30

CC = g++ -std=c++11 -Wall -g

MKOCT=mkoctfile

WaveletSampleImage.class: WaveletSampleImage.java
	javac $<


IS_CYGWIN=

ifeq "$(shell uname | head -c 9)" "CYGWIN_NT"

IS_CYGWIN=YES
CUDA_OBJS=dwt_cpu.obj dwt_gpu.obj data_io.obj transpose_gpu.obj nixtimer.obj
# NVCC_LIBS=-lws2_32
LIBS=-lstdc++
%.obj: %.cc
	$(NVCC) -c $<
%.obj: %.cu
	$(NVCC) -c $<
CLASSPATH_DIR="$(shell cygpath --windows `pwd`)"

else

CUDA_OBJS=dwt_cpu.o dwt_gpu.o data_io.o transpose_gpu.o nixtimer.o
LIBS=-lstdc++ -lrt
NVCC_OCT_OPT=--compiler-options -fPIC
%.o: %.cc
	$(NVCC) $(NVCC_OCT_OPT) -c $<
%.o: %.cu
	$(NVCC) $(NVCC_OCT_OPT) -c $<
CLASSPATH_DIR=$(CURDIR)

endif


haar: $(CUDA_OBJS) haar.cu
	$(NVCC) -g $^ -o $@

test_haar_cpu: test_haar_cpu.cc dwt_cpu.cc data_io.cc
	gcc -Wall -g $^ -o $@ $(LIBS)

test_haar_thresh_quantUnif_cpu: test_haar_thresh_quantUnif_cpu.cc \
  dwt_cpu.cc dwt_cpu.h data_io.cc data_io.h nixtimer.cc nixtimer.h \
  thresh_cpu.cc thresh_cpu.h quant_unif_cpu.cc quant_unif_cpu.h \
  dquant_unif_cpu.cc dquant_unif_cpu.h
	$(CC) test_haar_thresh_quantUnif_cpu.cc dwt_cpu.cc data_io.cc \
	  nixtimer.cc thresh_cpu.cc quant_unif_cpu.cc dquant_unif_cpu.cc \
	  $(LIBS) -o $@

test_haar_thresh_quantLog_cpu: test_haar_thresh_quantLog_cpu.cc \
  dwt_cpu.cc dwt_cpu.h data_io.cc data_io.h nixtimer.cc nixtimer.h \
  thresh_cpu.cc thresh_cpu.h quant_log_cpu.cc quant_log_cpu.h \
  dquant_log_cpu.cc dquant_log_cpu.h
	$(CC) test_haar_thresh_quantLog_cpu.cc dwt_cpu.cc data_io.cc \
	  nixtimer.cc thresh_cpu.cc quant_log_cpu.cc dquant_log_cpu.cc \
	  $(LIBS) -o $@

normalize: normalize.cc data_io.cc
	gcc -Wall -g $^ -o $@ $(LIBS)

test_rle: test_rle.cc rle.h data_io.cc data_io.h huffman.h huffman.cc
	$(CC) test_rle.cc huffman.cc data_io.cc -o $@ $(LIBS)

test_huffman: test_huffman.cc huffman.cc huffman.h
	$(CC) test_huffman.cc huffman.cc -o $@ $(LIBS)

test_bit_stream: test_bit_stream.cc bit_stream.h nixtimer.h nixtimer.cc
	$(CC) test_bit_stream.cc nixtimer.cc -o $@ $(LIBS)

test_quant_count: test_quant_count.cc quant_count.h quant_count.cc
	$(CC) test_quant_count.cc quant_count.cc data_io.cc -o $@ $(LIBS)

test_compress: test_compress.cc dwt_cpu.cc data_io.cc \
	quant_unif_cpu.cc quant_log_cpu.cc nixtimer.cc \
	dquant_unif_cpu.cc dquant_log_cpu.cc thresh_cpu.cc \
	dwt_cpu.h data_io.h \
	quant_unif_cpu.h quant_log_cpu.h \
	dquant_unif_cpu.h dquant_log_cpu.h thresh_cpu.h \
	bit_stream.h nixtimer.h rle.h param_string.h param_string.cc \
	quant_count.h quant_count.cc
	$(CC) test_compress.cc dwt_cpu.cc thresh_cpu.cc \
	  quant_unif_cpu.cc quant_log_cpu.cc quant_count.cc \
	  dquant_unif_cpu.cc dquant_log_cpu.cc \
	  data_io.cc nixtimer.cc param_string.cc \
	  -o $@ $(LIBS)

proto: wavelet_compress_pb.h
wavelet_compress.pb.h: wavelet_compress.proto
	protoc $< --cpp_out=.

list_data: list_data.cc data_io.cc data_io.h
	$(CC) list_data.cc data_io.cc -o $@ $(LIBS)

libwaveletcuda.so: $(CUDA_OBJS) Octave/octave_wrapper.cu
	$(NVCC) $(NVCC_OCT_OPT) -I. -c Octave/octave_wrapper.cu
	$(NVCC) $(NVCC_OCT_OPT) -o $@ --shared $(CUDA_OBJS) octave_wrapper.o

cudahaar.mex: Octave/cudahaar.cc libwaveletcuda.so
	$(MKOCT) -L. -lwaveletcuda Octave/cudahaar.cc

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
	rm -f *.class *.obj *.o *.exp *.lib *.pdb *~ \
	  convert haar test_haar_cpu normalize libwaveletcuda.so cudahaar.oct \
	  test_haar_thresh_quantUnif_cpu test_haar_thresh_quantLog_cpu \
	  test_bit_stream test_compress

