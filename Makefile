all: WaveletSampleImage.class test_haar_cpu normalize convert haar

# enable this to generate code for multiple card generations
NVCC=nvcc -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30

# enable one of these to generate code for just one generation of GPU
# (reduces compile time by 30%)
# NVCC=nvcc -arch sm_20
# NVCC=nvcc -arch sm_30

WaveletSampleImage.class: WaveletSampleImage.java
	javac $<

IS_CYGWIN=

ifeq "$(shell uname | head -c 9)" "CYGWIN_NT"

IS_CYGWIN=YES
HAAR_OBJS=haar.obj dwt_cpu.obj dwt_gpu.obj data_io.obj transpose_gpu.obj nixtimer.obj
# NVCC_LIBS=-lws2_32
%.obj: %.cc
	$(NVCC) -c $<
%.obj: %.cu
	$(NVCC) -c $<
CLASSPATH_DIR="$(shell cygpath --windows `pwd`)"

else

HAAR_OBJS=haar.o dwt_cpu.o dwt_gpu.o data_io.o transpose_gpu.o nixtimer.o
LIBS=-lrt
%.o: %.cc
	$(NVCC) -c $<
%.o: %.cu
	$(NVCC) -c $<
CLASSPATH_DIR=$(CURDIR)

endif

haar: $(HAAR_OBJS)
	$(NVCC) -g $^ -o $@

test_haar_cpu: test_haar_cpu.cc dwt_cpu.cc data_io.cc
	gcc -Wall -g $^ -o $@ -lstdc++ $(LIBS)

normalize: normalize.cc data_io.cc
	gcc -Wall -g $^ -o $@ -lstdc++ $(LIBS)

convert: Makefile
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
	rm -f *.class *.obj *.o *.exp *.lib *.pdb *~ convert \
	  haar{,.exe} test_haar_cpu{,.exe} normalize{,.exe}

