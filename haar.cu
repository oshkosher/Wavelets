#include <cmath>
#include "cuda.h"
#include "cucheck.h"
#include "data_io.h"
#include "nixtimer.h"
#include "dwt_cpu.h"
#include "dwt_gpu.h"
#include "transpose_gpu.h"

#define NUM float

bool similar(NUM a, NUM b) {
  return fabs(a-b) < .00001;
}


void printHelp() {
  fprintf(stderr, 
          "\n"
          "  haar [options] <input datafile> <output datafile> [steps]\n"
          "  Do a Haar discrete wavelet transform.\n"
          "  Options:\n"
          "    -inverse : invert the transform\n"
          "    -text : output data in text format rather than binary\n"
          "    -blocksize : specify the thread block size\n"
          "    -gpu <id>|list : specify the GPU to use, or list all\n"
          "  By default, one transformation step will be done.\n"
          "\n");
  exit(1);
}


int getBestThreadBlockSize(int imageSize) {
  if (imageSize >= 4096) {
    return 1024;
  } else if (imageSize < 512) {
    return 128;
  } else {
    // round imageSize/4 to the nearest power of 2
    return 1 << (int)(log2((double)imageSize) - 2 + .5);
  }
}


void listGpus() {
  int gpuCount;
  cudaDeviceProp prop;
  CUCHECK(cudaGetDeviceCount(&gpuCount));

  for (int gpuId=0; gpuId < gpuCount; gpuId++) {
    CUCHECK(cudaGetDeviceProperties(&prop, gpuId));
    printf("GPU %d: %s, %.1f MHz, %d MB\n", 
           gpuId, prop.name, prop.clockRate / 1000.0, 
           (int)(prop.totalGlobalMem / (1024*1024)));
  }
}


int main(int argc, char **argv) {
  if (argc < 3) printHelp();

  bool inverse = false, textOutput = false;
  int argNo = 1, blockSize = -1, gpuId = 0;

  int gpuCount;
  CUCHECK(cudaGetDeviceCount(&gpuCount));

  while (argNo < argc && argv[argNo][0] == '-') {
    if (!strcmp(argv[argNo], "-inverse")) {
      inverse = true;
      argNo++;
    }

    else if (!strcmp(argv[argNo], "-text")) {
      textOutput = true;
      argNo++;
    }

    else if (!strcmp(argv[argNo], "-blocksize")) {
      if (argNo >= argc) printHelp();
      if (1 != sscanf(argv[++argNo], "%d", &blockSize) ||
          blockSize < 1) {
        printf("Invalid block size \"%s\"\n", argv[argNo]);
        return 1;
      }
      argNo++;
    }

    else if (!strcmp(argv[argNo], "-gpu")) {
      if (argNo >= argc) printHelp();
      argNo++;
      if (!strcmp(argv[argNo], "list")) {
        listGpus();
        return 0;
      } else {
        if (1 != sscanf(argv[argNo], "%d", &gpuId)
            || gpuId < 0
            || gpuId >= gpuCount) {
          printf("Invalid gpu id \"%s\"\n", argv[argNo]);
          return 1;
        }
      }
      argNo++;
    }

    else printHelp();
  }

  // not enough arguments for the input file and output file
  if (argNo+2 > argc) printHelp();

  // read the input file
  const char *inputFilename = argv[argNo++];
  const char *outputFilename = argv[argNo++];
  int stepCount = 1;

  if (argNo < argc) {
    const char *stepsArg = argv[argNo++];
    if (1 != sscanf(stepsArg, "%d", &stepCount)) {
      printf("Invalid step count: \"%s\"\n", stepsArg);
      return 1;
    }
  }
  if (argNo < argc) printHelp();

  NUM *data_cpu, *data_gpu, elapsed;
  int width, height;
  printf("Reading %s...", inputFilename);
  fflush(stdout);
  if (!readDataFile(inputFilename, &data_cpu, &width, &height)) return 1;
  printf("%d x %d\n", width, height);
  fflush(stdout);

  if (width != height) {
    printf("Error: only square data is currently supported.\n");
    return 1;
  }

  int size = width;
  
  CUCHECK(cudaSetDevice(gpuId));
  cudaDeviceProp prop;
  CUCHECK(cudaGetDeviceProperties(&prop, gpuId));
  printf("GPU %d: %s\n", gpuId, prop.name);

  // Make a copy of the data for the GPU to use.
  // Allocate page-locked virtual memory (that won't be moved from its
  // position in physical memory) so the data can be copied to the GPU
  // via DMA This approximately double the throughput.  Just be sure
  // to free the data with cudaFreeHost() rather than delete[].
  CUCHECK(cudaMallocHost((void**)&data_gpu, size*size*sizeof(NUM)));
  memcpy(data_gpu, data_cpu, sizeof(NUM)*size*size);

  // run the CPU version of the algorithm
  printf("CPU: "); fflush(stdout);
  elapsed = haar_2d(data_cpu, size, size, inverse, stepCount);
  printf("%.3f ms\n", elapsed);

  // run the GPU version of the algorithm
  if (blockSize == -1) blockSize = getBestThreadBlockSize(size);

  elapsed = haar_2d_cuda(size, data_gpu, inverse, stepCount,
                                     blockSize, false);

  // Alternative implementation using surfaces.
  // For all inputs I tested, this is slightly slower.
  // elapsed = haar_not_lifting_2d_cuda_surfaces(size, data_gpu, inverse,
  // stepCount, blockSize);

  printf("CUDA: %.6f ms\n", elapsed);

  /*
    // try a variety of thread block sizes
  NUM *data_gpu_copy = new NUM[size*size];
  memcpy(data_gpu_copy, data_gpu, sizeof(NUM)*height*width);
  for (int threadBlockSize = 32; threadBlockSize <= 1024; threadBlockSize*=2) {
    printf("Thread block size: %d\n", threadBlockSize);
    memcpy(data_gpu, data_gpu_copy, sizeof(NUM)*height*width);
    elapsed = haar_not_lifting_2d_cuda(size, data_gpu, inverse, stepCount,
                                       threadBlockSize);
    printf("CUDA: %.6f ms\n", elapsed);
  }
  delete[] data_gpu_copy
  */

  double totalErr = 0;
  for (int i=0; i < size*size; i++) {
    totalErr += fabs(data_cpu[i] - data_gpu[i]);
  }
  
  double averageErr = totalErr / (size*size);

  if (averageErr < 0.000001) {

    // if the CPU version and the GPU version produced similar results,
    // output the requested file.
    writeDataFile(outputFilename, data_gpu, size, size, !textOutput);
    printf("Wrote %s\n", outputFilename);

  } else {

    // if the results look bad, output two files; one with the CPU results
    // and one with the GPU results.
    printf("Average error = %.7f\n", averageErr);
    
    /*
      printf("CPU:\n");
      printMatrix(width, height, data_cpu);
      
      printf("GPU:\n");
      printMatrix(width, height, data_gpu);
    */
    
    writeDataFile("err_cpu.data", data_cpu, size, size, !textOutput);
    writeDataFile("err_gpu.data", data_gpu, size, size, !textOutput);
    printf("Wrote err_cpu.data and err_gpu.data\n");
  }

  delete[] data_cpu;
  CUCHECK(cudaFreeHost(data_gpu));

  return 0;
}
