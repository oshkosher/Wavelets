#include <cmath>
#include "cuda.h"
#include "cucheck.h"
#include "dwt_gpu.h"

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

int haar(float *output, float *input, int width, int steps, bool inverse, int blockSize) {

  float *plmemory;
  float elapsed;
  
  cudaDeviceProp prop;
  CUCHECK(cudaGetDeviceProperties(&prop, 0));
  printf("GPU %d: %s\n", 0, prop.name);

  // Make a copy of the data for the GPU to use.
  // Allocate page-locked virtual memory (that won't be moved from its
  // position in physical memory) so the data can be copied to the GPU
  // via DMA This approximately double the throughput.  Just be sure
  // to free the data with cudaFreeHost() rather than delete[].
  CUCHECK(cudaMallocHost((void**)&plmemory, width*width*sizeof(float)));
  memcpy(plmemory, input, sizeof(float)*width*width);

  // run the GPU version of the algorithm
  if (blockSize == -1) blockSize = getBestThreadBlockSize(width);

  elapsed = haar_not_lifting_2d_cuda(width, plmemory, inverse, steps,
                                     blockSize, true);

  memcpy(output, plmemory, sizeof(float)*width*width);

  printf("CUDA: %.6f ms\n", elapsed);

  CUCHECK(cudaFreeHost(plmemory));

  return 0;
}

int haar(double *output, double *input, int width, int steps, bool inverse, int blockSize) {

  double *plmemory;
  float elapsed;
  
  cudaDeviceProp prop;
  CUCHECK(cudaGetDeviceProperties(&prop, 0));
  printf("GPU %d: %s\n", 0, prop.name);

  // Make a copy of the data for the GPU to use.
  // Allocate page-locked virtual memory (that won't be moved from its
  // position in physical memory) so the data can be copied to the GPU
  // via DMA This approximately double the throughput.  Just be sure
  // to free the data with cudaFreeHost() rather than delete[].
  CUCHECK(cudaMallocHost((void**)&plmemory, width*width*sizeof(double)));
  memcpy(plmemory, input, sizeof(double)*width*width);

  // run the GPU version of the algorithm
  if (blockSize == -1) blockSize = getBestThreadBlockSize(width);

  elapsed = haar_not_lifting_2d_cuda(width, plmemory, inverse, steps,
                                     blockSize, true);

  memcpy(output, plmemory, sizeof(double)*width*width);

  printf("CUDA: %.6f ms\n", elapsed);

  CUCHECK(cudaFreeHost(plmemory));

  return 0;
}

