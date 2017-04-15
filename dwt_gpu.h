#ifndef __DWT_GPU_H__
#define __DWT_GPU_H__

#include "cuda.h"
#include "cuda_timer.h"

/*
  CUDA implementation of Haar discrete wavelet transform.

  Ed Karrels, ed.karrels@gmail.com, June 2014
*/

class WaveletAtomic {

 public:

  // I borrowed this code from CUDA/lloyds/cudalloyds.cu, thanks David :-)
  __device__ static float max(float* address, float val) {
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
      assumed = old;
      old = ::atomicCAS(address_as_i, assumed,
                        __float_as_int(::fmaxf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
  }

  __device__ static float min(float* address, float val) {
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
      assumed = old;
      old = ::atomicCAS(address_as_i, assumed,
                        __float_as_int(::fminf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
  }
};


/*
  Forward transform on a square 2-d array.

  arrayWidth is the length of a row in the data.
  transformLength is the size of the upper-left part of the array
  that will be transformed.  For example, in the third pass on a
  1024x1024 array, arrayWidth would be 1024 and transformLength would
  be 256, and array elements [0..255][0..255] would be modified.
*/
template<typename NUM>
__global__ void haar_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *temp);

// Same as haar_not_lifting_2d_kernel, but the inverse transform
template<typename NUM>
__global__ void haar_inv_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *temp);

// Wrapper functions that call the CUDA functions above.
// Even though the functions above can be called with NUM as any type,
// use these wrappers to make it difficult to use them for anything other
// than floats or doubles, since those are they only things that have
// been tested.
float haar_2d_cuda
(int size, float *data, bool inverse = false, int stepCount = -1,
 int threadBlockSize = 128, bool useCombinedTranspose = true);

// double support was added in version 1.3
#if !defined(__CUDA_ARCH__) || (__CUDA_ARCH__ >= 130)
float haar_2d_cuda
(int size, double *data, bool inverse = false, int stepCount = -1,
 int threadBlockSize = 128, bool useCombinedTranspose = true);
#endif

template<typename NUM>
__global__ void haar_transpose_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result,
 int tileSize);

template<typename NUM>
__global__ void haar_inv_transpose_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result, int tileSize);

// transform the data and update 'size', since the dimensions will rotate
void haar_3d_cuda(float *data, float *tmpData,
                  scu_wavelet::int3 &size, scu_wavelet::int3 stepCount,
                  bool inverse = false,
                  CudaTimer *transformTimer = NULL,
                  CudaTimer *transposeTimer = NULL);

void cdf97_3d_cuda(float *data, float *tmpData,
                   scu_wavelet::int3 &size, scu_wavelet::int3 stepCount,
                   bool inverse = false,
                   CudaTimer *transformTimer = NULL,
                   CudaTimer *transposeTimer = NULL);

int bestHaarGPUTileSize();

void haar_v2(float *data_in, float *data_tmp,
             scu_wavelet::int3 size, int level_count,
             CudaTimer *transformTimer = NULL);

void haar_3d_cuda_v2(float *data, float *data_tmp, scu_wavelet::int3 &size,
                     scu_wavelet::int3 stepCount, bool inverse,
                     CudaTimer *transformTimer, CudaTimer *transposeTimer);

void cdf97_v2(float *data, float *data_tmp,
              scu_wavelet::int3 size, int level_count,
              CudaTimer *transformTimer);

void cdf97_3d_cuda_v2(float *data, float *tmpData,
                      scu_wavelet::int3 &size, scu_wavelet::int3 stepCount,
                      bool inverse,
                      CudaTimer *transformTimer, CudaTimer *transposeTimer);

void cdf97_3d_cuda_v3(float *data, float *tmpData,
                      scu_wavelet::int3 &size, scu_wavelet::int3 stepCount,
                      bool inverse,
                      CudaTimer *transformTimer, CudaTimer *transposeTimer);

#endif // __DWT_GPU_H__
