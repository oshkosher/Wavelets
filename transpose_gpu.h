#ifndef __TRANSPOSE_GPU_H__
#define __TRANSPOSE_GPU_H__

#include "wavelet.h"

// Transpose the upper left square corner of the 2-d array.

// These are implemented via templates, so these wrappers just instatiate
// the template for floats and doubles.
void gpuTransposeSquare(int fullWidth, int transposeSize,
                        float *matrix_d, float *matrixTx_d,
                        cudaStream_t stream);

// double support was added in version 1.3
#if !defined(__CUDA_ARCH__) || (__CUDA_ARCH__ >= 130)
void gpuTransposeSquare(int fullWidth, int transposeSize,
                        double *matrix_d, double *matrixTx_d,
                        cudaStream_t stream);
#endif

void gpuTranspose(float *dest, float *src, int width, int height, 
                  cudaStream_t stream = 0);

// 'size' will be modified to match the destination shape
void gpuTranspose3dFwd(float *dest, float *src, scu_wavelet::int3 &size,
                       cudaStream_t stream = 0);

void gpuTranspose3dBack(float *dest, float *src, scu_wavelet::int3 &size,
                        cudaStream_t stream = 0);


#endif // __TRANSPOSE_GPU_H__
