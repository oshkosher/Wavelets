#ifndef __TRANSPOSE_GPU_H__
#define __TRANSPOSE_GPU_H__

// Transpose the upper left square corner of the 2-d array.
void gpuTranspose(int fullWidth, int transposeSize,
                  float *matrix_d, float *matrixTx_d,
                  cudaStream_t stream = 0);

__global__ void gpuTransposeKernel(int fullWidth, int transposeSize,
                                   float *matrix, float *matrixTx);

#endif // __TRANSPOSE_GPU_H__
