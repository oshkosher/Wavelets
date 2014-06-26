#ifndef __TRANSPOSE_GPU_H__
#define __TRANSPOSE_GPU_H__

void gpuTranspose(int width, int height, float *matrix_d, float *matrixTx_d,
                  cudaStream_t stream = 0);

__global__ void gpuTransposeKernel(int width, int height, float *matrix,
                                   float *matrixTx);

#endif // __TRANSPOSE_GPU_H__
