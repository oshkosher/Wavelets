#ifndef __TRANSPOSE_GPU_H__
#define __TRANSPOSE_GPU_H__

// Transpose the upper left square corner of the 2-d array.

// These are implemented via templates, so these wrappers just instatiate
// the template for floats and doubles.
void gpuTranspose(int fullWidth, int transposeSize,
                  float *matrix_d, float *matrixTx_d,
                  cudaStream_t stream);

void gpuTranspose(int fullWidth, int transposeSize,
                  double *matrix_d, double *matrixTx_d,
                  cudaStream_t stream);

#endif // __TRANSPOSE_GPU_H__
