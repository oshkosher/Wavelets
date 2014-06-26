#ifndef __DWT_GPU_H__
#define __DWT_GPU_H__

/*
  CUDA implementation of Haar discrete wavelet transform.
*/

// CUDA code for forward transform
__global__ void haar_not_lifting_2d_kernel
  (int width, int height, float *data, float *temp);

// CUDA code for inverse transform
__global__ void haar_inv_not_lifting_2d_kernel
  (int width, int height, float *data, float *temp);

// Wrapper function that calls the CUDA functions above.
float haar_not_lifting_2d_cuda
  (int width, int height, float *data, bool inverse = false);




#endif // __DWT_GPU_H__
