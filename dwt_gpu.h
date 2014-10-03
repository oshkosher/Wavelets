#ifndef __DWT_GPU_H__
#define __DWT_GPU_H__

/*
  CUDA implementation of Haar discrete wavelet transform.

  Ed Karrels, ed.karrels@gmail.com, June 2014
*/

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


#endif // __DWT_GPU_H__
