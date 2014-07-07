#ifndef __DWT_GPU_H__
#define __DWT_GPU_H__

/*
  CUDA implementation of Haar discrete wavelet transform.
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
__global__ void haar_not_lifting_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *temp);

// Same as haar_not_lifting_2d_kernel, but the inverse transform
template<typename NUM>
__global__ void haar_inv_not_lifting_2d_kernel
  (int arrayWidth, int transformLength, NUM *data, NUM *temp);

// Wrapper functions that call the CUDA functions above.
// Even though the functions above can be called with NUM as any type,
// use these wrappers to make it difficult to use them for anything other
// than floats or doubles, since those are they only things that have
// been tested.
float haar_not_lifting_2d_cuda_float
(int size, float *data, bool inverse = false, int stepCount = -1,
 int threadBlockSize = 128);

float haar_not_lifting_2d_cuda_double
(int size, double *data, bool inverse = false, int stepCount = -1,
 int threadBlockSize = 128);

float haar_not_lifting_2d_cuda_surfaces
(int size, float *data, bool inverse, int stepCount, int threadBlockSize=128);


#endif // __DWT_GPU_H__
