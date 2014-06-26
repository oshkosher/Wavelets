#include "transpose_gpu.h"
#include "cucheck.h"

#define BLOCK_SIZE 128
#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f
#define DO_NORMALIZE 0

/*
  This does a Haar discrete wavelet transform on each row of
  a 2-d array. Each thread block processes one row.
  This version does not use lifting, and all data is in global memory.
*/
__global__ void haar_not_lifting_2d_kernel
(int width, int height, float *data, float *temp) {

  // each thread block processes one row of data
  int y = blockIdx.x;

  // adjust 'data' and 'temp' so they point at my row
  data += y * width;
  temp += y * width;

  // Set s to point to my row in the temporary data
  float *s = temp;

  int sampleCount = width;

  while (sampleCount > 1) {
    int half = sampleCount >> 1;

    // point d at the second half of the temporary row
    float *d = s + half;

    for (int i=threadIdx.x; i < half; i += blockDim.x) {
      // possible optimization: read these with a float2
      float a = data[2*i], b = data[2*i+1];

      d[i] = b - a;
      s[i] = a + .5f * d[i];
#if DO_NORMALIZE
      s[i] *= SQRT2;
      d[i] *= INV_SQRT2;
#endif
    }

    // sync before other threads read from s[i] and d[i]
    __syncthreads();

    // copy the results back to data[]
    for (int i=threadIdx.x; i < sampleCount; i += blockDim.x)
      data[i] = temp[i];
    
    // sync before the next round of reading
    __syncthreads();

    sampleCount = half;
  }

}


/* Inverse Haar wavelet transform. */
__global__ void haar_inv_not_lifting_2d_kernel
(int width, int height, float *data, float *temp) {

  // each thread block processes one row of data
  int y = blockIdx.x;

  // adjust 'data' and 'temp' so they point at my row
  data += y * width;
  temp += y * width;

  // Set s to point to my row in the temporary data
  float *s = data;

  int sampleCount = 2;

  while (sampleCount <= width) {
    int half = sampleCount >> 1;

    // point d at the second half of the temporary row
    float *d = s + half;

    for (int i=threadIdx.x; i < half; i += blockDim.x) {
#if DO_NORMALIZE
      s[i] *= INV_SQRT2;
      d[i] *= SQRT2;
#endif

      temp[2*i]     = s[i] - .5f * d[i];
      temp[2*i + 1] = temp[2*i] + d[i];
    }

    // sync before other threads read from s[i] and d[i]
    __syncthreads();

    // copy the results back to data[]
    for (int i=threadIdx.x; i < sampleCount; i += blockDim.x)
      data[i] = temp[i];
    
    // sync before the next round of reading
    __syncthreads();

    sampleCount <<= 1;
  }

}


// Wrapper function that handles the CUDA details.
float haar_not_lifting_2d_cuda
(int width, int height, float *data, bool inverse) {

  // create a timer
  cudaEvent_t startEvent, stopEvent;
  float elapsedTimeMs;
  CUCHECK(cudaEventCreate(&startEvent));
  CUCHECK(cudaEventCreate(&stopEvent));

  // allocate memory for the data and the temp space on the GPU
  float *data_dev, *temp_dev;
  size_t totalBytes = width * height * sizeof(float);
  CUCHECK(cudaMalloc((void**) &data_dev, totalBytes));
  CUCHECK(cudaMalloc((void**) &temp_dev, totalBytes));

  // Create a stream to enable asynchronous operation, to minimize
  // time between kernel calls.
  // cudaStream_t stream;
  // CUCHECK(cudaStreamCreate(&stream));

  // copy the data to the GPU
  CUCHECK(cudaMemcpy(data_dev, data, totalBytes, cudaMemcpyHostToDevice));

  // start the timer
  CUCHECK(cudaEventRecord(startEvent, 0));

  if (inverse) {
    gpuTranspose(width, height, data_dev, temp_dev);

    haar_inv_not_lifting_2d_kernel<<<height, BLOCK_SIZE>>>
      (height, width, temp_dev, data_dev);

    gpuTranspose(height, width, temp_dev, data_dev);
    
    haar_inv_not_lifting_2d_kernel<<<width, BLOCK_SIZE>>>
      (width, height, data_dev, temp_dev);

  } else {
    
    // do the wavelet transform on rows
    haar_not_lifting_2d_kernel<<<height, BLOCK_SIZE>>>
      (width, height, data_dev, temp_dev);
    
    // transpose the matrix into temp_dev
    gpuTranspose(width, height, data_dev, temp_dev);
    
    // do the wavelet transform on columns
    haar_not_lifting_2d_kernel<<<width, BLOCK_SIZE>>>
      (height, width, temp_dev, data_dev);
    
    // transpose the matrix back into data_dev
    gpuTranspose(height, width, temp_dev, data_dev);
  }

  // stop the timer
  CUCHECK(cudaEventRecord(stopEvent, 0));
  CUCHECK(cudaEventSynchronize(stopEvent));
  CUCHECK(cudaEventElapsedTime(&elapsedTimeMs, startEvent, stopEvent));

  // copy the data back from the GPU
  CUCHECK(cudaMemcpy(data, data_dev, totalBytes, cudaMemcpyDeviceToHost));

  // check for errors
  CUCHECK(cudaGetLastError());

  // deallocate GPU memory
  CUCHECK(cudaFree(data_dev));
  CUCHECK(cudaFree(temp_dev));

  return elapsedTimeMs;
}
