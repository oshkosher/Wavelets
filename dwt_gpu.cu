#include "dwt_cpu.h"
#include "dwt_gpu.h"
#include "transpose_gpu.h"
#include "cucheck.h"
#include "nixtimer.h"

#define BLOCK_SIZE 32
#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f

/*
  This does a Haar discrete wavelet transform on each row of
  a 2-d array. Each thread block processes one row.
  This version does not use lifting, and all data is in global memory.
*/
__global__ void haar_not_lifting_2d_kernel
(int arrayWidth, int transformLength, float *data, float *temp) {

  // each thread block processes one row of data
  int y = blockIdx.x;

  // adjust 'data' and 'temp' so they point at my row
  float *inputRow = data + y * arrayWidth;
  float *tempRow  = temp + y * arrayWidth;

  // if (transformLength < arrayWidth)
  // printf("[%d,%d] row data[%d]\n", blockIdx.x, threadIdx.x, inputRow - data);

  // Set s to point to my row in the temporary data
  float *s = tempRow;

  int half = transformLength >> 1;
  
  // point d at the second half of the temporary row
  float *d = s + half;
  
  for (int i=threadIdx.x; i < half; i += blockDim.x) {
    // possible optimization: read these with a float2
    float a = inputRow[2*i], b = inputRow[2*i + 1];
    d[i] = (a - b) * INV_SQRT2;
    s[i] = (a + b) * INV_SQRT2;
  }

  // sync before other threads read from s[i] and d[i]
  __syncthreads();

  // copy the results back to data[]
  for (int i=threadIdx.x; i < transformLength; i += blockDim.x)
    inputRow[i] = tempRow[i];
}


/* Inverse Haar wavelet transform. */
__global__ void haar_inv_not_lifting_2d_kernel
(int arrayWidth, int transformLength, float *data, float *temp) {

  // each thread block processes one row of data
  int y = blockIdx.x;

  // adjust 'data' and 'temp' so they point at my row
  data += y * arrayWidth;
  temp += y * arrayWidth;

  // Set s to point to my row in the temporary data
  float *s = data;

  int half = transformLength >> 1;

  // point d at the second half of the temporary row
  float *d = s + half;

  for (int i=threadIdx.x; i < half; i += blockDim.x) {
    temp[2*i]   = INV_SQRT2 * (s[i] + d[i]);
    temp[2*i+1] = INV_SQRT2 * (s[i] - d[i]);
  }

  // sync before other threads read from s[i] and d[i]
  __syncthreads();

  // copy the results back to data[]
  for (int i=threadIdx.x; i < transformLength; i += blockDim.x)
    data[i] = temp[i];
}


float elapsed(cudaEvent_t ev1, cudaEvent_t ev2) {
  float ms;
  CUCHECK(cudaEventElapsedTime(&ms, ev1, ev2));
  return ms;
}


// Wrapper function that handles the CUDA details.
float haar_not_lifting_2d_cuda
(int size, float *data, bool inverse, int stepCount) {

  float *dummy;
  CUCHECK(cudaMalloc((void**) &dummy, 100));

  int maxSteps = dwtMaximumSteps(size);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  // create timers
  cudaEvent_t startEvent, eventCopyTo, eventCopyFrom,
    eventTransform1, eventTransform2,
    eventTranspose1, eventTranspose2;

  CUCHECK(cudaEventCreate(&startEvent));
  CUCHECK(cudaEventCreate(&eventCopyTo));
  CUCHECK(cudaEventCreate(&eventTransform1));
  CUCHECK(cudaEventCreate(&eventTransform2));
  CUCHECK(cudaEventCreate(&eventTranspose1));
  CUCHECK(cudaEventCreate(&eventTranspose2));
  CUCHECK(cudaEventCreate(&eventCopyFrom));

  // allocate memory for the data and the temp space on the GPU
  float *data_dev, *temp_dev;
  size_t totalBytes = size * size * sizeof(float);
  CUCHECK(cudaMalloc((void**) &data_dev, totalBytes));
  CUCHECK(cudaMalloc((void**) &temp_dev, totalBytes));

  // Create a stream to enable asynchronous operation, to minimize
  // time between kernel calls.
  cudaStream_t stream = 0;
  CUCHECK(cudaStreamCreate(&stream));

  // start the timer
  double startTimeCPU = NixTimer::time();
  CUCHECK(cudaEventRecord(startEvent, stream));

  // copy the data to the GPU
  // CUCHECK(cudaMemcpy(data_dev, data, totalBytes, cudaMemcpyHostToDevice));
  CUCHECK(cudaMemcpyAsync(data_dev, data, totalBytes, cudaMemcpyHostToDevice,
                          stream));
  CUCHECK(cudaEventRecord(eventCopyTo, stream));

  int transformLength;

  if (inverse) {

    // inverse
    transformLength = 2;
    for (int i=0; i < stepCount; i++) {

      gpuTranspose(size, transformLength, data_dev, temp_dev, stream);
      CUCHECK(cudaEventRecord(eventTranspose1, stream));

      haar_inv_not_lifting_2d_kernel<<<transformLength, BLOCK_SIZE, 0, stream>>>
        (size, transformLength, temp_dev, data_dev);
      CUCHECK(cudaEventRecord(eventTransform1, stream));

      gpuTranspose(size, transformLength, temp_dev, data_dev, stream);
      CUCHECK(cudaEventRecord(eventTranspose2, stream));
    
      haar_inv_not_lifting_2d_kernel<<<transformLength, BLOCK_SIZE, 0, stream>>>
        (size, transformLength, data_dev, temp_dev);
      CUCHECK(cudaEventRecord(eventTransform2, stream));

      transformLength <<= 1;
    }

  } else {

    // forward
    transformLength = size;
    for (int i=0; i < stepCount; i++) {

      printf("%d of %d\n", transformLength, size);
      
      // do the wavelet transform on rows
      haar_not_lifting_2d_kernel<<<transformLength, BLOCK_SIZE, 0, stream>>>
        (size, transformLength, data_dev, temp_dev);
      CUCHECK(cudaEventRecord(eventTransform1, stream));
    
      // transpose the matrix into temp_dev
      gpuTranspose(size, transformLength, data_dev, temp_dev, stream);
      CUCHECK(cudaEventRecord(eventTranspose1, stream));
    
      // do the wavelet transform on columns
      haar_not_lifting_2d_kernel<<<transformLength, BLOCK_SIZE, 0, stream>>>
        (size, transformLength, temp_dev, data_dev);
      CUCHECK(cudaEventRecord(eventTransform2, stream));
    
      // transpose the matrix back into data_dev
      gpuTranspose(size, transformLength, temp_dev, data_dev, stream);
      CUCHECK(cudaEventRecord(eventTranspose2, stream));

      transformLength >>= 1;
    }

  }

  // copy the data back from the GPU
  CUCHECK(cudaMemcpyAsync(data, data_dev, totalBytes, cudaMemcpyDeviceToHost,
                          stream));
  double endTimeCPU = NixTimer::time();
  printf("Time elapsed during GPU processing: %.6f ms\n",
         1000*(endTimeCPU - startTimeCPU));

  // stop the timer
  CUCHECK(cudaEventRecord(eventCopyFrom, stream));
  CUCHECK(cudaEventSynchronize(eventCopyFrom));

  // check for errors
  CUCHECK(cudaGetLastError());

  printf("Times:\n");
  printf("  Copy data to GPU: %.6f ms\n", elapsed(startEvent, eventCopyTo));
  cudaEvent_t lastEvent;
  if (inverse) {
    printf("  Transpose: %.6f ms\n", elapsed(eventCopyTo, eventTranspose1));
    printf("  Transform: %.6f ms\n", elapsed(eventTranspose1, eventTransform1));
    printf("  Transpose: %.6f ms\n", elapsed(eventTransform1, eventTranspose2));
    printf("  Transform: %.6f ms\n", elapsed(eventTranspose2, eventTransform2));
    lastEvent = eventTransform2;
  } else {
    printf("  Transform: %.6f ms\n", elapsed(eventCopyTo, eventTransform1));
    printf("  Transpose: %.6f ms\n", elapsed(eventTransform1, eventTranspose1));
    printf("  Transform: %.6f ms\n", elapsed(eventTranspose1, eventTransform2));
    printf("  Transpose: %.6f ms\n", elapsed(eventTransform2, eventTranspose2));
    lastEvent = eventTranspose2;
  }
  printf("  Copy data from GPU: %.6f ms\n", elapsed(lastEvent, eventCopyFrom));

  // deallocate GPU memory
  CUCHECK(cudaFree(data_dev));
  CUCHECK(cudaFree(temp_dev));

  float totalTime = elapsed(startEvent, eventCopyFrom);

  // deallocate timer event
  CUCHECK(cudaEventDestroy(startEvent));
  CUCHECK(cudaEventDestroy(eventCopyTo));
  CUCHECK(cudaEventDestroy(eventTransform1));
  CUCHECK(cudaEventDestroy(eventTransform2));
  CUCHECK(cudaEventDestroy(eventTranspose1));
  CUCHECK(cudaEventDestroy(eventTranspose2));
  CUCHECK(cudaEventDestroy(eventCopyFrom));

  return totalTime;
}
