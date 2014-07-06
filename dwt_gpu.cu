#include "dwt_cpu.h"
#include "dwt_gpu.h"
#include "transpose_gpu.h"
#include "cucheck.h"
#include "nixtimer.h"

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
(int size, float *data, bool inverse, int stepCount, int threadBlockSize) {

  int maxSteps = dwtMaximumSteps(size);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  // create timers
  cudaEvent_t eventStart, eventEnd, eventCopyToStart, eventCopyToEnd,
    eventCopyFromStart, eventCopyFromEnd;
  cudaEvent_t *transformEvents, *transposeEvents;
  transformEvents = new cudaEvent_t[stepCount*4];
  transposeEvents = new cudaEvent_t[stepCount*4];

  CUCHECK(cudaEventCreate(&eventStart));
  CUCHECK(cudaEventCreate(&eventEnd));
  CUCHECK(cudaEventCreate(&eventCopyToStart));
  CUCHECK(cudaEventCreate(&eventCopyToEnd));
  CUCHECK(cudaEventCreate(&eventCopyFromStart));
  CUCHECK(cudaEventCreate(&eventCopyFromEnd));

  for (int i=0; i < stepCount*4; i++) {
    CUCHECK(cudaEventCreate(transformEvents+i));
    CUCHECK(cudaEventCreate(transposeEvents+i));
  }

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
  // double startTimeCPU = NixTimer::time();
  CUCHECK(cudaEventRecord(eventStart, stream));

  // copy the data to the GPU
  // CUCHECK(cudaMemcpy(data_dev, data, totalBytes, cudaMemcpyHostToDevice));
  CUCHECK(cudaEventRecord(eventCopyToStart, stream));
  CUCHECK(cudaMemcpyAsync(data_dev, data, totalBytes, cudaMemcpyHostToDevice,
                          stream));
  CUCHECK(cudaEventRecord(eventCopyToEnd, stream));

  int transformLength;

  if (inverse) {

    // inverse
    transformLength = size >> (stepCount - 1);
    for (int i=0; i < stepCount; i++) {

      // transpose the matrix into temp_dev
      CUCHECK(cudaEventRecord(transposeEvents[i*4], stream));
      gpuTranspose(size, transformLength, data_dev, temp_dev, stream);
      CUCHECK(cudaEventRecord(transposeEvents[i*4+1], stream));

      // transform columns in temp_dev
      CUCHECK(cudaEventRecord(transformEvents[i*4], stream));
      haar_inv_not_lifting_2d_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, temp_dev, data_dev);
      CUCHECK(cudaEventRecord(transformEvents[i*4+1], stream));

      // transpose the matrix into data_dev
      CUCHECK(cudaEventRecord(transposeEvents[i*4+2], stream));
      gpuTranspose(size, transformLength, temp_dev, data_dev, stream);
      CUCHECK(cudaEventRecord(transposeEvents[i*4+3], stream));
    
      // transform rows in data_dev
      CUCHECK(cudaEventRecord(transformEvents[i*4+2], stream));
      haar_inv_not_lifting_2d_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, data_dev, temp_dev);
      CUCHECK(cudaEventRecord(transformEvents[i*4+3], stream));

      transformLength <<= 1;
    }

  } else {

    // forward
    transformLength = size;
    for (int i=0; i < stepCount; i++) {
      
      // do the wavelet transform on rows
      CUCHECK(cudaEventRecord(transformEvents[i*4], stream));
      haar_not_lifting_2d_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, data_dev, temp_dev);
      CUCHECK(cudaEventRecord(transformEvents[i*4+1], stream));
    
      // transpose the matrix into temp_dev
      CUCHECK(cudaEventRecord(transposeEvents[i*4], stream));
      gpuTranspose(size, transformLength, data_dev, temp_dev, stream);
      CUCHECK(cudaEventRecord(transposeEvents[i*4+1], stream));
    
      // do the wavelet transform on columns
      CUCHECK(cudaEventRecord(transformEvents[i*4+2], stream));
      haar_not_lifting_2d_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, temp_dev, data_dev);
      CUCHECK(cudaEventRecord(transformEvents[i*4+3], stream));
    
      // transpose the matrix back into data_dev
      CUCHECK(cudaEventRecord(transposeEvents[i*4+2], stream));
      gpuTranspose(size, transformLength, temp_dev, data_dev, stream);
      CUCHECK(cudaEventRecord(transposeEvents[i*4+3], stream));

      transformLength >>= 1;
    }

  }

  // copy the data back from the GPU
  CUCHECK(cudaEventRecord(eventCopyFromStart, stream));
  CUCHECK(cudaMemcpyAsync(data, data_dev, totalBytes, cudaMemcpyDeviceToHost,
                          stream));
  CUCHECK(cudaEventRecord(eventCopyFromEnd, stream));

  // double endTimeCPU = NixTimer::time();
  // printf("Time elapsed during GPU processing: %.6f ms\n",
  // 1000*(endTimeCPU - startTimeCPU));

  // stop the timer
  CUCHECK(cudaEventRecord(eventEnd, stream));
  CUCHECK(cudaEventSynchronize(eventEnd));

  // check for errors
  CUCHECK(cudaGetLastError());

  printf("Times:\n");
  float transformTime = 0, transposeTime = 0;
  for (int i=0; i < stepCount*2; i++) {
    transformTime += elapsed(transformEvents[i*2], transformEvents[i*2+1]);
    transposeTime += elapsed(transposeEvents[i*2], transposeEvents[i*2+1]);
  }

  printf("  Copy data to GPU:   %9.3f ms\n", 
         elapsed(eventCopyToStart, eventCopyToEnd));
  printf("  Transform time:     %9.3f ms (%d calls)\n",
         transformTime, stepCount*2);
  printf("  Transpose time:     %9.3f ms (%d calls)\n", 
         transposeTime, stepCount*2);
  printf("  Copy data from GPU: %9.3f ms\n", 
         elapsed(eventCopyFromStart, eventCopyFromEnd));

  // deallocate GPU memory
  CUCHECK(cudaFree(data_dev));
  CUCHECK(cudaFree(temp_dev));

  float totalTime = elapsed(eventStart, eventEnd);

  // deallocate timer events
  CUCHECK(cudaEventDestroy(eventStart));
  CUCHECK(cudaEventDestroy(eventEnd));
  CUCHECK(cudaEventDestroy(eventCopyToStart));
  CUCHECK(cudaEventDestroy(eventCopyToEnd));
  CUCHECK(cudaEventDestroy(eventCopyFromStart));
  CUCHECK(cudaEventDestroy(eventCopyFromEnd));
  for (int i=0; i < stepCount*4; i++) {
    CUCHECK(cudaEventDestroy(transformEvents[i]));
    CUCHECK(cudaEventDestroy(transposeEvents[i]));
  }

  return totalTime;
}
