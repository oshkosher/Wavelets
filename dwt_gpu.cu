#include "dwt_cpu.h"
#include "dwt_gpu.h"
#include "transpose_gpu.h"
#include "cucheck.h"
#include "nixtimer.h"

#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f

surface<void, cudaSurfaceType2D> surfRef;

template<typename NUM>
static float haar_not_lifting_2d_cuda_internal
  (int size, NUM *data, bool inverse, int stepCount, int threadBlockSize);

// Returns the time in milliseconds between ev1 and ev2 (where ev1
// happened first.
static float elapsed(cudaEvent_t ev1, cudaEvent_t ev2) {
  float ms;
  CUCHECK(cudaEventElapsedTime(&ms, ev1, ev2));
  return ms;
}


/*
  This does a Haar discrete wavelet transform on each row of
  a 2-d array. Each thread block processes one row.
  This version does not use lifting, and all data is in global memory.
  Input data is in data[], results will be in result[].
*/
template<typename NUM>
__global__ void haar_not_lifting_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result) {

  // each thread block processes one row of data
  int y = blockIdx.x;

  // make pointers to my row of data
  NUM *inputRow = data + y * arrayWidth;
  NUM *outputRow = result + y * arrayWidth;

  // Set s to point to my row in the output data
  NUM *s = outputRow;

  int half = transformLength >> 1;
  
  // point d at the second half of the temporary row
  NUM *d = s + half;
  
  for (int i=threadIdx.x; i < half; i += blockDim.x) {
    NUM a = inputRow[2*i], b = inputRow[2*i + 1];
    d[i] = (a - b) * INV_SQRT2;
    s[i] = (a + b) * INV_SQRT2;
  }
}


__global__ void haar_not_lifting_2d_surfacedata_kernel
(int arrayWidth, int transformLength, float *temp,
 bool transposed) {
  // each thread block processes one row of data
  int y = blockIdx.x;

  // Using shared memory slows this down for some reason,
  // from 144ms to 460ms
  // __shared__ float tempRow[8192];

  // adjust 'data' and 'temp' so they point at my row
  float *tempRow  = temp + y * arrayWidth;

  // Set s to point to my row in the temporary data
  float *s = tempRow;

  int half = transformLength >> 1;
  
  // point d at the second half of the temporary row
  float *d = s + half;
  
  for (int i=threadIdx.x; i < half; i += blockDim.x) {
    float a, b;
    if (transposed) {
      surf2Dread(&a, surfRef, y*sizeof(float), 2*i);
      surf2Dread(&b, surfRef, y*sizeof(float), (2*i+1));
    } else {
      surf2Dread(&a, surfRef, 2*i*sizeof(float), y);
      surf2Dread(&b, surfRef, (2*i+1)*sizeof(float), y);
    }

    float diff = (a - b) * INV_SQRT2;
    float sum = (a + b) * INV_SQRT2;

    d[i] = diff;
    s[i] = sum;
  }

  // sync before other threads read from s[i] and d[i]
  __syncthreads();

  // copy the results back to the surface
  for (int i=threadIdx.x; i < transformLength; i += blockDim.x) {
    if (transposed)
      surf2Dwrite(tempRow[i], surfRef, y * sizeof(float), i);
    else
      surf2Dwrite(tempRow[i], surfRef, i * sizeof(float), y);
  }
}


/* Inverse Haar wavelet transform. */
template<typename NUM>
__global__ void haar_inv_not_lifting_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result) {

  // each thread block processes one row of data
  int y = blockIdx.x;

  // make pointers to my row of data
  NUM *inputRow = data + y * arrayWidth;
  NUM *outputRow = result + y * arrayWidth;

  // Set s to point to my row in the input data
  NUM *s = inputRow;

  int half = transformLength >> 1;

  // point d at the second half of the temporary row
  NUM *d = s + half;

  for (int i=threadIdx.x; i < half; i += blockDim.x) {
    outputRow[2*i]   = INV_SQRT2 * (s[i] + d[i]);
    outputRow[2*i+1] = INV_SQRT2 * (s[i] - d[i]);
  }
}


float haar_not_lifting_2d_cuda
  (int size, float *data, bool inverse, int stepCount, int threadBlockSize) {
  return haar_not_lifting_2d_cuda_internal
    (size, data, inverse, stepCount, threadBlockSize);
}     


float haar_not_lifting_2d_cuda
  (int size, double *data, bool inverse, int stepCount, int threadBlockSize) {
  return haar_not_lifting_2d_cuda_internal
    (size, data, inverse, stepCount, threadBlockSize);
}     


// Wrapper function that handles the CUDA details.
template<typename NUM>
float haar_not_lifting_2d_cuda_internal
  (int size, NUM *data, bool inverse, int stepCount, int threadBlockSize) {

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
  NUM *data1_dev, *data2_dev;
  size_t totalBytes = size * size * sizeof(NUM);
  CUCHECK(cudaMalloc((void**) &data1_dev, totalBytes));
  CUCHECK(cudaMalloc((void**) &data2_dev, totalBytes));

  // Create a stream to enable asynchronous operation, to minimize
  // time between kernel calls.
  cudaStream_t stream = 0;
  CUCHECK(cudaStreamCreate(&stream));

  // start the timer
  double startTimeCPU = NixTimer::time();
  CUCHECK(cudaEventRecord(eventStart, stream));

  // copy the data to the GPU
  CUCHECK(cudaEventRecord(eventCopyToStart, stream));
  CUCHECK(cudaMemcpyAsync(data1_dev, data, totalBytes, cudaMemcpyHostToDevice,
                          stream));
  CUCHECK(cudaEventRecord(eventCopyToEnd, stream));

  int transformLength;

  if (inverse) {

    // inverse
    transformLength = size >> (stepCount - 1);
    for (int i=0; i < stepCount; i++) {

      // transpose the matrix into temp_dev
      CUCHECK(cudaEventRecord(transposeEvents[i*4], stream));
      gpuTranspose(size, transformLength, data1_dev, data2_dev, stream);
      CUCHECK(cudaEventRecord(transposeEvents[i*4+1], stream));

      // transform columns in temp_dev
      CUCHECK(cudaEventRecord(transformEvents[i*4], stream));
      haar_inv_not_lifting_2d_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, data2_dev, data1_dev);
      CUCHECK(cudaEventRecord(transformEvents[i*4+1], stream));

      // transpose the matrix into data_dev
      CUCHECK(cudaEventRecord(transposeEvents[i*4+2], stream));
      gpuTranspose(size, transformLength, data1_dev, data2_dev, stream);
      CUCHECK(cudaEventRecord(transposeEvents[i*4+3], stream));
    
      // transform rows in data_dev
      CUCHECK(cudaEventRecord(transformEvents[i*4+2], stream));
      haar_inv_not_lifting_2d_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, data2_dev, data1_dev);
      CUCHECK(cudaEventRecord(transformEvents[i*4+3], stream));

      // results are in data1_dev

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
        (size, transformLength, data1_dev, data2_dev);
      CUCHECK(cudaEventRecord(transformEvents[i*4+1], stream));

      // transpose the matrix into temp_dev
      CUCHECK(cudaEventRecord(transposeEvents[i*4], stream));
      gpuTranspose(size, transformLength, data2_dev, data1_dev, stream);
      CUCHECK(cudaEventRecord(transposeEvents[i*4+1], stream));
    
      // do the wavelet transform on columns
      CUCHECK(cudaEventRecord(transformEvents[i*4+2], stream));
      haar_not_lifting_2d_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, data1_dev, data2_dev);
      CUCHECK(cudaEventRecord(transformEvents[i*4+3], stream));
    
      // transpose the matrix back into data_dev
      CUCHECK(cudaEventRecord(transposeEvents[i*4+2], stream));
      gpuTranspose(size, transformLength, data2_dev, data1_dev, stream);
      CUCHECK(cudaEventRecord(transposeEvents[i*4+3], stream));

      // results are in data1_dev

      transformLength >>= 1;
    }

  }

  // copy the data back from the GPU
  CUCHECK(cudaEventRecord(eventCopyFromStart, stream));
  CUCHECK(cudaMemcpyAsync(data, data1_dev, totalBytes, cudaMemcpyDeviceToHost,
                          stream));
  CUCHECK(cudaEventRecord(eventCopyFromEnd, stream));

  // Since all the GPU tasks were started asynchronously, control should
  // flow to this point very quickly. The cudaEventSynchronize() call will
  // wait until the GPU is finished.
  double endTimeCPU = NixTimer::time();
  printf("Time elapsed creating GPU tasks: %.6f ms\n",
         1000*(endTimeCPU - startTimeCPU));
  fflush(stdout);

  // stop the timer
  CUCHECK(cudaEventRecord(eventEnd, stream));
  CUCHECK(cudaEventSynchronize(eventEnd));
  cudaStreamDestroy(stream);

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
  CUCHECK(cudaFree(data1_dev));
  CUCHECK(cudaFree(data2_dev));

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


// like haar_not_lifting_2d_cuda, but using surfaces
float haar_not_lifting_2d_cuda_surfaces
(int size, float *data, bool inverse, int stepCount, int threadBlockSize) {

  int maxSteps = dwtMaximumSteps(size);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  // create timers
  cudaEvent_t eventStart, eventEnd, eventCopyToStart, eventCopyToEnd,
    eventCopyFromStart, eventCopyFromEnd;
  cudaEvent_t *transformEvents;
  transformEvents = new cudaEvent_t[stepCount*4];

  CUCHECK(cudaEventCreate(&eventStart));
  CUCHECK(cudaEventCreate(&eventEnd));
  CUCHECK(cudaEventCreate(&eventCopyToStart));
  CUCHECK(cudaEventCreate(&eventCopyToEnd));
  CUCHECK(cudaEventCreate(&eventCopyFromStart));
  CUCHECK(cudaEventCreate(&eventCopyFromEnd));

  for (int i=0; i < stepCount*4; i++) {
    CUCHECK(cudaEventCreate(transformEvents+i));
  }

  // allocate memory for the data and the temp space on the GPU
  float *data_dev = 0, *temp_dev = 0;
  size_t totalBytes = size * size * sizeof(float);
  // CUCHECK(cudaMalloc((void**) &data_dev, totalBytes));
  CUCHECK(cudaMalloc((void**) &temp_dev, totalBytes));

  // try out arrays
  cudaChannelFormatDesc texFormat = cudaCreateChannelDesc<float>();
  cudaArray_t cuArray;
  CUCHECK(cudaMallocArray(&cuArray, &texFormat, size, size,
                          cudaArraySurfaceLoadStore));
  CUCHECK(cudaBindSurfaceToArray(surfRef, cuArray));

  // Create a stream to enable asynchronous operation, to minimize
  // time between kernel calls.
  cudaStream_t stream = 0;
  CUCHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));

  // start the timer
  double startTimeCPU = NixTimer::time();
  CUCHECK(cudaEventRecord(eventStart, stream));

  // copy the data to the GPU
  CUCHECK(cudaEventRecord(eventCopyToStart, stream));
  CUCHECK(cudaMemcpy2DToArrayAsync(cuArray, 0, 0, data, size*sizeof(float),
                                   size * sizeof(float), size, 
                                   cudaMemcpyHostToDevice, stream));
  CUCHECK(cudaEventRecord(eventCopyToEnd, stream));

  int transformLength;

  if (inverse) {

    // inverse
    transformLength = size >> (stepCount - 1);
    for (int i=0; i < stepCount; i++) {

      // transform columns in temp_dev
      CUCHECK(cudaEventRecord(transformEvents[i*4], stream));
      /*
      haar_inv_not_lifting_2d_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, temp_dev, data_dev);
      */
      CUCHECK(cudaEventRecord(transformEvents[i*4+1], stream));
    
      // transform rows in data_dev
      CUCHECK(cudaEventRecord(transformEvents[i*4+2], stream));
      /*
      haar_inv_not_lifting_2d_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, data_dev, temp_dev);
      */
      CUCHECK(cudaEventRecord(transformEvents[i*4+3], stream));

      transformLength <<= 1;
    }

  } else {

    // forward
    transformLength = size;
    for (int i=0; i < stepCount; i++) {
      
      // do the wavelet transform on rows
      CUCHECK(cudaEventRecord(transformEvents[i*4], stream));
      haar_not_lifting_2d_surfacedata_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, temp_dev, false);
      CUCHECK(cudaEventRecord(transformEvents[i*4+1], stream));
    
      // do the wavelet transform on columns
      CUCHECK(cudaEventRecord(transformEvents[i*4+2], stream));
      haar_not_lifting_2d_surfacedata_kernel
        <<<transformLength, threadBlockSize, 0, stream>>>
        (size, transformLength, temp_dev, true);
      CUCHECK(cudaEventRecord(transformEvents[i*4+3], stream));

      transformLength >>= 1;
    }

  }

  // copy the data back from the GPU
  CUCHECK(cudaEventRecord(eventCopyFromStart, stream));
  CUCHECK(cudaMemcpy2DFromArrayAsync(data, size*sizeof(float), cuArray, 0, 0, 
                                     size * sizeof(float), size, 
                                     cudaMemcpyDeviceToHost, stream));
  CUCHECK(cudaEventRecord(eventCopyFromEnd, stream));

  double endTimeCPU = NixTimer::time();
  printf("Time elapsed creating GPU tasks: %.6f ms\n",
         1000*(endTimeCPU - startTimeCPU));
  fflush(stdout);

  // stop the timer
  CUCHECK(cudaEventRecord(eventEnd, stream));
  CUCHECK(cudaEventSynchronize(eventEnd));
  cudaStreamDestroy(stream);

  // check for errors
  CUCHECK(cudaGetLastError());

  printf("Times:\n");
  float transformTime = 0, transposeTime = 0;
  for (int i=0; i < stepCount*2; i++) {
    transformTime += elapsed(transformEvents[i*2], transformEvents[i*2+1]);
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
  }

  return totalTime;
}


/*
  Do Haar DWT on columns rather than rows.
  This is equivalent to:
    transpose
    haar_not_lifting_2d_kernel
    transpose
  This just demonstrates the importance of the memory access pattern.
  It is approximately 20x slower than accessing the data by rows.
    Row transform:         20.382 ms (1 call)
    Column transform:     407.445 ms (1 call)

  Compared to rows-transpose-rows-transpose:
    Transform time:        40.730 ms (2 calls)
    Transpose time:        57.759 ms (2 calls)
 */
template<typename NUM>
__global__ void haar_not_lifting_2d_columns_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result) {

  // each thread block processes one column of data
  int x = blockIdx.x;

  // make pointers to the first element in my column of data
  NUM *inputCol = data + x;
  NUM *outputCol = result + x;

  // Set s to point to my row in the output data
  NUM *s = outputCol;

  int half = transformLength >> 1;
  
  // point d at the second half of the temporary row
  NUM *d = s + half*arrayWidth;
  
  for (int i=threadIdx.x; i < half; i += blockDim.x) {
    NUM a = inputCol[2*i*arrayWidth], b = inputCol[(2*i+1)*arrayWidth];
    /*
    int os = (s+i*arrayWidth) - result;
    int od = (d+i*arrayWidth) - result;
    printf("[%d,%d]  %d.%d %d.%d\n", blockIdx.x, threadIdx.x,
           os/arrayWidth, os%arrayWidth,
           od/arrayWidth, od%arrayWidth);
    */
    d[i*arrayWidth] = (a - b) * INV_SQRT2;
    s[i*arrayWidth] = (a + b) * INV_SQRT2;
  }
}
