#include "dwt_cpu.h"
#include "dwt_gpu.h"
#include "transpose_gpu.h"
#include "cucheck.h"
#include "nixtimer.h"
#include "cuda_timer.h"

#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f

// width and height of one tile in the Haar-transpose algorithm
// Best setting:
//   GTX 480: 16
//   GTX 570: 16
//   GTX 680: 32
//   K2000M: 16
//   Tesla K20c: 32
#ifndef HT_TILE_SIZE
#define HT_TILE_SIZE 16
#endif

surface<void, cudaSurfaceType2D> surfRef;

template<typename NUM>
float haar_not_lifting_2d_cuda_internal
(int size, NUM *data, bool inverse, int stepCount, int threadBlockSize,
 bool useCombinedTranspose);

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


/*
  Do one pass of a Haar discrete wavelet transform and transpose the
  matrix at the same time.  This splits the data into tiles, computing
  the transform and writing the results into a transposed matrix.

  In these diagrams, the number indicates the thread that reads
  or writes each element. It depicts each thread block as a 4x4 block
  of 16 threads, but in reality it will be larger--16x16 or 32x32.

  Input matrix (global memory)
  +----------------------------------------
  |  0  0  1  1  2  2  3  3  ...next tile...
  |  4  4  5  5  6  6  7  7
  |  8  8  9  9 10 10 11 11
  | 12 12 13 13 14 14 15 15
  | ...next tile...
   
  Shared memory temporary storage
  The results from each "sum" operation are stored in one 2d array,
  and the results form each "difference" operation are stored in a
  different 2d array.

    Write in this order:
     0  1  2  3
     4  5  6  7
     8  9 10 11
    12 13 14 15

    Read in this order:
    0  4  8 12
    1  5  9 13
    2  6 10 14
    3  7 11 15

  Output matrix (global memory)
  +----------------------------------------
  |  0  1  2  3    ...next tile...
  |  4  5  6  7
  |  8  9 10 11
  | 12 13 14 15
  | ...
  | ...
  | ...lower half....
  |  0  1  2  3
  |  4  5  6  7
  |  8  9 10 11
  | 12 13 14 15
  | ...

Success!
Before:
  Transform time:        40.707 ms (2 calls)
  Transpose time:        57.732 ms (2 calls)
After:
  Transform time:        52.512 ms (2 calls)
  Transpose time:         0.004 ms (2 calls)
 */

template<typename NUM>
__global__ void haar_transpose_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result) {
  __shared__ NUM sums [HT_TILE_SIZE][HT_TILE_SIZE+1];
  __shared__ NUM diffs[HT_TILE_SIZE][HT_TILE_SIZE+1];

  int inputx = (blockIdx.x*blockDim.x + threadIdx.x) * 2;
  int inputy = blockIdx.y*blockDim.y + threadIdx.y;
  
  // read a tile 2*HT_TILE_SIZE wide, HT_TILE_SIZE tall, compute
  // the sum and difference coefficients, and store those coefficients
  // transposed in the sums and diffs shared memory arrays.
  int readIdx = inputy * arrayWidth + inputx;
  // (inputx < (transformLength<<1) ??
  if (inputx+1 < transformLength && inputy < transformLength) {
    NUM a = data[readIdx], b = data[readIdx+1];
    sums [threadIdx.x][threadIdx.y] = (a + b) * INV_SQRT2;
    diffs[threadIdx.x][threadIdx.y] = (a - b) * INV_SQRT2;
  }

  __syncthreads();

  // Read the transposed sums and diffs shared memory arrays,
  // and write the data to a tile whose position has been transposed
  int writey = blockIdx.x*blockDim.x + threadIdx.y;
  int writex = blockIdx.y*blockDim.y + threadIdx.x;
  if (writex < transformLength && writey*2 < transformLength) {
    int writeIdx = writey * arrayWidth + writex;
    result[writeIdx] = sums[threadIdx.y][threadIdx.x];
    writeIdx += arrayWidth*(transformLength>>1);
    result[writeIdx] = diffs[threadIdx.y][threadIdx.x];
  }
}

template<typename NUM>
__global__ void haar_inv_transpose_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result) {
  __shared__ NUM sums [HT_TILE_SIZE][HT_TILE_SIZE];
  __shared__ NUM diffs[HT_TILE_SIZE][HT_TILE_SIZE];

  int inputx = blockIdx.x*blockDim.x + threadIdx.x;
  int inputy = blockIdx.y*blockDim.y + threadIdx.y;
  
  // read a tile 2*HT_TILE_SIZE wide, HT_TILE_SIZE tall, compute
  // the sum and difference coefficients, and store those coefficients
  // transposed in the sums and diffs shared memory arrays.
  int readIdx1 = inputy * arrayWidth + inputx;
  int readIdx2 = readIdx1 + (transformLength>>1);
  if (inputx < (transformLength>>1) && inputy < transformLength) {
    NUM s = data[readIdx1], d = data[readIdx2];
    sums [threadIdx.x][threadIdx.y] = (s + d) * INV_SQRT2;
    diffs[threadIdx.x][threadIdx.y] = (s - d) * INV_SQRT2;
  }

  __syncthreads();

  // Read the transposed sums and diffs shared memory arrays,
  // and write the data to a tile whose position has been transposed
  int writex = blockIdx.y*blockDim.y + threadIdx.x;
  int writey = (blockIdx.x*blockDim.x + threadIdx.y) * 2;
  if (writex < transformLength && writey+1 < transformLength) {
    int writeIdx1 = writey * arrayWidth + writex;
    int writeIdx2 = writeIdx1 + arrayWidth;

    result[writeIdx1] = sums[threadIdx.y][threadIdx.x];
    result[writeIdx2] = diffs[threadIdx.y][threadIdx.x];
  }
}


float haar_not_lifting_2d_cuda
  (int size, float *data, bool inverse, int stepCount, int threadBlockSize,
   bool useCombinedTranspose) {
  return haar_not_lifting_2d_cuda_internal
    (size, data, inverse, stepCount, threadBlockSize, useCombinedTranspose);
}     


float haar_not_lifting_2d_cuda
  (int size, double *data, bool inverse, int stepCount, int threadBlockSize,
   bool useCombinedTranspose) {
  return haar_not_lifting_2d_cuda_internal
    (size, data, inverse, stepCount, threadBlockSize, useCombinedTranspose);
}     


// Wrapper function that handles the CUDA details.
template<typename NUM>
float haar_not_lifting_2d_cuda_internal
(int size, NUM *data, bool inverse, int stepCount, int threadBlockSize,
 bool useCombinedTranspose) {

  // printf("Tile size %d\n", HT_TILE_SIZE);
  int maxSteps = dwtMaximumSteps(size);
  if (stepCount < 1 || stepCount > maxSteps)
    stepCount = maxSteps;

  // create timers
  CudaTimer overallTimer, copyToTimer, copyFromTimer, 
    transformTimer, transposeTimer;

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
  overallTimer.start(stream);

  // copy the data to the GPU
  copyToTimer.start(stream);
  CUCHECK(cudaMemcpyAsync(data1_dev, data, totalBytes, cudaMemcpyHostToDevice,
                          stream));
  copyToTimer.end(stream);

  int transformLength;

  if (inverse) {

    // inverse
    transformLength = size >> (stepCount - 1);
    for (int i=0; i < stepCount; i++) {

      dim3 gridDim((transformLength - 1) / (HT_TILE_SIZE*2) + 1,
                   (transformLength - 1) / (HT_TILE_SIZE) + 1);
      dim3 blockDim(HT_TILE_SIZE, HT_TILE_SIZE);

      if (useCombinedTranspose) {

        // transform columns and transpose
        transformTimer.start(stream);
        haar_inv_transpose_2d_kernel
          <<<gridDim, blockDim, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);
    
        // transform rows and transpose
        transformTimer.start(stream);
        haar_inv_transpose_2d_kernel
          <<<gridDim, blockDim, 0, stream>>>
          (size, transformLength, data2_dev, data1_dev);
        transformTimer.end(stream);

      } else {

        // transform columns
        transformTimer.start(stream);
        haar_inv_not_lifting_2d_kernel
          <<<transformLength, threadBlockSize, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);

        // transpose the matrix into temp_dev
        transposeTimer.start(stream);
        gpuTranspose(size, transformLength, data2_dev, data1_dev, stream);
        transposeTimer.end(stream);
    
        // transform rows
        transformTimer.start(stream);
        haar_inv_not_lifting_2d_kernel
          <<<transformLength, threadBlockSize, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);

        // transpose the matrix into data_dev
        transposeTimer.start(stream);
        gpuTranspose(size, transformLength, data2_dev, data1_dev, stream);
        transposeTimer.end(stream);

        // results are in data1_dev

      }

      transformLength <<= 1;
    }

  } else {

    // forward
    transformLength = size;
    
    for (int i=0; i < stepCount; i++) {

      dim3 gridDim((transformLength - 1) / (HT_TILE_SIZE*2) + 1,
                   (transformLength - 1) / (HT_TILE_SIZE) + 1);
      dim3 blockDim(HT_TILE_SIZE, HT_TILE_SIZE);
    
      if (useCombinedTranspose) {

        // do the wavelet transform on rows
        transformTimer.start(stream);
        haar_transpose_2d_kernel
          <<<gridDim, blockDim, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);

        // do the wavelet transform on columns
        transformTimer.start(stream);
        haar_transpose_2d_kernel
          <<<gridDim, blockDim, 0, stream>>>
          (size, transformLength, data2_dev, data1_dev);
        transformTimer.end(stream);

      } else {

        // do the wavelet transform on rows
        transformTimer.start(stream);
        haar_not_lifting_2d_kernel
          <<<transformLength, threadBlockSize, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);

        // transpose the matrix into temp_dev
        transposeTimer.start(stream);
        gpuTranspose(size, transformLength, data2_dev, data1_dev, stream);
        transposeTimer.end(stream);
    
        // do the wavelet transform on columns
        transformTimer.start(stream);
        haar_not_lifting_2d_kernel
          <<<transformLength, threadBlockSize, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);
    
        // transpose the matrix back into data_dev
        transposeTimer.start(stream);
        gpuTranspose(size, transformLength, data2_dev, data1_dev, stream);
        transposeTimer.end(stream);

      }

      transformLength >>= 1;
    }

  }

  // copy the data back from the GPU
  copyFromTimer.start(stream);
  CUCHECK(cudaMemcpyAsync(data, data1_dev, totalBytes, cudaMemcpyDeviceToHost,
                          stream));
  copyFromTimer.end(stream);


  // Since all the GPU tasks were started asynchronously, control should
  // flow to this point very quickly. The cudaEventSynchronize() call will
  // wait until the GPU is finished.
  double endTimeCPU = NixTimer::time();
  printf("Time elapsed creating GPU tasks: %.3f ms\n",
         1000*(endTimeCPU - startTimeCPU));
  fflush(stdout);

  // stop the timer
  overallTimer.end(stream);
  CUCHECK(cudaEventSynchronize(overallTimer.getLastEvent()));
  cudaStreamDestroy(stream);

  // check for errors
  CUCHECK(cudaGetLastError());

  printf("Times:\n");
  printf("  Copy data to GPU:   %9.3f ms\n", copyToTimer.time());
  printf("  Transform time:     %9.3f ms (%d calls)\n",
         transformTimer.time(), transformTimer.count());
  if (transposeTimer.count() > 0) {
    printf("  Transpose time:     %9.3f ms (%d calls)\n", 
           transposeTimer.time(), transposeTimer.count());
  }
  printf("  Copy data from GPU: %9.3f ms\n", copyFromTimer.time());

  // deallocate GPU memory
  CUCHECK(cudaFree(data1_dev));
  CUCHECK(cudaFree(data2_dev));

  return overallTimer.time();
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
  printf("Time elapsed creating GPU tasks: %.3f ms\n",
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
