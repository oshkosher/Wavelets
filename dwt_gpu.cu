/*
  CUDA implementation of Haar discrete wavelet transform.

  Ed Karrels, ed.karrels@gmail.com, June 2014
*/

#include "dwt_cpu.h"
#include "dwt_gpu.h"
#include "transpose_gpu.h"
#include "cucheck.h"
#include "nixtimer.h"
#include "cuda_timer.h"


#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f

#define HAAR_3D_BLOCK_SIZE 128

/*
  To see a previous version of this code that tried out surfaces and
  iterating down columns (rather than rows), see the version as checked
  into Git history as checkin "736c52c":
    git show 736c52c:dwt_gpu.cu | less
*/


/*
  Call structure:

  haar_2d_cuda
    haar_2d_cuda_internal
      haar_2d_kernel
      haar_inv_2d_kernel
      haar_transpose_2d_kernel
      haar_inv_transpose_2d_kernel

*/


template<typename NUM>
float haar_2d_cuda_internal
(int size, NUM *data, bool inverse, int stepCount, int threadBlockSize,
 bool useCombinedTranspose);


/**
   Perform one pass of the Haar transform on the first 'length' 
   elements of inputRow using outputRow as temp space.
*/
template<typename NUM>
__device__ void haar_kernel_row(int length, NUM *inputRow, NUM *outputRow) {

  int half = length >> 1;

  // point d at the first half of the temporary row
  NUM *s = outputRow;
  
  // point d at the second half of the temporary row
  NUM *d = s + half;
  
  for (int i=threadIdx.x; i < half; i += blockDim.x) {
    NUM a = inputRow[2*i], b = inputRow[2*i + 1];
    d[i] = (a - b) * INV_SQRT2;
    s[i] = (a + b) * INV_SQRT2;
  }
}


/**
   Perform one pass of the inverse Haar transform on the first 'length' 
   elements of inputRow using outputRow as temp space.
*/
template<typename NUM>
__device__ void haar_kernel_row_inverse(int length, NUM *inputRow,
                                        NUM *outputRow) {

  // Set s to point to my row in the input data
  NUM *s = inputRow;

  int half = length >> 1;

  // point d at the second half of the temporary row
  NUM *d = s + half;

  for (int i=threadIdx.x; i < half; i += blockDim.x) {
    outputRow[2*i]   = INV_SQRT2 * (s[i] + d[i]);
    outputRow[2*i+1] = INV_SQRT2 * (s[i] - d[i]);
  }
}



/*
  This does a Haar discrete wavelet transform on each row of
  a 2-d array. Each thread block processes one row.
  All data is in global memory.
  Input data is in data[], results will be in result[].
*/
template<typename NUM>
__global__ void haar_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result) {

  // each thread block processes one row of data
  int y = blockIdx.x;

  // make pointers to my row of data
  NUM *inputRow = data + y * arrayWidth;
  NUM *outputRow = result + y * arrayWidth;

  haar_kernel_row(transformLength, inputRow, outputRow);
}


/* Inverse Haar wavelet transform. */
template<typename NUM>
__global__ void haar_inv_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result) {

  // each thread block processes one row of data
  int y = blockIdx.x;

  // make pointers to my row of data
  NUM *inputRow = data + y * arrayWidth;
  NUM *outputRow = result + y * arrayWidth;

  haar_kernel_row_inverse(transformLength, inputRow, outputRow);
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
(int arrayWidth, int transformLength, NUM *data, NUM *result,
 int tileSize) {

  // dynamically-sized shared memory
  extern __shared__ int shared[];

  // assign parts of shared memory to my arrays
  NUM *sums, *diffs;
  sums = (NUM*) shared;
  diffs = sums + tileSize * (tileSize+1);

  int inputx = (blockIdx.x*blockDim.x + threadIdx.x) * 2;
  int inputy = blockIdx.y*blockDim.y + threadIdx.y;
  
  // read a tile 2*tileSize wide, tileSize tall, compute
  // the sum and difference coefficients, and store those coefficients
  // transposed in the sums and diffs shared memory arrays.
  int readIdx = inputy * arrayWidth + inputx;

  if (inputx+1 < transformLength && inputy < transformLength) {
    NUM a = data[readIdx], b = data[readIdx+1];
    int shidx = threadIdx.x + threadIdx.y*(tileSize+1);
    sums [shidx] = (a + b) * INV_SQRT2;
    diffs[shidx] = (a - b) * INV_SQRT2;
  }

  __syncthreads();

  // Read the transposed sums and diffs shared memory arrays,
  // and write the data to a tile whose position has been transposed
  int writey = blockIdx.x*blockDim.x + threadIdx.y;
  int writex = blockIdx.y*blockDim.y + threadIdx.x;
  if (writex < transformLength && writey*2 < transformLength) {
    int writeIdx = writey * arrayWidth + writex;
    int shidx = threadIdx.y + threadIdx.x*(tileSize+1);
    result[writeIdx] = sums[shidx];
    writeIdx += arrayWidth*(transformLength>>1);
    result[writeIdx] = diffs[shidx];
  }
}

template<typename NUM>
__global__ void haar_inv_transpose_2d_kernel
(int arrayWidth, int transformLength, NUM *data, NUM *result, int tileSize) {

  // dynamically-sized shared memory
  extern __shared__ int shared[];

  // assign parts of shared memory to my arrays
  NUM *v1, *v2;
  v1 = (NUM*) shared;
  v2 = v1 + tileSize * (tileSize+1);

  int inputx = blockIdx.x*blockDim.x + threadIdx.x;
  int inputy = blockIdx.y*blockDim.y + threadIdx.y;

  // Read the sum and difference coefficients, where the difference coeff
  // is in the second half of the array. Compute the original values v1 and v2,
  // and store them in two shared memory arrays.
  int readIdx1 = inputy * arrayWidth + inputx;
  int readIdx2 = readIdx1 + (transformLength>>1);
  if (inputx < (transformLength>>1) && inputy < transformLength) {
    NUM s = data[readIdx1], d = data[readIdx2];
    int shidx = threadIdx.x * (tileSize+1) + threadIdx.y;
    v1[shidx] = (s + d) * INV_SQRT2;
    v2[shidx] = (s - d) * INV_SQRT2;
  }

  __syncthreads();

  // Read the transposed pair of values v1 and v2 from the transposed
  // shared memory arrays, and write the values to a tile tileSize wide
  // and tileSize*2 tall.
  int writex = blockIdx.y*blockDim.y + threadIdx.x;
  int writey = (blockIdx.x*blockDim.x + threadIdx.y) * 2;
  if (writex < transformLength && writey+1 < transformLength) {
    int writeIdx1 = writey * arrayWidth + writex;
    int writeIdx2 = writeIdx1 + arrayWidth;
    int shidx = threadIdx.y * (tileSize+1) + threadIdx.x;
    result[writeIdx1] = v1[shidx];
    result[writeIdx2] = v2[shidx];
  }
}


float haar_2d_cuda
  (int size, float *data, bool inverse, int stepCount, int threadBlockSize,
   bool useCombinedTranspose) {
  return haar_2d_cuda_internal
    (size, data, inverse, stepCount, threadBlockSize, useCombinedTranspose);
}     


// double support was added in version 1.3
#if !defined(__CUDA_ARCH__) || (__CUDA_ARCH__ >= 130)
float haar_2d_cuda
  (int size, double *data, bool inverse, int stepCount, int threadBlockSize,
   bool useCombinedTranspose) {
  return haar_2d_cuda_internal
    (size, data, inverse, stepCount, threadBlockSize, useCombinedTranspose);
}
#endif


// haar_transpose_2d_kernel and haar_inv_transpose_2d_kernel use tiles
// to optimize the memory access pattern. After testing tile sizes from
// 8x8 to 32x32 on a few different GPUs, here are the sizes that produced
// the best performance:
//
//   GTX 480: 16     (compute level 2.0)
//   GTX 570: 16     (compute level 2.0)
//   K2000M: 16      (compute level 3.0, laptop)
//   GTX 680: 32     (compute level 3.0)
//   GTX 690: 32     (compute level 3.0)
//   Tesla K20c: 32  (compute level 3.5)
//
int bestHaarGPUTileSize() {
  int gpuId;
  cudaDeviceProp prop;
  CUCHECK(cudaGetDevice(&gpuId));
  CUCHECK(cudaGetDeviceProperties(&prop, gpuId));

  // Based on the tests listed above, older (Fermi) and smaller (laptop)
  // GPUs seem to work better with 16x16 tiles, but newer regular GPUs
  // are faster with 32x32 tiles.
  if (prop.major <= 2 || prop.multiProcessorCount <= 2)
    return 16;
  else
    return 32;
}


// Wrapper function that handles the CUDA details.
template<typename NUM>
float haar_2d_cuda_internal
(int size, NUM *data, bool inverse, int stepCount, int threadBlockSize,
 bool useCombinedTranspose) {

  int tileSize = bestHaarGPUTileSize();

  if (useCombinedTranspose) printf("Tile size %dx%d\n", tileSize, tileSize);

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

  size_t sharedMemSize = tileSize * (tileSize+1)
    * 2 * sizeof(float);
  
  int transformLength;

  if (inverse) {

    // inverse
    transformLength = size >> (stepCount - 1);
    for (int i=0; i < stepCount; i++) {

      dim3 gridDim((transformLength - 1) / (tileSize*2) + 1,
                   (transformLength - 1) / (tileSize) + 1);
      dim3 blockDim(tileSize, tileSize);

      if (useCombinedTranspose) {

        // transform columns and transpose
        transformTimer.start(stream);
        haar_inv_transpose_2d_kernel
          <<<gridDim, blockDim, sharedMemSize, stream>>>
          (size, transformLength, data1_dev, data2_dev, tileSize);
        transformTimer.end(stream);
    
        // transform rows and transpose
        transformTimer.start(stream);
        haar_inv_transpose_2d_kernel
          <<<gridDim, blockDim, sharedMemSize, stream>>>
          (size, transformLength, data2_dev, data1_dev, tileSize);
        transformTimer.end(stream);

      } else {

        // transform columns
        transformTimer.start(stream);
        haar_inv_2d_kernel
          <<<transformLength, threadBlockSize, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);

        // transpose the matrix into temp_dev
        transposeTimer.start(stream);
        gpuTransposeSquare(size, transformLength, data2_dev, data1_dev, stream);
        transposeTimer.end(stream);
    
        // transform rows
        transformTimer.start(stream);
        haar_inv_2d_kernel
          <<<transformLength, threadBlockSize, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);

        // transpose the matrix into data_dev
        transposeTimer.start(stream);
        gpuTransposeSquare(size, transformLength, data2_dev, data1_dev, stream);
        transposeTimer.end(stream);

        // results are in data1_dev

      }

      transformLength <<= 1;
    }

  } else {

    // forward
    transformLength = size;
    
    for (int i=0; i < stepCount; i++) {

      dim3 gridDim((transformLength - 1) / (tileSize*2) + 1,
                   (transformLength - 1) / (tileSize) + 1);
      dim3 blockDim(tileSize, tileSize);
    
      if (useCombinedTranspose) {

        // do the wavelet transform on rows
        transformTimer.start(stream);
        haar_transpose_2d_kernel
          <<<gridDim, blockDim, sharedMemSize, stream>>>
          (size, transformLength, data1_dev, data2_dev, tileSize);
        transformTimer.end(stream);

        // do the wavelet transform on columns
        transformTimer.start(stream);
        haar_transpose_2d_kernel
          <<<gridDim, blockDim, sharedMemSize, stream>>>
          (size, transformLength, data2_dev, data1_dev, tileSize);
        transformTimer.end(stream);

      } else {

        // do the wavelet transform on rows
        transformTimer.start(stream);
        haar_2d_kernel
          <<<transformLength, threadBlockSize, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);

        // transpose the matrix into temp_dev
        transposeTimer.start(stream);
        gpuTransposeSquare(size, transformLength, data2_dev, data1_dev, stream);
        transposeTimer.end(stream);
    
        // do the wavelet transform on columns
        transformTimer.start(stream);
        haar_2d_kernel
          <<<transformLength, threadBlockSize, 0, stream>>>
          (size, transformLength, data1_dev, data2_dev);
        transformTimer.end(stream);
    
        // transpose the matrix back into data_dev
        transposeTimer.start(stream);
        gpuTransposeSquare(size, transformLength, data2_dev, data1_dev, stream);
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


__device__ void copy_row(float *dest, float *src, int count) {
  int i = threadIdx.x;
  while (i < count) {
    dest[i] = src[i];
    i += blockDim.x;
  }
}

__global__ void haar_3d_kernel_allglobal
(float *data1, float *data2, int rowLength, int stepCount) {
  
  int offset = rowLength * (blockIdx.y + blockIdx.z*gridDim.y);
  float *row = data1 + offset;
  float *tempRow = data2 + offset;
  
  while (stepCount > 0) {

    // copy data from tempRow to row
    copy_row(tempRow, row, rowLength);

    __syncthreads();

    // read tempRow, write to row
    haar_kernel_row(rowLength, tempRow, row);

    stepCount--;
    rowLength >>= 1;

    __syncthreads();
  }
}


__global__ void haar_3d_kernel_allglobal_inverse
(float *data1, float *data2, int rowLength, int stepCount) {
  
  int offset = rowLength * (blockIdx.y + blockIdx.z*gridDim.y);
  float *row = data1 + offset;
  float *tempRow = data2 + offset;

  rowLength >>= (stepCount-1);
  
  while (stepCount > 0) {

    // copy data from tempRow to row
    copy_row(tempRow, row, rowLength);

    __syncthreads();

    // read tempRow, write to row
    haar_kernel_row_inverse(rowLength, tempRow, row);

    stepCount--;
    rowLength <<= 1;

    __syncthreads();
  }
}

/**
  Input data is in data1. data2 is temp space.
  Update 'size', since the dimensions will rotate.

  On each pass, the input row is copied to a temp row.
  Entries in the temp row are read, processed, and written to the input row.
  
  If the row is short enough to fit in shared memory, use shared memory
  as the temp row.

  If the row is short enough to fit two copies in shared memory
    1. Copy input row from global to shared memory.
    2. Do regular processing with shared memory as temp row.
    3. Copy input row from shared to global memory.

  If the row is short enough to fit in registers, use them as the temp row.
*/
void haar_3d_cuda(float *data, float *tmpData, scu_wavelet::int3 &size,
                  scu_wavelet::int3 stepCount, bool inverse,
                  bool isStandardTranspose,
                  CudaTimer *transformTimer,
                  CudaTimer *transposeTimer) {

  if (!(size <= scu_wavelet::int3(65535,65535,65535))) {
    fprintf(stderr, "Cubelet too large: max is 65535x65535x65535\n");
    return;
  }

  if (!isStandardTranspose) {
    fprintf(stderr, "Nonstandard transpose not implemented.\n");
    return;
  }

  if (inverse) stepCount.rotateBack();
  if (!is_padded_for_wavelet(size, stepCount)) {
    fprintf(stderr, "%dx%dx%d data is not properly padded for %d,%d,%d "
            "transform steps.\n",
            size.x, size.y, size.z, stepCount.x, stepCount.y, stepCount.z);
    return;
  }
  if (inverse) stepCount.rotateFwd();

  dim3 gridDim;
  dim3 blockDim(HAAR_3D_BLOCK_SIZE);

  if (!inverse) {

    // X transform
    gridDim.y = size.y;
    gridDim.z = size.z;
    if (transformTimer) transformTimer->start();
    haar_3d_kernel_allglobal<<<gridDim, blockDim>>>
      (data, tmpData, size.x, stepCount.x);
    if (transformTimer) transformTimer->end();

    // transpose XYZ -> YZX
    if (transposeTimer) transposeTimer->start();
    gpuTranspose3dFwd(tmpData, data, size);
    if (transposeTimer) transposeTimer->end();

    // data is now in tmpData

    // Y transform
    gridDim.y = size.y;
    gridDim.z = size.z;
    if (transformTimer) transformTimer->start();
    haar_3d_kernel_allglobal<<<gridDim, blockDim>>>
      (tmpData, data, size.x, stepCount.y);
    if (transformTimer) transformTimer->end();

    // transpose YZX -> ZXY
    if (transposeTimer) transposeTimer->start();
    gpuTranspose3dFwd(data, tmpData, size);
    if (transposeTimer) transposeTimer->end();

    // data is back in data

    // Z transform
    gridDim.y = size.y;
    gridDim.z = size.z;
    if (transformTimer) transformTimer->start();
    haar_3d_kernel_allglobal<<<gridDim, blockDim>>>
      (data, tmpData, size.x, stepCount.z);
    if (transformTimer) transformTimer->end();

  } else { // is inverse

    // inverse Z transform
    gridDim.y = size.y;
    gridDim.z = size.z;
    if (transformTimer) transformTimer->start();
    haar_3d_kernel_allglobal_inverse<<<gridDim, blockDim>>>
      (data, tmpData, size.x, stepCount.z);
    if (transformTimer) transformTimer->end();

    // transpose ZXY -> YZX
    if (transposeTimer) transposeTimer->start();
    gpuTranspose3dBack(tmpData, data, size);
    if (transposeTimer) transposeTimer->end();

    // data is in 'tmpData'

    // inverse Y transform
    gridDim.y = size.y;
    gridDim.z = size.z;
    if (transformTimer) transformTimer->start();
    haar_3d_kernel_allglobal_inverse<<<gridDim, blockDim>>>
      (tmpData, data, size.x, stepCount.y);
    if (transformTimer) transformTimer->end();

    // transpose YZX -> XYZ
    if (transposeTimer) transposeTimer->start();
    gpuTranspose3dBack(data, tmpData, size);
    if (transposeTimer) transposeTimer->end();

    // data is in 'data'

    // inverse X transform
    gridDim.y = size.y;
    gridDim.z = size.z;
    if (transformTimer) transformTimer->start();
    haar_3d_kernel_allglobal_inverse<<<gridDim, blockDim>>>
      (data, tmpData, size.x, stepCount.x);
    if (transformTimer) transformTimer->end();
  }

    /*
      CUCHECK(cudaMemcpy(data, tmpData, sizeof(float)*size.x*size.y*size.z,
      cudaMemcpyDeviceToDevice));

Mean squared error: 469.826, peak SNR: 21.411
Huffman build table 0.779 ms
Huffman encoding: 288 bytes, 4.48 bits/pixel, longest encoding = 10 bits
Write data file: 3.32 ms
Total: 19.02 ms

    */
  
}
