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


#include "test_compress_gpu.h"

#define SQRT2     1.4142135623730950488f
#define INV_SQRT2 0.70710678118654752440f

#define HAAR_3D_BLOCK_SIZE 128
#define CDF97_3D_BLOCK_SIZE 128

#define HAAR_V2_BLOCK_SIZE 128
#define CDF97_V2_BLOCK_SIZE 128

#define CDF97_V3_BLOCK_WIDTH 32
#define CDF97_V3_BLOCK_HEIGHT 10
#define CDF97_V3_TILE_WIDTH 32
#define CDF97_V3_TILE_HEIGHT 32

#define MAX_KERNEL_ROW_LEN 12000
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


void printArray(const float *array, int width, int height, int depth, 
                const char *name) {
  if (name) printf("%s\n", name);
  for (int level=0; level < depth; level++) {
    printf("z=%d\n", level);
    for (int row=0; row < height; row++) {
      for (int col=0; col < width; col++) {
        printf("%8.4f ", array[(level*height + row)*width + col]);
      }
      putchar('\n');
    }
    putchar('\n');
  }
  putchar('\n');
}


void printArray(const int *array, int width, int height, int depth, 
                const char *name) {
  if (name) printf("%s\n", name);
  for (int level=0; level < depth; level++) {
    printf("z=%d\n", level);
    for (int row=0; row < height; row++) {
      for (int col=0; col < width; col++) {
        printf("%5d ", array[(level*height + row)*width + col]);
      }
      putchar('\n');
    }
    putchar('\n');
  }
  putchar('\n');
}


void printDeviceArray(const float *array_dev, scu_wavelet::int3 size,
                      const char *name) {
  printDeviceArray(array_dev, size.x, size.y, size.z, name);
}

void printDeviceArray(const float *array_dev, int width, int height, int depth, 
                      const char *name) {

  float *array = new float[width*height*depth];
  CUCHECK(cudaMemcpy(array, array_dev, sizeof(float)*width*height*depth,
                     cudaMemcpyDeviceToHost));

  printArray(array, width, height, depth, name);

  delete[] array;
}

void printDeviceArray(const int *array_dev, scu_wavelet::int3 size,
                      const char *name) {
  printDeviceArray(array_dev, size.x, size.y, size.z, name);
}

void printDeviceArray(const int *array_dev, int width, int height, int depth, 
                      const char *name) {

  int *array = new int[width*height*depth];
  CUCHECK(cudaMemcpy(array, array_dev, sizeof(int)*width*height*depth,
                     cudaMemcpyDeviceToHost));

  printArray(array, width, height, depth, name);

  delete[] array;
}




/**
   Perform one pass of the Haar transform on the first 'length'
   elements of inputRow using outputRow as temp space.
*/
template<typename NUM>
__device__ void haar_kernel_row(int length, NUM *outputRow,
                                const NUM *inputRow) {

  const int half = length >> 1;

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
__device__ void haar_kernel_row_inverse(int length, NUM *outputRow,
                                        NUM *inputRow) {

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

  haar_kernel_row(transformLength, outputRow, inputRow);
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

  haar_kernel_row_inverse(transformLength, outputRow, inputRow);
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


__device__ void copy_row(int count, float *dest, const float *src) {
  int i = threadIdx.x;
  while (i < count) {
    dest[i] = src[i];
    i += blockDim.x;
  }
}


__device__ void pad_row_mirrored(int length, float *row, int padLen) {
  assert(blockDim.x >= padLen*2);
  if (threadIdx.x < padLen*2) {
    int factor, offset;
    if (threadIdx.x < padLen) {  // first few threads pad the beginning
      factor = -1;
      offset = -1;
    } else {  // other half of the threads pad the end
      factor = 1;
      offset = length-padLen;
    }

    MirroredArray<float> mirrored(length, row);
    int i = threadIdx.x * factor + offset;
    row[i] = mirrored[i];
  }
}

__global__ void haar_3d_kernel_allglobal
(float *data1, float *data2, int rowLength, int stepCount) {

  int offset = rowLength * (blockIdx.y + blockIdx.z*gridDim.y);
  float *row = data1 + offset;
  float *tempRow = data2 + offset;

  while (stepCount > 0) {

    // copy data from row to tempRow
    copy_row(rowLength, tempRow, row);

    __syncthreads();

    // read tempRow, write to row
    haar_kernel_row(rowLength, row, tempRow);

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

    // copy data from row to tempRow
    copy_row(rowLength, tempRow, row);

    __syncthreads();

    // read tempRow, write to row
    haar_kernel_row_inverse(rowLength, row, tempRow);

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
                  CudaTimer *transformTimer,
                  CudaTimer *transposeTimer) {

  // limitation on the number of thread blocks
  if (!(size <= scu_wavelet::int3(65535,65535,65535))) {
    fprintf(stderr, "Cubelet too large: max is 65535x65535x65535\n");
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

  CudaTimer mytimer("haar_3d_cuda");
  mytimer.start();

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

  mytimer.end();
  mytimer.sync();
  mytimer.print();

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


template<class ARRAY>
__device__ void cdf97_row(int rowLength, float *outputRow,
                          const ARRAY  inputRow) {

  const int half = rowLength >> 1;

  int writeIdx = threadIdx.x;

  while (writeIdx < half) {
    int readIdx = writeIdx << 1;

    // Apply the sums convolution, write result to lower half of the row

    float t0, t1, t2, t3, t4, t5, t6, t7, t8;
    t0 = inputRow[readIdx - 4];
    t1 = inputRow[readIdx - 3];
    t2 = inputRow[readIdx - 2];
    t3 = inputRow[readIdx - 1];
    t4 = inputRow[readIdx];
    t5 = inputRow[readIdx + 1];
    t6 = inputRow[readIdx + 2];
    t7 = inputRow[readIdx + 3];
    t8 = inputRow[readIdx + 4];

    outputRow[writeIdx] =
      CDF97_ANALYSIS_LOWPASS_FILTER_0 * t4
      + CDF97_ANALYSIS_LOWPASS_FILTER_1 * (t3+t5)
      + CDF97_ANALYSIS_LOWPASS_FILTER_2 * (t2+t6)
      + CDF97_ANALYSIS_LOWPASS_FILTER_3 * (t1+t7)
      + CDF97_ANALYSIS_LOWPASS_FILTER_4 * (t0+t8);

    // Apply the differences convolution, write result to upper half of the row

    outputRow[writeIdx+half] =
      CDF97_ANALYSIS_HIGHPASS_FILTER_0 * t5
      + CDF97_ANALYSIS_HIGHPASS_FILTER_1 * (t4+t6)
      + CDF97_ANALYSIS_HIGHPASS_FILTER_2 * (t3+t7)
      + CDF97_ANALYSIS_HIGHPASS_FILTER_3 * (t2+t8);

    writeIdx += blockDim.x;
  }

}


// interleave: 01234567 -> 04152636
__device__ void cdf97_row_interleave(int rowLength, float *outputRow,
                                            const float *inputRow) {
  const int half = rowLength >> 1;
  int i = threadIdx.x;

  while (i < half) {
    outputRow[i*2] = inputRow[i];
    outputRow[i*2+1] = inputRow[i+half];
    i += blockDim.x;
  }
}


template<class ARRAY>
__device__ void cdf97_row_inverse(int rowLength, float *outputRow,
                                  const ARRAY inputRow) {

  int i = threadIdx.x * 2;

  while (i < rowLength) {
    float t0, t1, t2, t3, t4, t5, t6, t7, t8;

    t0 = inputRow[i-3];
    t1 = inputRow[i-2];
    t2 = inputRow[i-1];
    t3 = inputRow[i];
    t4 = inputRow[i+1];
    t5 = inputRow[i+2];
    t6 = inputRow[i+3];
    t7 = inputRow[i+4];
    t8 = inputRow[i+5];

    // evens
    outputRow[i] =
      CDF97_SYNTHESIS_LOWPASS_FILTER_0 * t3
      + CDF97_SYNTHESIS_LOWPASS_FILTER_1 * (t2+t4)
      + CDF97_SYNTHESIS_LOWPASS_FILTER_2 * (t1+t5)
      + CDF97_SYNTHESIS_LOWPASS_FILTER_3 * (t0+t6);

    // odds
    outputRow[i+1] =
      CDF97_SYNTHESIS_HIGHPASS_FILTER_0 * t4
      + CDF97_SYNTHESIS_HIGHPASS_FILTER_1 * (t3+t5)
      + CDF97_SYNTHESIS_HIGHPASS_FILTER_2 * (t2+t6)
      + CDF97_SYNTHESIS_HIGHPASS_FILTER_3 * (t1+t7)
      + CDF97_SYNTHESIS_HIGHPASS_FILTER_4 * (t0+t8);

    i += blockDim.x*2;
  }

}


__global__ void cdf97_3d_kernel
(float *data, int rowLength, int stepCount) {

  int offset = rowLength * (blockIdx.y + blockIdx.z*gridDim.y);
  float *row = data+offset;

  extern __shared__ float sharedData[];
  float *tempRow = sharedData+4;

  while (stepCount > 0) {

    // copy data from global memory to tempRow
    copy_row(rowLength, tempRow, row);
    __syncthreads();

    // pad the data
    pad_row_mirrored(rowLength, tempRow, 4);
    __syncthreads();

    // read tempRow, write to row in global memory
    cdf97_row(rowLength, row, tempRow);

    stepCount--;
    rowLength >>= 1;

    __syncthreads();
  }

}


__global__ void cdf97_3d_kernel_inverse
(float *data, float *tempData, int rowLength, int stepCount) {

  int offset = rowLength * (blockIdx.y + blockIdx.z*gridDim.y);
  float *row = data + offset;
  // float *tempRow = tempData + offset;

  extern __shared__ float sharedData[];
  float *tempRow = sharedData+4;

  rowLength >>= (stepCount-1);

  // MirroredArray inputRow(rowLength, tempRow);

  while (stepCount > 0) {
    // inputRow.setLength(rowLength);

    // copy data to tempData, interleaving 01234567 -> 04152636
    cdf97_row_interleave(rowLength, tempRow, row);
    __syncthreads();

    pad_row_mirrored(rowLength, tempRow, 4);
    __syncthreads();

    // read tempRow, write to row
    cdf97_row_inverse(rowLength, row, tempRow);

    stepCount--;
    rowLength <<= 1;

    __syncthreads();
  }
}


void cdf97_3d_cuda(float *data, float *tmpData,
                   scu_wavelet::int3 &size, scu_wavelet::int3 stepCount,
                   bool inverse,
                   CudaTimer *transformTimer, CudaTimer *transposeTimer) {

  // limitation based on the size of shared memory
  if (!(size <= scu_wavelet::int3(MAX_KERNEL_ROW_LEN,MAX_KERNEL_ROW_LEN,
                                  MAX_KERNEL_ROW_LEN))) {
    fprintf(stderr, "Cubelet too large: max is %d pixels on any side.\n",
            MAX_KERNEL_ROW_LEN);
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
  dim3 blockDim(CDF97_3D_BLOCK_SIZE);
  size_t sharedMemSize;

  if (!inverse) {

    // X transform
    gridDim.y = size.y;
    gridDim.z = size.z;
    if (transformTimer) transformTimer->start();

    // add 8 for 4 pad entries on either side of the second copy
    sharedMemSize = sizeof(float) * (size.x + 8);
    cdf97_3d_kernel<<<gridDim, blockDim, sharedMemSize>>>
      (data, size.x, stepCount.x);

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

    sharedMemSize = sizeof(float) * (size.x + 8);
    cdf97_3d_kernel<<<gridDim, blockDim, sharedMemSize>>>
      (tmpData, size.x, stepCount.y);

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

    sharedMemSize = sizeof(float) * (size.x + 8);
    cdf97_3d_kernel<<<gridDim, blockDim, sharedMemSize>>>
      (data, size.x, stepCount.z);

    if (transformTimer) transformTimer->end();

  } else { // is inverse

    // inverse Z transform
    gridDim.y = size.y;
    gridDim.z = size.z;
    if (transformTimer) transformTimer->start();
    sharedMemSize = sizeof(float) * (size.x + 8);
    cdf97_3d_kernel_inverse<<<gridDim, blockDim, sharedMemSize>>>
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
    sharedMemSize = sizeof(float) * (size.x + 8);
    cdf97_3d_kernel_inverse<<<gridDim, blockDim, sharedMemSize>>>
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
    sharedMemSize = sizeof(float) * (size.x + 8);
    cdf97_3d_kernel_inverse<<<gridDim, blockDim, sharedMemSize>>>
      (data, tmpData, size.x, stepCount.x);
    if (transformTimer) transformTimer->end();

  }
}

/**
   Compute Haar wavelet transform in the Y direction

   Each thread block reads two adjacent rows from data_in[] and writes
   the results to data_out[].
*/
__global__ void haar_v2_kernel(int size_x, int size_y,
                               float *data_out, const float *data_in) {

  const float *row_in_1 = data_in +
    size_x * (blockIdx.y*2 + size_y * blockIdx.z);
  const float *row_in_2 = row_in_1 + size_x;
  float *sum_row_out = data_out +
    size_x * (blockIdx.y + size_y * blockIdx.z);
  float *diff_row_out = sum_row_out + (size_x * gridDim.y);

  for (int x = threadIdx.x; x < size_x; x += blockDim.x) {
    float a = row_in_1[x];
    float b = row_in_2[x];
    sum_row_out[x]  = (a + b) * INV_SQRT2;
    diff_row_out[x] = (a - b) * INV_SQRT2;
  }
}


// Copy the top blockDim.y rows from each plane from data_in to data_out.
__global__ void haar_v2_copyrows(int size_x, int size_y,
                                 float *data_out, const float *data_in) {

  size_t offset = size_x * (blockIdx.y + size_y * blockIdx.z);
  const float *src_row = data_in + offset;
  float *dest_row = data_out + offset;

  copy_row(size_x, dest_row, src_row);
}



// do 'level_count' transforms in the y direction
void haar_v2(float *data, float *data_tmp,
             scu_wavelet::int3 size, int level_count,
             CudaTimer *transformTimer) {

  dim3 gridDim, blockDim(HAAR_V2_BLOCK_SIZE);
  gridDim.x = 1;
  gridDim.y = size.y / 2;
  gridDim.z = size.z;

  if (transformTimer) transformTimer->start();

  for (int i=0; i < level_count; i++) {
    haar_v2_kernel<<<gridDim, blockDim>>>
      (size.x, size.y, data_tmp, data);

    // copy the results back into 'data'
    gridDim.y *= 2;
    haar_v2_copyrows<<<gridDim, blockDim>>>(size.x, size.y, data, data_tmp);

    gridDim.y /= 4;
  }

  if (transformTimer) transformTimer->end();
}


void haar_3d_cuda_v2(float *data, float *data_tmp, scu_wavelet::int3 &size,
                     scu_wavelet::int3 stepCount, bool inverse,
                     CudaTimer *transformTimer, CudaTimer *transposeTimer) {

  if (inverse) {
    haar_3d_cuda(data, data_tmp, size, stepCount, inverse,
                 transformTimer, transposeTimer);
    return;
  }

  CudaTimer mytimer("haar_3d_cuda_v2");
  mytimer.start();

  // change XYZ -> ZXY
  if (transposeTimer) transposeTimer->start();
  gpuTranspose3dBack(data_tmp, data, size);
  if (transposeTimer) transposeTimer->end();

  haar_v2(data_tmp, data, size, stepCount.x, transformTimer);

  // ZXY -> XYZ
  if (transposeTimer) transposeTimer->start();
  gpuTranspose3dFwd(data, data_tmp, size);
  if (transposeTimer) transposeTimer->end();

  haar_v2(data, data_tmp, size, stepCount.y, transformTimer);

  // XYZ -> YZX
  if (transposeTimer) transposeTimer->start();
  gpuTranspose3dFwd(data_tmp, data, size);
  if (transposeTimer) transposeTimer->end();

  haar_v2(data_tmp, data, size, stepCount.z, transformTimer);

  // YZX -> ZXY
  if (transposeTimer) transposeTimer->start();
  gpuTranspose3dFwd(data, data_tmp, size);
  if (transposeTimer) transposeTimer->end();

  mytimer.end();
  mytimer.sync();
  mytimer.print();
}


/**
   Compute CDF97 wavelet transform in the Y direction by reading an
   entire row of data at a time. Each thread block reads all the rows
   it needs for input and writes one row as output.
*/
__global__ void cdf97_v2_kernel(int size_x, int size_y,
                                float *data_out, const float *data_in) {

  data_in += size_x * size_y * blockIdx.z;
  data_out += size_x * size_y * blockIdx.z;

  const float *read_row;

  for (int x = threadIdx.x; x < size_x; x += blockDim.x) {
    MirroredIterator row_id;
    float sum = 0, diff = 0, t;

    // gridDim.y = height/2, center of the convolution is blockIdx.y*2
    row_id.init(gridDim.y*2, blockIdx.y*2 - 4);

    /*
    float sum =
      CDF97_ANALYSIS_LOWPASS_FILTER_0 * t4
      + CDF97_ANALYSIS_LOWPASS_FILTER_1 * (t3+t5)
      + CDF97_ANALYSIS_LOWPASS_FILTER_2 * (t2+t6)
      + CDF97_ANALYSIS_LOWPASS_FILTER_3 * (t1+t7)
      + CDF97_ANALYSIS_LOWPASS_FILTER_4 * (t0+t8);

    float diff =
      CDF97_ANALYSIS_HIGHPASS_FILTER_0 * t5
      + CDF97_ANALYSIS_HIGHPASS_FILTER_1 * (t4+t6)
      + CDF97_ANALYSIS_HIGHPASS_FILTER_2 * (t3+t7)
      + CDF97_ANALYSIS_HIGHPASS_FILTER_3 * (t2+t8);
    */

    read_row = data_in + size_x * (row_id++);
    t = read_row[x];
    sum += CDF97_ANALYSIS_LOWPASS_FILTER_4 * t;

    read_row = data_in + size_x * (row_id++);
    t = read_row[x];
    sum += CDF97_ANALYSIS_LOWPASS_FILTER_3 * t;

    read_row = data_in + size_x * (row_id++);
    t = read_row[x];
    sum += CDF97_ANALYSIS_LOWPASS_FILTER_2 * t;
    diff += CDF97_ANALYSIS_HIGHPASS_FILTER_3 * t;

    read_row = data_in + size_x * (row_id++);
    t = read_row[x];
    sum += CDF97_ANALYSIS_LOWPASS_FILTER_1 * t;
    diff += CDF97_ANALYSIS_HIGHPASS_FILTER_2 * t;

    read_row = data_in + size_x * (row_id++);
    t = read_row[x];
    sum += CDF97_ANALYSIS_LOWPASS_FILTER_0 * t;
    diff += CDF97_ANALYSIS_HIGHPASS_FILTER_1 * t;

    read_row = data_in + size_x * (row_id++);
    t = read_row[x];
    sum += CDF97_ANALYSIS_LOWPASS_FILTER_1 * t;
    diff += CDF97_ANALYSIS_HIGHPASS_FILTER_0 * t;

    read_row = data_in + size_x * (row_id++);
    t = read_row[x];
    sum += CDF97_ANALYSIS_LOWPASS_FILTER_2 * t;
    diff += CDF97_ANALYSIS_HIGHPASS_FILTER_1 * t;

    read_row = data_in + size_x * (row_id++);
    t = read_row[x];
    sum += CDF97_ANALYSIS_LOWPASS_FILTER_3 * t;
    diff += CDF97_ANALYSIS_HIGHPASS_FILTER_2 * t;

    read_row = data_in + size_x * (row_id++);
    t = read_row[x];
    sum += CDF97_ANALYSIS_LOWPASS_FILTER_4 * t;
    diff += CDF97_ANALYSIS_HIGHPASS_FILTER_3 * t;

    /*
    float sum =
      CDF97_ANALYSIS_LOWPASS_FILTER_0 * t4
      + CDF97_ANALYSIS_LOWPASS_FILTER_1 * (t3+t5)
      + CDF97_ANALYSIS_LOWPASS_FILTER_2 * (t2+t6)
      + CDF97_ANALYSIS_LOWPASS_FILTER_3 * (t1+t7)
      + CDF97_ANALYSIS_LOWPASS_FILTER_4 * (t0+t8);

    float diff =
      CDF97_ANALYSIS_HIGHPASS_FILTER_0 * t5
      + CDF97_ANALYSIS_HIGHPASS_FILTER_1 * (t4+t6)
      + CDF97_ANALYSIS_HIGHPASS_FILTER_2 * (t3+t7)
      + CDF97_ANALYSIS_HIGHPASS_FILTER_3 * (t2+t8);
    */

    data_out[size_x * blockIdx.y + x] = sum;
    data_out[size_x * (blockIdx.y + gridDim.y) + x] = diff;

  }
}


// do 'level_count' transforms in the y direction
void cdf97_v2(float *data, float *data_tmp,
              scu_wavelet::int3 size, int level_count,
              CudaTimer *transformTimer) {

  dim3 gridDim, blockDim(CDF97_V2_BLOCK_SIZE);
  gridDim.x = 1;
  gridDim.y = size.y / 2;
  gridDim.z = size.z;

  if (transformTimer) transformTimer->start();

  for (int i=0; i < level_count; i++) {
    cdf97_v2_kernel<<<gridDim, blockDim>>>
      (size.x, size.y, data_tmp, data);

    // copy the results back into 'data'
    gridDim.y *= 2;
    haar_v2_copyrows<<<gridDim, blockDim>>>(size.x, size.y, data, data_tmp);

    gridDim.y /= 4;
  }

  if (transformTimer) transformTimer->end();
}


void cdf97_3d_cuda_v2(float *data, float *data_tmp,
                      scu_wavelet::int3 &size, scu_wavelet::int3 stepCount,
                      bool inverse,
                      CudaTimer *transformTimer, CudaTimer *transposeTimer) {

  if (inverse) {
    cdf97_3d_cuda(data, data_tmp, size, stepCount, inverse,
                  transformTimer, transposeTimer);
    return;
  }

  CudaTimer mytimer("cdf97_3d_cuda_v2");
  mytimer.start();

  // change XYZ -> ZXY

  if (transposeTimer) transposeTimer->start();
  gpuTranspose3dBack(data_tmp, data, size);
  if (transposeTimer) transposeTimer->end();

  /*
  if (transformTimer) transformTimer->start();
  dim3 gridDim(1, size.y, size.z);
  dim3 blockDim(CDF97_3D_BLOCK_SIZE);
  size_t sharedMemSize = sizeof(float) * (size.x + 8);
    cdf97_3d_kernel<<<gridDim, blockDim, sharedMemSize>>>
      (data, size.x, stepCount.x);
  if (transformTimer) transformTimer->end();
  */
  cdf97_v2(data_tmp, data, size, stepCount.x, transformTimer);


  // ZXY -> XYZ
  if (transposeTimer) transposeTimer->start();
  gpuTranspose3dFwd(data, data_tmp, size);
  if (transposeTimer) transposeTimer->end();


  cdf97_v2(data, data_tmp, size, stepCount.y, transformTimer);

  // XYZ -> YZX
  if (transposeTimer) transposeTimer->start();
  gpuTranspose3dFwd(data_tmp, data, size);
  if (transposeTimer) transposeTimer->end();

  cdf97_v2(data_tmp, data, size, stepCount.z, transformTimer);

  // YZX -> ZXY
  if (transposeTimer) transposeTimer->start();
  gpuTranspose3dFwd(data, data_tmp, size);
  if (transposeTimer) transposeTimer->end();

  mytimer.end();
  mytimer.sync();
  mytimer.print();
}


/**
   Compute CDF97 wavelet transform in the Y direction by reading a tile
   of data and computing as many output rows as the tile can provide:

             4
     row 0   | 5
     row 1   | | 6
     row 2   | | | 7
     row 3   | | | | 8
     row 4   | | | | | 9
     row 5   | | | | | |
     row 6   | | | | | |
     row 7   | | | | | |
     row 8   | | | | | |
     row 9     | | | | |
     row 10      | | | |
     row 11        | | |
     row 12          | |
     row 13            |

     Next tile starts by computing row 10, so it starts reading at row 6.

   With 14 rows read, 6 rows (rows-8) can be computed,
   so this kernel will compute (rows-8) rows starting at
   blockIdx.y * (rows-8).
*/
__global__ void cdf97_v3_kernel(int size_x, int size_y,
                                float *data_out, const float *data_in) {

  // the tile may be larger than the thread block
  __shared__ float tile[CDF97_V3_TILE_WIDTH * CDF97_V3_TILE_HEIGHT];

  int first_col = blockIdx.x * CDF97_V3_TILE_WIDTH;
  int first_output_row = blockIdx.y * (CDF97_V3_TILE_HEIGHT - 8);

  // position data_in and data_out to the correct x and z offset
  data_in += size_x * size_y * blockIdx.z + first_col;
  data_out += size_x * size_y * blockIdx.z + first_col;

  int tile_width = min(size_x - first_col, CDF97_V3_TILE_WIDTH);
  int tile_height = min(size_y - first_output_row + 8, CDF97_V3_TILE_HEIGHT);
  
  int first_input_row = first_output_row - 4;
  int last_input_row =  first_output_row + 4;

  // read the tile into shared memory
  for (int tile_y = threadIdx.y; tile_y < tile_height; tile_y += blockDim.y) {

    // find the right row, since they will be mirrored
    int row_id = first_input_row + tile_y;

    // no thread divergence because these are just based on blockIdx.y
    if (first_input_row < 0 || last_input_row >= size_y) {
      row_id = MirroredIterator::getOffset(row_id, size_y);
    }
    
    for (int tile_x = threadIdx.x; tile_x < tile_width; tile_x += blockDim.x) {
      tile[tile_y * CDF97_V3_TILE_WIDTH + tile_x] =
        data_in[row_id * size_x + tile_x + first_col];
      /*
      printf("[%d,%d,%d %d,%d] tile %d,%d = %f\n",
             blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y,
             tile_x, tile_y, tile[tile_y * CDF97_V3_TILE_WIDTH + tile_x]);
      */
    }
  }

  /*
  if (blockIdx.x+blockIdx.y+blockIdx.z+threadIdx.x+threadIdx.y+threadIdx.z==0) {
    printf("Tile:\n");
    for (int y=0; y < tile_height; y++) {
      for (int x=0; x < tile_width; x++) {
        printf("%8.4f ", tile[y*CDF97_V3_TILE_WIDTH + x]);
      }
      printf("\n");
    }
  }
  */

  __syncthreads();

  if (threadIdx.x + threadIdx.y == 0) {
    printf("tile size %dx%d\n", tile_width, tile_height);
  }

  // compute each of the output rows
  for (int output_row_id = threadIdx.y;
       output_row_id < tile_height - 8;
       output_row_id += blockDim.y) {

    float *readp = tile + output_row_id*CDF97_V3_TILE_WIDTH + threadIdx.x;

    float sum = *readp * CDF97_ANALYSIS_LOWPASS_FILTER_4;

    readp += CDF97_V3_TILE_WIDTH;
    sum += *readp * CDF97_ANALYSIS_LOWPASS_FILTER_3;

    readp += CDF97_V3_TILE_WIDTH;
    sum += *readp * CDF97_ANALYSIS_LOWPASS_FILTER_2;
    float diff = *readp * CDF97_ANALYSIS_HIGHPASS_FILTER_3;

    readp += CDF97_V3_TILE_WIDTH;
    sum += *readp * CDF97_ANALYSIS_LOWPASS_FILTER_1;
    diff += *readp * CDF97_ANALYSIS_HIGHPASS_FILTER_2;

    readp += CDF97_V3_TILE_WIDTH;
    sum += *readp * CDF97_ANALYSIS_LOWPASS_FILTER_0;
    diff += *readp * CDF97_ANALYSIS_HIGHPASS_FILTER_1;

    readp += CDF97_V3_TILE_WIDTH;
    sum += *readp * CDF97_ANALYSIS_LOWPASS_FILTER_1;
    diff += *readp * CDF97_ANALYSIS_HIGHPASS_FILTER_0;

    readp += CDF97_V3_TILE_WIDTH;
    sum += *readp * CDF97_ANALYSIS_LOWPASS_FILTER_2;
    diff += *readp * CDF97_ANALYSIS_HIGHPASS_FILTER_1;

    readp += CDF97_V3_TILE_WIDTH;
    sum += *readp * CDF97_ANALYSIS_LOWPASS_FILTER_3;
    diff += *readp * CDF97_ANALYSIS_HIGHPASS_FILTER_2;

    readp += CDF97_V3_TILE_WIDTH;
    sum += *readp * CDF97_ANALYSIS_LOWPASS_FILTER_4;
    diff += *readp * CDF97_ANALYSIS_HIGHPASS_FILTER_3;

    size_t offset = (first_output_row + output_row_id) * size_x + threadIdx.x;
    if (threadIdx.x == 0)
      printf("Results for row %d (%f %f) to %d / %d\n", threadIdx.y, sum,
             diff, (int)offset, (int)(offset + (size_y>>1) * size_x));

    data_out[offset] = sum;
    data_out[offset + (size_y>>1) * size_x] = diff;

    /*
    if (blockIdx.x+blockIdx.y+blockIdx.z+threadIdx.x+threadIdx.z==0) {
      printf("[%d] sum = %f, diff = %f\n", threadIdx.y, sum, diff);
    }
    */

  }

}


// do 'level_count' transforms in the y direction
void cdf97_v3(float *data, float *data_tmp,
              scu_wavelet::int3 size, int level_count,
              CudaTimer *transformTimer) {

  dim3 gridDim, blockDim(CDF97_V3_BLOCK_WIDTH, CDF97_V3_BLOCK_HEIGHT);
  
  gridDim.x = (size.x-1) / blockDim.x + 1;
  gridDim.y = (size.y-1) / (CDF97_V3_TILE_HEIGHT - 8) + 1;
  gridDim.z = size.z;
  
  if (transformTimer) transformTimer->start();

  printDeviceArray(data, size.x, size.y, 1, "before cdf97_v3_kernel");
  cdf97_v3_kernel<<<gridDim, blockDim>>>
    (size.x, size.y, data_tmp, data);

  if (transformTimer) transformTimer->end();
}


void cdf97_3d_cuda_v3(float *data, float *data_tmp,
                      scu_wavelet::int3 &size, scu_wavelet::int3 stepCount,
                      bool inverse,
                      CudaTimer *transformTimer, CudaTimer *transposeTimer) {
  
  cdf97_v3(data, data_tmp, size, 1, transformTimer);
  printDeviceArray(data_tmp, size, "data_tmp after v3");

  cdf97_v2(data, data_tmp, size, 1, transformTimer);
  printDeviceArray(data_tmp, size, "data_tmp after v2");
}
