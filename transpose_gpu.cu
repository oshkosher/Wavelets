#include <cstdio>
#include <cassert>
#include "cuda.h"
#include "transpose_gpu.h"

#define TX_BLOCK_SIZE 16


/* gpuTransposeTiledKernel
   Copy one 16x16 tile (one element per thread) from the first matrix
   into shared memory, then from shared memory to the second matrix.
   Be sure to structure the global memory accesses so that consecutive
   threads access consecutive memory.
*/
template<typename NUM>
__global__ void gpuTransposeSquareKernel(int fullWidth, int transposeSize,
                                         NUM *matrix, NUM *matrixTx) {
                                   

  __shared__ NUM cache[TX_BLOCK_SIZE+1][TX_BLOCK_SIZE+1];
  int tileTop = blockIdx.y * blockDim.y;
  int tileLeft = blockIdx.x * blockDim.x;

  // Since the matrix size is a multiple of TILE_SIZE, as long as the 
  // upper left corner of the tile is in the matrix, the entire tile will be.
  // if (tileTop >= height || tileLeft >= width) return;

  int row = tileTop + threadIdx.y;
  int col = tileLeft + threadIdx.x;
  
  if (row < transposeSize && col < transposeSize)
    cache[threadIdx.y][threadIdx.x] = matrix[row * fullWidth + col];

  // Sync is necessary because the thread that wrote to a shared memory
  // entry won't always be the one that reads from it.
  __syncthreads();

  row = tileLeft + threadIdx.y;
  col = tileTop + threadIdx.x;
  if (row < transposeSize && col < transposeSize)
    matrixTx[row * fullWidth + col] = cache[threadIdx.x][threadIdx.y];
}


/* Given the width and height of the input matrix, transpose it
   into the given output matrix.

   blockSize is the size of the square thread block that will process
   each tile.
*/
template<typename NUM>
void gpuTransposeSquareInternal(int fullWidth, int transposeSize,
                          NUM *matrix_d, NUM *matrixTx_d,
                          cudaStream_t stream) {

  dim3 gridSize, blockSize(TX_BLOCK_SIZE, TX_BLOCK_SIZE);
  gridSize.x = (unsigned)(ceil((float)(transposeSize-1)) / TX_BLOCK_SIZE + 1);
  gridSize.y = gridSize.x;

  gpuTransposeSquareKernel<<<gridSize, blockSize, 0, stream>>>
    (fullWidth, transposeSize, matrix_d, matrixTx_d);
}


void gpuTransposeSquare(int fullWidth, int transposeSize,
                        float *matrix_d, float *matrixTx_d,
                        cudaStream_t stream) {
  gpuTransposeSquareInternal(fullWidth, transposeSize, matrix_d, matrixTx_d,
                             stream);
}


// double support was added in version 1.3
#if !defined(__CUDA_ARCH__) || (__CUDA_ARCH__ >= 130)
void gpuTransposeSquare(int fullWidth, int transposeSize,
                        double *matrix_d, double *matrixTx_d,
                        cudaStream_t stream) {
  gpuTransposeSquareInternal(fullWidth, transposeSize, matrix_d, matrixTx_d,
                             stream);
}
#endif



/* gpuTransposeTiledKernel
   Copy one 16x16 tile (one element per thread) from the first matrix
   into shared memory, then from shared memory to the second matrix.

   Unlike gpuTransposeSquareKernel, this one supports rectangular arrays.
*/
__global__ void gpuTransposeKernel(float *src, float *dest,
                                   int width, int height) {

  __shared__ float cache[TX_BLOCK_SIZE+1][TX_BLOCK_SIZE+1];

  // get the upper left corner of this tile (same for all threads in this block)
  int tileLeft = blockIdx.x * blockDim.x;
  int tileTop = blockIdx.y * blockDim.y;

  /*
  if (threadIdx.x == 0 && threadIdx.y== 0) {
    printf("Block %d,%d at %d,%d\n", 
           blockIdx.x, blockIdx.y, tileLeft, tileTop);
  }
  */

  // Since the number of thread blocks was set to void creating unused
  // ones, this shouldn't happen.
  if (tileTop >= height || tileLeft >= width) {
    printf("Block %d,%d starts out of range at %d,%d\n",
           blockIdx.x, blockIdx.y, tileLeft, tileTop);
    return;
  }
  
  // get read coordinate for this thread
  int col = tileLeft + threadIdx.x;
  int row = tileTop + threadIdx.y;
  
  if (row < height && col < width) {
    cache[threadIdx.y][threadIdx.x] = src[row * width + col];
    /*
    printf("Thread %d,%d,%d,%d read %d,%d to %d,%d\n", 
           blockIdx.y, blockIdx.x, threadIdx.y, threadIdx.x,
           row, col, threadIdx.y, threadIdx.x);
    */
  }

  // Sync is necessary because the thread that wrote to a shared memory
  // entry won't always be the one that reads from it.
  __syncthreads();

  // swap what width and height mean

  row = tileLeft + threadIdx.y;
  col = tileTop + threadIdx.x;
  if (row < width && col < height) {
    dest[row * height + col] = cache[threadIdx.x][threadIdx.y];
    /*
    printf("Thread %d,%d,%d,%d write %d,%d to %d,%d\n", 
           blockIdx.y, blockIdx.x, threadIdx.y, threadIdx.x,
           threadIdx.x, threadIdx.y, row, col);
    */
  }
}


/* Transpose rectangular 2-d arrays. */
void gpuTranspose(float *src_d, float *dest_d, int width, int height, 
                  cudaStream_t stream) {
  
  dim3 gridSize, blockSize(TX_BLOCK_SIZE, TX_BLOCK_SIZE);
  gridSize.x = (width + TX_BLOCK_SIZE - 1) / TX_BLOCK_SIZE;
  gridSize.y = (height + TX_BLOCK_SIZE - 1) / TX_BLOCK_SIZE;

  assert(gridSize.x <= 65535 && gridSize.y <= 65535);

  // printf("Transpose %dx%d with %dx%d tiles of size %dx%d\n",
  // width, height, gridSize.x, gridSize.y, TX_BLOCK_SIZE, TX_BLOCK_SIZE);

  gpuTransposeKernel<<<gridSize, blockSize, 0, stream>>>
    (src_d, dest_d, width, height);

}

