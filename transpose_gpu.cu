#include <cstdio>
#include <cassert>
#include "cuda.h"
#include "transpose_gpu.h"

#define TX_BLOCK_SIZE 32
#define TX_TILE_SIZE 32


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


/*
  Try using tiles larger than a thread block.
  This is slower due to large shared memory requirements and
  lower occupancy.
*/
__global__ void gpuTransposeKernelLargeTiles(float *dest, const float *src,
                                             int width, int height) {

  __shared__ float cache[TX_TILE_SIZE][TX_TILE_SIZE+1];

  int tileTop, tileLeft, tileRow, tileCol, row, col;
  
  // get the upper left corner of this tile (same for all threads in this block)
  tileTop = blockIdx.y * blockDim.y;
  int endRow = min(TX_TILE_SIZE, height - tileTop);

  while (tileTop < height) {

    tileLeft = blockIdx.x * blockDim.x;
    int endCol = min(TX_TILE_SIZE, width - tileLeft);

    while (tileLeft < width) {
  
      // tileRow is the y-coordinate within the tile
      tileRow = threadIdx.y;
  
      // within a tile, copy multiple elements to the cache
      while (tileRow < endRow) {

        // tileCol is the x-coordinate within the tile
        tileCol = threadIdx.x;

        while (tileCol < endCol) {
          row = tileTop + tileRow;
          col = tileLeft + tileCol;
          cache[tileRow][tileCol] = src[row * width + col];
          tileCol += blockDim.x;
        }

        tileRow += blockDim.y;
      }

      // Sync is necessary because the thread that wrote to a shared memory
      // entry won't always be the one that reads from it.
      __syncthreads();

      // swap what width and height mean
      
      tileCol = threadIdx.y;
      while (tileCol < endRow) {  // tileLeft + y < width  vs  y < height - tileTop
        
        tileRow = threadIdx.x;
        while (tileRow < endCol) {
          
          row = tileLeft + tileCol;
          col = tileTop + tileRow;
          dest[row * height + col] = cache[tileRow][tileCol];
          
          tileRow += blockDim.x;
        }

        tileCol += blockDim.y;
      }
          
      /*
      row = tileLeft + threadIdx.y;
      col = tileTop + threadIdx.x;
      if (row < width && col < height) {
        dest[row * height + col] = cache[threadIdx.x][threadIdx.y];
      }
      */

      tileLeft += blockDim.x * gridDim.x;
    }

    tileTop += blockDim.y * gridDim.y;
  }

}



/* gpuTransposeTiledKernel
   Copy one 16x16 tile (one element per thread) from the first matrix
   into shared memory, then from shared memory to the second matrix.

   Unlike gpuTransposeSquareKernel, this one supports rectangular arrays.

   Naive: 49.7ms
   Tiled/shared memory: 14.5ms
   Avoid bank conflicts with [n][n+1]: 8.0ms
*/
__global__ void gpuTransposeKernel(float *dest, const float *src,
                                   int width, int height) {

  __shared__ float cache[TX_BLOCK_SIZE][TX_BLOCK_SIZE+1];

  // get the upper left corner of this tile (same for all threads in this block)
  int tileTop = blockIdx.y * blockDim.y;

  while (tileTop < height) {

    int tileLeft = blockIdx.x * blockDim.x;

    while (tileLeft < width) {
  
      // get read coordinate for this thread
      int col = tileLeft + threadIdx.x;
      int row = tileTop + threadIdx.y;
  
      if (row < height && col < width) {
        cache[threadIdx.y][threadIdx.x] = src[row * width + col];
      }

      // Sync is necessary because the thread that wrote to a shared memory
      // entry won't always be the one that reads from it.
      __syncthreads();

      // swap what width and height mean

      row = tileLeft + threadIdx.y;
      col = tileTop + threadIdx.x;
      if (row < width && col < height) {
        dest[row * height + col] = cache[threadIdx.x][threadIdx.y];
      }

      tileLeft += blockDim.x * gridDim.x;
    }

    tileTop += blockDim.y * gridDim.y;
  }

}


__global__ void gpuTransposeKernelNaive(float *dest_d, float *src_d, int width, int height) {

  int y = blockIdx.y * blockDim.y + threadIdx.y;

  while (y < height) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    while (x < width) {
      dest_d[x * height + y] = src_d[y * width + x];
      x += blockDim.x * gridDim.x;
    }
    y += blockDim.y * gridDim.y;
  }
}


/* Transpose rectangular 2-d arrays.
 */
void gpuTranspose(float *dest_d, float *src_d, int width, int height, 
                  cudaStream_t stream) {
  
  dim3 gridSize, blockSize(TX_BLOCK_SIZE, TX_BLOCK_SIZE);
  gridSize.x = 1;
  // gridSize.x = (width - 1) / TX_BLOCK_SIZE + 1;
  gridSize.y = (height - 1) / TX_BLOCK_SIZE + 1;

  assert(gridSize.x <= 65535 && gridSize.y <= 65535);

  // printf("Transpose %dx%d with %dx%d tiles of size %dx%d\n",
  // width, height, gridSize.x, gridSize.y, TX_BLOCK_SIZE, TX_BLOCK_SIZE);

  gpuTransposeKernel<<<gridSize, blockSize, 0, stream>>>
    (dest_d, src_d, width, height);
  // gpuTransposeKernelNaive<<<gridSize, blockSize, 0, stream>>>
  // (dest_d, src_d, width, height);

}


void gpuTranspose3dFwd(float *dest, float *src, scu_wavelet::int3 &size,
                       cudaStream_t stream) {
  gpuTranspose(dest, src, size.x, size.y*size.z, stream);
  size.rotateFwd();
}


void gpuTranspose3dBack(float *dest, float *src, scu_wavelet::int3 &size,
                       cudaStream_t stream) {

  gpuTranspose(dest, src, size.x*size.y, size.z, stream);
  size.rotateBack();

}

