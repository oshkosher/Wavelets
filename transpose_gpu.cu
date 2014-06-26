#include "transpose_gpu.h"

#define TX_BLOCK_SIZE 16

/* Given the width and height of the input matrix, transpose it
   into the given output matrix.

   blockSize is the size of the square thread block that will process
   each tile.
*/
void gpuTranspose(int width, int height, float *matrix_d, float *matrixTx_d,
                  cudaStream_t stream) {
                  

  dim3 gridSize, blockSize(TX_BLOCK_SIZE, TX_BLOCK_SIZE);
  gridSize.x = ceil(width-1) / TX_BLOCK_SIZE + 1;
  gridSize.y = ceil(height-1) / TX_BLOCK_SIZE + 1;

  gpuTransposeKernel<<<gridSize, blockSize, 0, stream>>>
    (width, height, matrix_d, matrixTx_d);
}


/* gpuTransposeTiledKernel
   Copy one 16x16 tile (one element per thread) from the first matrix
   into shared memory, then from shared memory to the second matrix.
   Be sure to structure the global memory accesses so that consecutive
   threads access consecutive memory.
*/
__global__ void gpuTransposeKernel(int width, int height, float *matrix,
                                   float *matrixTx) {

  __shared__ float cache[TX_BLOCK_SIZE+1][TX_BLOCK_SIZE+1];
  int tileTop = blockIdx.y * blockDim.y;
  int tileLeft = blockIdx.x * blockDim.x;

  // Since the matrix size is a multiple of TILE_SIZE, as long as the 
  // upper left corner of the tile is in the matrix, the entire tile will be.
  // if (tileTop >= height || tileLeft >= width) return;

  int row = tileTop + threadIdx.y;
  int col = tileLeft + threadIdx.x;
  
  if (row < width && col < height)
    cache[threadIdx.y][threadIdx.x] = matrix[row * width + col];

  // Sync is necessary because the thread that wrote to a shared memory
  // entry won't always be the one that reads from it.
  __syncthreads();

  row = tileLeft + threadIdx.y;
  col = tileTop + threadIdx.x;
  if (row < width && col < height)
    matrixTx[row * height + col] = cache[threadIdx.x][threadIdx.y];
}
