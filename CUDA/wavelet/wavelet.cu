/*
 *
 *  Based on Nvidia convolution separable example.
 *
 *
 * Copyright 1993-2014 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#include <cstdio>
#include <cuda.h>
#include <assert.h>
#include "../../cucheck.h"
#include "wavelet.h"

////////////////////////////////////////////////////////////////////////////////
// Convolution kernel storage
////////////////////////////////////////////////////////////////////////////////
__constant__ float c_Kernel[KERNEL_LENGTH*2];


////////////////////////////////////////////////////////////////////////////////
// Row convolution with Low and Hi pass filter
////////////////////////////////////////////////////////////////////////////////
#define ROWS_BLOCKDIM_X 16
#define ROWS_BLOCKDIM_Y 16
#define	ROWS_RESULT_STEPS 1 //8
#define ROWS_HALO_STEPS 1

__global__ void convolutionRowsMirrorHiLoKernel(float *d_Dst, float *d_Src, int imageW, int imageH, int pitch) {
    __shared__ float s_Data[ROWS_BLOCKDIM_Y][(ROWS_RESULT_STEPS + 2 * ROWS_HALO_STEPS) * ROWS_BLOCKDIM_X];

    //Offset to the left halo edge
    const int baseX = (blockIdx.x * ROWS_RESULT_STEPS - ROWS_HALO_STEPS) * ROWS_BLOCKDIM_X + threadIdx.x;
    const int baseY = blockIdx.y * ROWS_BLOCKDIM_Y + threadIdx.y;
	const int fidx  = threadIdx.x % 2;

    d_Src += baseY * pitch + baseX;
	const int half = (baseY * pitch + (blockIdx.x * ROWS_BLOCKDIM_X + threadIdx.x))/2;
	d_Dst += half+(fidx*(imageH*pitch)/2) - ROWS_HALO_STEPS * ROWS_BLOCKDIM_X;

    //Load main data
#pragma unroll

    for (int i = ROWS_HALO_STEPS; i < ROWS_HALO_STEPS + ROWS_RESULT_STEPS; i++)
    {
        s_Data[threadIdx.y][threadIdx.x + i * ROWS_BLOCKDIM_X] = d_Src[i * ROWS_BLOCKDIM_X];
    }

    //Load left halo
#pragma unroll

    for (int i = 0; i < ROWS_HALO_STEPS; i++)
    {
		// If HALO is > 1 maybe d_Src[i * ROWS_BLOCKDIM_X - baseX*2]; is not correct for every ROW_HALO_STEP
        s_Data[threadIdx.y][threadIdx.x + i * ROWS_BLOCKDIM_X] = (baseX >= -i * ROWS_BLOCKDIM_X) ? d_Src[i * ROWS_BLOCKDIM_X] : d_Src[i * ROWS_BLOCKDIM_X - baseX*2];
    }

    //Load right halo
#pragma unroll

    for (int i = ROWS_HALO_STEPS + ROWS_RESULT_STEPS; i < ROWS_HALO_STEPS + ROWS_RESULT_STEPS + ROWS_HALO_STEPS; i++)
    {
		s_Data[threadIdx.y][threadIdx.x + i * ROWS_BLOCKDIM_X] = (imageW - baseX > i * ROWS_BLOCKDIM_X) ? d_Src[i * ROWS_BLOCKDIM_X] : d_Src[i * ROWS_BLOCKDIM_X - (threadIdx.x+1)*2];
    }

    //Compute and store results
    __syncthreads();
#pragma unroll
    for (int i = ROWS_HALO_STEPS; i < ROWS_HALO_STEPS + ROWS_RESULT_STEPS; i++)
    {
        float sum = 0;

#pragma unroll

        for (int j = -KERNEL_RADIUS; j <= KERNEL_RADIUS; j++)
        {
			sum += c_Kernel[fidx * KERNEL_LENGTH + KERNEL_RADIUS + j] * s_Data[threadIdx.y][threadIdx.x + i * ROWS_BLOCKDIM_X + j];
        }
		d_Dst[i * ROWS_BLOCKDIM_X] = sum;
    }
}

__global__ void invConvolutionRowsMirrorHiLoKernel(float *d_Dst, float *d_Src, int imageW, int imageH, int pitch) {
    __shared__ float s_Data[ROWS_BLOCKDIM_Y][(ROWS_RESULT_STEPS + 2 * ROWS_HALO_STEPS) * ROWS_BLOCKDIM_X];

    //Offset to the left halo edge
    const int baseX = (blockIdx.x * ROWS_RESULT_STEPS - ROWS_HALO_STEPS) * ROWS_BLOCKDIM_X + threadIdx.x;
    const int baseY = blockIdx.y * ROWS_BLOCKDIM_Y + threadIdx.y;
	const int fidx  = threadIdx.x % 2;

	const int half = (baseY * pitch + (blockIdx.x * ROWS_BLOCKDIM_X + threadIdx.x))/2;
	d_Src += half+(fidx*(imageH*pitch)/2) - ROWS_HALO_STEPS * ROWS_BLOCKDIM_X;
	d_Dst += baseY * pitch + baseX;

    //Load main data
#pragma unroll

    for (int i = ROWS_HALO_STEPS; i < ROWS_HALO_STEPS + ROWS_RESULT_STEPS; i++)
    {
        s_Data[threadIdx.y][threadIdx.x + i * ROWS_BLOCKDIM_X] = d_Src[i * ROWS_BLOCKDIM_X];
    }

    //Load left halo
#pragma unroll

    for (int i = 0; i < ROWS_HALO_STEPS; i++)
    {
		// If HALO is > 1 maybe d_Src[i * ROWS_BLOCKDIM_X - baseX*2]; is not correct for every ROW_HALO_STEP
        s_Data[threadIdx.y][threadIdx.x + i * ROWS_BLOCKDIM_X] = (baseX >= -i * ROWS_BLOCKDIM_X) ? d_Src[i * ROWS_BLOCKDIM_X] : d_Src[i * ROWS_BLOCKDIM_X - baseX*2];
    }

    //Load right halo
#pragma unroll

    for (int i = ROWS_HALO_STEPS + ROWS_RESULT_STEPS; i < ROWS_HALO_STEPS + ROWS_RESULT_STEPS + ROWS_HALO_STEPS; i++)
    {
		s_Data[threadIdx.y][threadIdx.x + i * ROWS_BLOCKDIM_X] = (imageW - baseX > i * ROWS_BLOCKDIM_X) ? d_Src[i * ROWS_BLOCKDIM_X] : d_Src[i * ROWS_BLOCKDIM_X - (threadIdx.x+1)*2];
    }

    //Compute and store results
    __syncthreads();
#pragma unroll

    for (int i = ROWS_HALO_STEPS; i < ROWS_HALO_STEPS + ROWS_RESULT_STEPS; i++)
    {
        float sum = 0;

#pragma unroll

        for (int j = -KERNEL_RADIUS; j <= KERNEL_RADIUS; j++)
        {
			sum += c_Kernel[fidx * KERNEL_LENGTH + KERNEL_RADIUS + j] * s_Data[threadIdx.y][threadIdx.x + i * ROWS_BLOCKDIM_X + j];
        }

		d_Dst[i * ROWS_BLOCKDIM_X] = sum;
    }
}

void fwt_1D(float **data, const unsigned level, const unsigned nx, const unsigned ny) {
    assert(ROWS_BLOCKDIM_X * ROWS_HALO_STEPS >= KERNEL_RADIUS);
    assert(nx % (ROWS_RESULT_STEPS * ROWS_BLOCKDIM_X) == 0);
    assert(ny % ROWS_BLOCKDIM_Y == 0);

	const int mem_size = nx*ny*sizeof(float);

	float *data1, *data2, *aux;
	data1 = *data;
	cudaMalloc(&data2, mem_size);

	unsigned w = nx;

    dim3 blocks(nx / (ROWS_RESULT_STEPS * ROWS_BLOCKDIM_X), ny / ROWS_BLOCKDIM_Y);
    dim3 threads(ROWS_BLOCKDIM_X, ROWS_BLOCKDIM_Y);
	convolutionRowsMirrorHiLoKernel<<<blocks, threads>>>(data2, data1, w, ny, w);
	CUCHECK(cudaGetLastError());


	for (unsigned i = 1; i < level; i++) {
		blocks.x /= 2;
		w /= 2;

		aux = data2;
		data2 = data1;
		data1 = aux;

		cudaMemcpy(data2+w*ny, data1+w*ny, w*ny*sizeof(float), cudaMemcpyDeviceToDevice);

		convolutionRowsMirrorHiLoKernel<<<blocks, threads>>>(data2, data1, w, ny, w);
                CUCHECK(cudaGetLastError());
	}
    
	*data = data2;
	cudaFree(data1);

	printf("Rows fwt_1D: %s\n",cudaGetErrorString(cudaGetLastError()));
}

void iwt_1D(float **data, const unsigned level, const unsigned nx, const unsigned ny) {
    assert(ROWS_BLOCKDIM_X * ROWS_HALO_STEPS >= KERNEL_RADIUS);
    assert(nx % (ROWS_RESULT_STEPS * ROWS_BLOCKDIM_X) == 0);
    assert(ny % ROWS_BLOCKDIM_Y == 0);

	const int mem_size = nx*ny*sizeof(float);

	float *data1, *data2, *aux;
	data1 = *data;
	cudaMalloc(&data2, mem_size);

	unsigned w = nx >> (level-1);

    dim3 blocks(w / (ROWS_RESULT_STEPS * ROWS_BLOCKDIM_X), ny / ROWS_BLOCKDIM_Y);
    dim3 threads(ROWS_BLOCKDIM_X, ROWS_BLOCKDIM_Y);
	invConvolutionRowsMirrorHiLoKernel<<<blocks, threads>>>(data2, data1, w, ny, w);
        CUCHECK(cudaGetLastError());

	for (unsigned i = 1; i < level; i++) {
		cudaMemcpy(data2+w*ny, data1+w*ny, (nx-w)*ny*sizeof(float), cudaMemcpyDeviceToDevice);
		
		blocks.x *= 2;
		w *= 2;

		aux = data2;
		data2 = data1;
		data1 = aux;

		invConvolutionRowsMirrorHiLoKernel<<<blocks, threads>>>(data2, data1, w, ny, w);
                CUCHECK(cudaGetLastError());
	}
    
	*data = data2;
	cudaFree(data1);

	printf("Rows iwt_1D: %s\n",cudaGetErrorString(cudaGetLastError()));
}

////////////////////////////////////////////////////////////////////////////////
// Transpose
////////////////////////////////////////////////////////////////////////////////

const int TILE_DIM = 32;
const int BLOCK_ROWS = 8;

__global__ void transposeDiagonal(float *odata, const float *idata, int width, int height) {
	__shared__ float tile[TILE_DIM][TILE_DIM+1];

	int blockIdx_x, blockIdx_y;

	// diagonal reordering
	if (width == height) {
		blockIdx_y = blockIdx.x;
		blockIdx_x = (blockIdx.x+blockIdx.y)%gridDim.x;
	} else {
		int bid = blockIdx.x + gridDim.x*blockIdx.y;
		blockIdx_y = bid%gridDim.y;
		blockIdx_x = ((bid/gridDim.y)+blockIdx_y)%gridDim.x;
	}

	int xIndex = blockIdx_x*TILE_DIM + threadIdx.x;
	int yIndex = blockIdx_y*TILE_DIM + threadIdx.y;
	int index_in = xIndex + (yIndex)*width;

	xIndex = blockIdx_y*TILE_DIM + threadIdx.x;
	yIndex = blockIdx_x*TILE_DIM + threadIdx.y;
	int index_out = xIndex + (yIndex)*height;

	for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
		tile[threadIdx.y+i][threadIdx.x] = idata[index_in+i*width];
	}

	__syncthreads();

	for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
		odata[index_out+i*height] = tile[threadIdx.x][threadIdx.y+i];
	}
} 
void transpose(float *tdata, const float *idata, const unsigned nx, const unsigned ny) {

	dim3 grid(nx/TILE_DIM, ny/TILE_DIM);
	dim3 threads(TILE_DIM,BLOCK_ROWS); 
	transposeDiagonal<<<grid, threads>>>(tdata, idata, nx, ny);
        CUCHECK(cudaDeviceSynchronize());
}

extern "C" void setUpFilter(const float *filter){

    cudaMemcpyToSymbol(c_Kernel, filter, KERNEL_LENGTH*2 * sizeof(float));

	printf("Setup: %s\n",cudaGetErrorString(cudaGetLastError()));
}

extern "C" void fwt_1D_GPU(float *data, const unsigned level, const unsigned nx, const unsigned ny) {
	const int mem_size = nx*ny*sizeof(float);

	float *d_idata;
	cudaMalloc(&d_idata, mem_size);

	cudaMemcpy(d_idata, data, mem_size, cudaMemcpyHostToDevice);

	fwt_1D(&d_idata, level, nx, ny);

	cudaMemcpy(data, d_idata, mem_size, cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();
	printf("FWT_GPU: %s\n",cudaGetErrorString(cudaGetLastError()));

	cudaFree(d_idata);
}

extern "C" void iwt_1D_GPU(float *data, const unsigned level, const unsigned nx, const unsigned ny) {
	const int mem_size = nx*ny*sizeof(float);

	float *d_idata;
	cudaMalloc(&d_idata, mem_size);

	cudaMemcpy(d_idata, data, mem_size, cudaMemcpyHostToDevice);

	iwt_1D(&d_idata, level, nx, ny);

	cudaMemcpy(data, d_idata, mem_size, cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();
	printf("IWT_GPU: %s\n",cudaGetErrorString(cudaGetLastError()));

	cudaFree(d_idata);
}

// data_dest: set this to the device address where the output data resides.
// If this is not equal to 'data', then 'data' has been freed.
extern "C" void wavelet_cuda_3d_fwd(float *data, const unsigned nx, const unsigned ny, const unsigned nz, const unsigned lvlx, const unsigned lvly, const unsigned lvlz, bool data_is_on_gpu) {
	const int mem_size = nx*ny*nz*sizeof(float);

	float *d_idata, *d_tdata;
	cudaMalloc(&d_tdata, mem_size);
        cudaMalloc(&d_idata, mem_size);

        if (data_is_on_gpu) {
          cudaMemcpy(d_idata, data, mem_size, cudaMemcpyDeviceToDevice);
        } else {
          cudaMemcpy(d_idata, data, mem_size, cudaMemcpyHostToDevice);
        }

	fwt_1D(&d_idata, lvlx, nx, ny*nz);

	transpose(d_tdata, d_idata, nx, ny*nz);

	fwt_1D(&d_tdata, lvly, ny, nz*nx);

	transpose(d_idata, d_tdata, ny, nz*nx);

	fwt_1D(&d_idata, lvlz, nz, nx*ny);

        if (data_is_on_gpu) {
          cudaMemcpy(data, d_idata, mem_size, cudaMemcpyDeviceToDevice);
        } else {
          cudaMemcpy(data, d_idata, mem_size, cudaMemcpyDeviceToHost);
        }
        cudaFree(d_idata);

	//cudaDeviceSynchronize();
	//printf("comp: %s\n",cudaGetErrorString(cudaGetLastError()));

	cudaFree(d_tdata);
}

extern "C" void wavelet_cuda_3d_back(float *data, const unsigned nx, const unsigned ny, const unsigned nz, const unsigned lvlx, const unsigned lvly, const unsigned lvlz) {
	const int mem_size = nx*ny*nz*sizeof(float);

	float *d_idata, *d_tdata;
	cudaMalloc(&d_idata, mem_size);
	cudaMalloc(&d_tdata, mem_size);

	cudaMemcpy(d_idata, data, mem_size, cudaMemcpyHostToDevice);

	iwt_1D(&d_idata, lvlz, nz, nx*ny);

	transpose(d_tdata, d_idata, nz*nx, ny);

	iwt_1D(&d_tdata, lvly, ny, nz*nx);

	transpose(d_idata, d_tdata, ny*nz, nx);

	iwt_1D(&d_idata, lvlx, nx, ny*nz);


	cudaMemcpy(data, d_idata, mem_size, cudaMemcpyDeviceToHost);

	//cudaDeviceSynchronize();
	//printf("ucomp: %s\n",cudaGetErrorString(cudaGetLastError()));

	cudaFree(d_idata);
	cudaFree(d_tdata);
}
