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

#include "wavelet.h"

////////////////////////////////////////////////////////////////////////////////
// Convolution kernel storage
////////////////////////////////////////////////////////////////////////////////
__constant__ float c_Kernel[KERNEL_LENGTH*2];


////////////////////////////////////////////////////////////////////////////////
// Row convolution with Low and Hi pass filter
////////////////////////////////////////////////////////////////////////////////
#define ROWS_BLOCKDIM_X 16
#define ROWS_BLOCKDIM_Y 2
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

void fwt_1D(float *data, const unsigned level, const unsigned nx, const unsigned ny) {
    assert(ROWS_BLOCKDIM_X * ROWS_HALO_STEPS >= KERNEL_RADIUS);
    assert(d_w % (ROWS_RESULT_STEPS * ROWS_BLOCKDIM_X) == 0);
    assert(d_h % ROWS_BLOCKDIM_Y == 0);

    dim3 blocks(nx / (ROWS_RESULT_STEPS * ROWS_BLOCKDIM_X), ny / ROWS_BLOCKDIM_Y);
    dim3 threads(ROWS_BLOCKDIM_X, ROWS_BLOCKDIM_Y);
	unsigned w = nx;
	for (unsigned i = 0; i < level; i++) {
		convolutionRowsMirrorHiLoKernel<<<blocks, threads>>>(data, data, w, ny, w);

		blocks.x /= 2;
		w /= 2;
	}
    
	//cudaMemcpy(output, data, d_w * d_h * sizeof(float), cudaMemcpyDeviceToHost);

	//printf("Rows: %s\n",cudaGetErrorString(cudaGetLastError()));
}

void iwt_1D(float *data, const unsigned level, const unsigned nx, const unsigned ny) {
    assert(ROWS_BLOCKDIM_X * ROWS_HALO_STEPS >= KERNEL_RADIUS);
    assert(d_w % (ROWS_RESULT_STEPS * ROWS_BLOCKDIM_X) == 0);
    assert(d_h % ROWS_BLOCKDIM_Y == 0);

    dim3 blocks(nx / (ROWS_RESULT_STEPS * ROWS_BLOCKDIM_X), ny / ROWS_BLOCKDIM_Y);
    dim3 threads(ROWS_BLOCKDIM_X, ROWS_BLOCKDIM_Y);
	unsigned w = nx;
	for (unsigned i = 0; i < level; i++) {
		invConvolutionRowsMirrorHiLoKernel<<<blocks, threads>>>(data, data, w, ny, w);

		blocks.x /= 2;
		w /= 2;
	}
    
	//cudaMemcpy(output, data, d_w * d_h * sizeof(float), cudaMemcpyDeviceToHost);

	//printf("Rows: %s\n",cudaGetErrorString(cudaGetLastError()));
}

////////////////////////////////////////////////////////////////////////////////
// Transpose
////////////////////////////////////////////////////////////////////////////////

// Copyright 2012 NVIDIA Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

const int TILE_DIM = 32;
const int BLOCK_ROWS = 8;

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline cudaError_t checkCuda(cudaError_t result) {
#if defined(DEBUG) || defined(_DEBUG)
	if (result != cudaSuccess) {
		fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
		assert(result == cudaSuccess);
	}
#endif
	return result;
}

// No bank-conflict transpose
// Same as transposeCoalesced except the first tile dimension is padded
// to avoid shared memory bank conflicts.
__global__ void transposeNoBankConflicts(float *odata, const float *idata) {
	__shared__ float tile[TILE_DIM][TILE_DIM+1];
	
	int x = blockIdx.x * TILE_DIM + threadIdx.x;
	int y = blockIdx.y * TILE_DIM + threadIdx.y;
	int width = gridDim.x * TILE_DIM;

	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
		tile[threadIdx.y+j][threadIdx.x] = idata[(y+j)*width + x];

	__syncthreads();

	x = blockIdx.y * TILE_DIM + threadIdx.x; // transpose block offset
	y = blockIdx.x * TILE_DIM + threadIdx.y;

	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
		odata[(y+j)*width + x] = tile[threadIdx.x][threadIdx.y + j];
}


void transpose(float *tdata, const float *idata, const unsigned nx, const unsigned ny) {
	//const int mem_size = nx*ny*sizeof(float);

	//float *d_idata, *d_tdata;
	//checkCuda( cudaMalloc(&d_idata, mem_size) );
	//checkCuda( cudaMalloc(&d_tdata, mem_size) );

	//checkCuda( cudaMemcpy(d_idata, idata, mem_size, cudaMemcpyHostToDevice) );

	dim3 dimGrid(nx/TILE_DIM, ny/TILE_DIM, 1);
	dim3 dimBlock(TILE_DIM, BLOCK_ROWS, 1);

	transposeNoBankConflicts<<<dimGrid, dimBlock>>>(tdata, idata);
	//checkCuda( cudaMemcpy(tdata, d_tdata, mem_size, cudaMemcpyDeviceToHost) );
}

extern "C" void setUpFilter(const float *filter){

    cudaMemcpyToSymbol(c_Kernel, filter, KERNEL_LENGTH*2 * sizeof(float));

	//printf("Setup: %s\n",cudaGetErrorString(cudaGetLastError()));
}


extern "C" void comp(float *data, const unsigned nx, const unsigned ny, const unsigned nz, const unsigned lvlx, const unsigned lvly, const unsigned lvlz) {
	const int mem_size = nx*ny*nz*sizeof(float);

	float *d_idata, *d_tdata;
	checkCuda( cudaMalloc(&d_idata, mem_size) );
	checkCuda( cudaMalloc(&d_tdata, mem_size) );

	checkCuda( cudaMemcpy(d_idata, data, mem_size, cudaMemcpyHostToDevice) );

	fwt_1D(d_idata, lvlx, nx, ny*nz);

	transpose(d_tdata, d_idata, nx, ny*nz);

	fwt_1D(d_tdata, lvly, ny, nz*nx);

	transpose(d_idata, d_tdata, ny, nz*nx);

	fwt_1D(d_idata, lvlz, nz, nx*ny);

	checkCuda( cudaMemcpy(data, d_idata, mem_size, cudaMemcpyDeviceToHost) );

	checkCuda( cudaFree(d_idata) );
	checkCuda( cudaFree(d_tdata) );
}

extern "C" void invComp(float *data, const unsigned nx, const unsigned ny, const unsigned nz, const unsigned lvlx, const unsigned lvly, const unsigned lvlz) {
	const int mem_size = nx*ny*nz*sizeof(float);

	float *d_idata, *d_tdata;
	checkCuda( cudaMalloc(&d_idata, mem_size) );
	checkCuda( cudaMalloc(&d_tdata, mem_size) );

	checkCuda( cudaMemcpy(d_idata, data, mem_size, cudaMemcpyHostToDevice) );
	//falta copiar los filtros
	iwt_1D(d_idata, lvlx, nz, nx*ny);

	transpose(d_tdata, d_idata, nz*nx, ny);

	fwt_1D(d_tdata, lvly, ny, nz*nx);

	transpose(d_idata, d_tdata, ny*nz, nx);

	fwt_1D(d_idata, lvlz, nx, ny*nz);

	checkCuda( cudaMemcpy(data, d_idata, mem_size, cudaMemcpyDeviceToHost) );

	checkCuda( cudaFree(d_idata) );
	checkCuda( cudaFree(d_tdata) );
}