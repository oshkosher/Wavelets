#include "cuda.h"

#include <cstdio>

static float *h_points;
static float *d_points;
static double *h_pointsd;
static double *d_pointsd;
static unsigned int *d_groups;
static unsigned int d_psize;
static float d_pmax;
static float d_pmin;
static double d_pmaxd;
static double d_pmind;

__global__ void groupKernel(float *points, unsigned int psize, float *codebook, unsigned int csize, unsigned int *groups) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < psize) {
		float d = abs(points[idx] - codebook[0]);
		float min = d;
		unsigned int g = 0;
		for(int i = 1; i < csize; i++) {
			d = abs(points[idx] - codebook[i]);
			g += (d < min) * (i - g);
			min -= (d < min) * (min - d);
		}
		groups[idx] = g;
	}
}

__global__ void groupKernel(double *points, unsigned int psize, double *codebook, unsigned int csize, unsigned int *groups) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < psize) {
		double d = abs(points[idx] - codebook[0]);
		double min = d;
		unsigned int g = 0;
		for(int i = 1; i < csize; i++) {
			d = abs(points[idx] - codebook[i]);
			g += (d < min) * (i - g);
			min -= (d < min) * (min - d);
		}
		groups[idx] = g;
	}
}

void setupLloyds(float *points, unsigned int psize, float pmax, float pmin) {
	d_psize = psize;
	d_pmax = pmax;
	d_pmin = pmin;

	h_points = points;
  cudaMalloc((void**)&d_points, sizeof(float)*psize);
  cudaMalloc((void**)&d_groups, sizeof(unsigned int)*psize);

	cudaMemcpy(d_points, points, sizeof(float) * psize, cudaMemcpyHostToDevice);
}

void setupLloyds(double *points, unsigned int psize, double pmax, double pmin) {
	d_psize = psize;
	d_pmaxd = pmax;
	d_pmind = pmin;

	h_pointsd = points;
  cudaMalloc((void**)&d_pointsd, sizeof(double)*psize);
  cudaMalloc((void**)&d_groups, sizeof(unsigned int)*psize);

	cudaMemcpy(d_pointsd, points, sizeof(double) * psize, cudaMemcpyHostToDevice);
}

void lloyd(float *codebook, unsigned int csize, float stop_criteria, float *table, float &dist, float &reldist) {

	float *d_codebook;
	dist = 0;
	reldist = 0;

	unsigned int *groups = (unsigned int*)malloc(sizeof(unsigned int)*d_psize);

	// Calculate initial table
	//can be done in gpu if csize is big enough
	for(int i = 0; i < csize-1; i++) {
		table[i] = (codebook[i] + codebook[i+1]) / 2;
	}

  cudaMalloc((void**)&d_codebook, sizeof(float)*csize);
	cudaMemcpy(d_codebook, codebook, sizeof(float)*csize, cudaMemcpyHostToDevice);

	// Assign each point its codebook group
	unsigned int threads = 256;
	unsigned int blocks = (d_psize - 1) / threads + 1;
	groupKernel<<<blocks, threads>>>(d_points, d_psize, d_codebook, csize, d_groups);

	cudaMemcpy(groups, d_groups, sizeof(unsigned int) * d_psize, cudaMemcpyDeviceToHost);
	
	unsigned int *incode = (unsigned int*)calloc(csize, sizeof(unsigned int));
	float *meancode = (float*)calloc(csize, sizeof(float));

	// Calculate the mean of each codebook group and distortion
	for(int i = 0; i < d_psize; i++) {
		incode[groups[i]]++;
		meancode[groups[i]] += h_points[i];
		dist += codebook[groups[i]] - h_points[i];
	}
	dist /= d_psize;

	reldist = abs(dist);

	while(reldist > stop_criteria) {

		// Update codebook
		//can be done in gpu if csize is big enough
		for(int i = 0; i < csize; i++) {
			if(incode[i] != 0) {
				codebook[i] = meancode[i] / incode[i];
				incode[i] = 0;
			} else if (i == 0) {
				codebook[i] = (table[0] + d_pmin) / 2;
			} else if (i == csize-1) {
				codebook[i] = (table[i-1] + d_pmax) / 2;
			} else {
				codebook[i] = (table[i-1] + table[i]) / 2;
			}
			meancode[i] = 0;
		}

		cudaMemcpy(d_codebook, codebook, sizeof(float)*csize, cudaMemcpyHostToDevice);

		// How to do this in CUDA efficiently?
		for(int i = 0; i < d_psize; i++) {
			for(int j = 0; j < csize-1; j++) {
				if((h_points[i] > codebook[j])
				&& (h_points[i] < codebook[j+1])) {
					incode[j]++;
					meancode[j] += h_points[i];
					break;
				}
			}
		}

		// Update table
		//can be done in gpu if csize is big enough
		for(int i = 0; i < csize-1; i++) {
			if(incode[i] != 0) {
				table[i] = meancode[i] / incode[i];
				incode[i] = 0;
			} else {
				table[i] = (codebook[i] + codebook[i+1]) / 2;
			}
			meancode[i] = 0;
		}

		// Assign each point its codebook group
		groupKernel<<<blocks, threads>>>(d_points, d_psize, d_codebook, csize, d_groups);

		cudaMemcpy(groups, d_groups, sizeof(unsigned int) * d_psize, cudaMemcpyDeviceToHost);

		reldist = dist;
		dist = 0;

		// Calculate the mean of each codebook group and distortion
		for(int i = 0; i < d_psize; i++) {
			incode[groups[i]]++;
			meancode[groups[i]] += h_points[i];
			dist += codebook[groups[i]] - h_points[i];
		}
		dist /= d_psize;

		reldist = abs(reldist - dist);
	}  // END WHILE

	cudaFree(d_codebook);
	free(groups);
	free(incode);
	free(meancode);
}

void lloyd(double *codebook, unsigned int csize, double stop_criteria, double *table, double &dist, double &reldist) {

	double *d_codebook;
	dist = 0;
	reldist = 0;

	unsigned int *groups = (unsigned int*)malloc(sizeof(unsigned int)*d_psize);

	// Calculate initial table
	//can be done in gpu if csize is big enough
	for(int i = 0; i < csize-1; i++) {
		table[i] = (codebook[i] + codebook[i+1]) / 2;
	}

  cudaMalloc((void**)&d_codebook, sizeof(double)*csize);
	cudaMemcpy(d_codebook, codebook, sizeof(double)*csize, cudaMemcpyHostToDevice);

	// Assign each point its codebook group
	unsigned int threads = 256;
	unsigned int blocks = (d_psize - 1) / threads + 1;
	groupKernel<<<blocks, threads>>>(d_pointsd, d_psize, d_codebook, csize, d_groups);

	cudaMemcpy(groups, d_groups, sizeof(unsigned int) * d_psize, cudaMemcpyDeviceToHost);
	
	unsigned int *incode = (unsigned int*)calloc(csize, sizeof(unsigned int));
	double *meancode = (double*)calloc(csize, sizeof(double));

	// Calculate the mean of each codebook group and distortion
	for(int i = 0; i < d_psize; i++) {
		incode[groups[i]]++;
		meancode[groups[i]] += h_pointsd[i];
		dist += codebook[groups[i]] - h_pointsd[i];
	}
	dist /= d_psize;

	reldist = abs(dist);

	while(reldist > stop_criteria) {

		// Update codebook
		//can be done in gpu if csize is big enough
		for(int i = 0; i < csize; i++) {
			if(incode[i] != 0) {
				codebook[i] = meancode[i] / incode[i];
				incode[i] = 0;
			} else if (i == 0) {
				codebook[i] = (table[0] + d_pmind) / 2;
			} else if (i == csize-1) {
				codebook[i] = (table[i-1] + d_pmaxd) / 2;
			} else {
				codebook[i] = (table[i-1] + table[i]) / 2;
			}
			meancode[i] = 0;
		}

		cudaMemcpy(d_codebook, codebook, sizeof(double)*csize, cudaMemcpyHostToDevice);

		// How to do this in CUDA efficiently?
		for(int i = 0; i < d_psize; i++) {
			for(int j = 0; j < csize-1; j++) {
				if((h_pointsd[i] > codebook[j])
				&& (h_pointsd[i] < codebook[j+1])) {
					incode[j]++;
					meancode[j] += h_pointsd[i];
					break;
				}
			}
		}

		// Update table
		//can be done in gpu if csize is big enough
		for(int i = 0; i < csize-1; i++) {
			if(incode[i] != 0) {
				table[i] = meancode[i] / incode[i];
				incode[i] = 0;
			} else {
				table[i] = (codebook[i] + codebook[i+1]) / 2;
			}
			meancode[i] = 0;
		}

		// Assign each point its codebook group
		groupKernel<<<blocks, threads>>>(d_pointsd, d_psize, d_codebook, csize, d_groups);

		cudaMemcpy(groups, d_groups, sizeof(unsigned int) * d_psize, cudaMemcpyDeviceToHost);

		reldist = dist;
		dist = 0;

		// Calculate the mean of each codebook group and distortion
		for(int i = 0; i < d_psize; i++) {
			incode[groups[i]]++;
			meancode[groups[i]] += h_pointsd[i];
			dist += codebook[groups[i]] - h_pointsd[i];
		}
		dist /= d_psize;

		reldist = abs(reldist - dist);
	}  // END WHILE

	cudaFree(d_codebook);
	free(groups);
	free(incode);
	free(meancode);
}


void finalize() {
	cudaFree(d_points);
	cudaFree(d_groups);
}

void finalized() {
	cudaFree(d_pointsd);
	cudaFree(d_groups);
}

