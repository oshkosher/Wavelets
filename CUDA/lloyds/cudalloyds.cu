/*****************************************************************************/
/* CUDA implementation of Lloyd's quantization algorithm                     */
/* lloyd(points, psize, codebook, csize, stop_criteria =  10e-7)             */
/*		                                                                     */
/*  - IN     - points:        values to be quantized                         */
/*  - IN     - psize:         size of the points array                       */
/*  - IN/OUT - codebook:      initial codebook and final codebook            */
/*  - IN     - csize:         size of the codebook array                     */
/*  - IN     - stop_criteria: typically 10e-7                                */
/*  - IN     - points_are_on_gpu: if true, point data is already on the GPU, */
/*                            and 'points' is a device address.              */
/*                                                                           */
/*  Doesn't return the partition because it can be easily calculated         */
/*  as the mid-point between codebooks.                                      */
/*****************************************************************************/
__device__ static float atomicMax(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_i, assumed,
            __float_as_int(::fmaxf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

__device__ static float atomicMin(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_i, assumed,
            __float_as_int(::fminf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

__global__ void groupKernelMaxMin(const float *points, unsigned int psize, float *codebook, unsigned int csize, unsigned int *groups, unsigned int *counts, float *sum, float *dist, float *max, float *min) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < psize) {
		float p = points[idx];
		float d = abs(p - codebook[0]);
		float m = d;
		unsigned int g = 0;
		for(int i = 1; i < csize; i++) {
			d = abs(p - codebook[i]);
			bool isminus = (d < m);
			g += isminus * (i - g);
			m -= isminus * (m - d);
		}
		groups[idx] = g;
		atomicAdd(&(counts[g]),1);
		atomicAdd(&(sum[g]),p);
		atomicAdd(dist,codebook[g]-p);
		atomicMax(max,p);
		atomicMin(min,p);
	}
}

__global__ void groupKernel(const float *points, unsigned int psize, float *codebook, unsigned int csize, unsigned int *groups, unsigned int *counts, float *sum, float *dist) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < psize) {
		float p = points[idx];
		float d = abs(p - codebook[0]);
		float min = d;
		unsigned int g = 0;
		for(int i = 1; i < csize; i++) {
			d = abs(p - codebook[i]);
			bool isminus = (d < min);
			g += isminus * (i - g);
			min -= isminus * (min - d);
		}
		groups[idx] = g;
		atomicAdd(&(counts[g]),1);
		atomicAdd(&(sum[g]),p);
		atomicAdd(dist,fabs(codebook[g]-p));
	}
}


void cudaLloyd(const float *points, unsigned int psize, float *codebook, unsigned int csize, float stop_criteria, bool points_are_on_gpu, int *iterations) {
	const float *d_points;
	unsigned *d_groups;
	float *d_codebook;
	float *partition = (float*)malloc((csize-1)*sizeof(float));;
	unsigned *counts = (unsigned*)calloc(csize,sizeof(unsigned));
	unsigned *d_counts;
	float *sum = (float*)calloc(csize,sizeof(float));
	float *d_sum;
	float dist;
	float reldist;
	float *d_dist;
	float *d_max,max;
	float *d_min,min;

        // check if the point data is already on the GPU
        if (points_are_on_gpu) {
          d_points = points;
        } else {
          cudaMalloc((void**)&d_points, sizeof(float)*psize);
          cudaMemcpy((void*)d_points, points, sizeof(float) * psize, cudaMemcpyHostToDevice);
        }

	cudaMalloc((void**)&d_groups, sizeof(unsigned)*psize);
	cudaMalloc((void**)&d_codebook, sizeof(float)*csize);
	cudaMalloc((void**)&d_counts, sizeof(unsigned)*csize);
	cudaMalloc((void**)&d_sum, sizeof(float)*csize);
	cudaMalloc((void**)&d_dist, sizeof(float));
	cudaMalloc((void**)&d_max, sizeof(float));
	cudaMalloc((void**)&d_min, sizeof(float));

	cudaMemcpy(d_codebook, codebook, sizeof(float)*csize, cudaMemcpyHostToDevice);
	cudaMemset(d_counts,0,sizeof(unsigned)*csize);
	cudaMemset(d_sum,0,sizeof(float)*csize);
	cudaMemset(d_dist,0,sizeof(float));
	cudaMemcpy(d_max, d_points, sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_min, d_points, sizeof(float), cudaMemcpyDeviceToDevice);

        if (iterations) *iterations = 0;

	// Initial Table
	for(int i = 0; i < csize-1; i++) {
		partition[i] = (codebook[i] + codebook[i+1]) / 2;
	}

	// Assign each point its codebook group
	unsigned int threads = 256;
	unsigned int blocks = (psize - 1) / threads + 1;
	groupKernelMaxMin<<<blocks, threads>>>(d_points, psize, d_codebook, csize, d_groups, d_counts, d_sum, d_dist,d_max, d_min);

	cudaMemcpy(&max, d_max, sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(&min, d_min, sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(&dist, d_dist, sizeof(float), cudaMemcpyDeviceToHost);
	reldist = abs(dist);

	while(reldist > stop_criteria) {
		if (iterations) ++*iterations;
		cudaMemcpy(counts, d_counts, sizeof(unsigned)*csize, cudaMemcpyDeviceToHost);
		cudaMemcpy(sum, d_sum, sizeof(float)*csize, cudaMemcpyDeviceToHost);

		// Update codebook
		//can be done in gpu if csize is big enough
		for(int i = 0; i < csize; i++) {
			if(counts[i] != 0) {
				codebook[i] = sum[i] / counts[i];
				counts[i] = 0;
			} else if (i == 0) {
				codebook[i] = (partition[0] + min) / 2;
			} else if (i == csize-1) {
				codebook[i] = (partition[i-1] + max) / 2;
			} else {
				codebook[i] = (partition[i-1] + partition[i]) / 2;
			}
			sum[i] = 0;
		}
		cudaMemcpy(d_codebook, codebook, sizeof(float)*csize, cudaMemcpyHostToDevice);
		cudaMemset(d_counts,0,sizeof(unsigned)*csize);
		cudaMemset(d_sum,0,sizeof(float)*csize);

		// Update Table
		for(int i = 0; i < csize-1; i++) {
			partition[i] = (codebook[i] + codebook[i+1]) / 2;
		}

		reldist = dist;
		cudaMemset(d_dist,0,sizeof(float));

		// Assign each point its codebook group
		groupKernel<<<blocks, threads>>>(d_points, psize, d_codebook, csize, d_groups, d_counts, d_sum, d_dist);

		//cudaMemcpy(groups, d_groups, sizeof(unsigned int) * psize, cudaMemcpyDeviceToHost);

		cudaMemcpy(&dist, d_dist, sizeof(float), cudaMemcpyDeviceToHost);
		dist /= psize;

		reldist = abs(reldist - dist);
	}  // END WHILE

	free(partition);
	free(counts);
	free(sum);
	if (!points_are_on_gpu) cudaFree((void*)d_points);
	cudaFree(d_groups);
	cudaFree(d_codebook);
	cudaFree(d_counts);
	cudaFree(d_sum);
	cudaFree(d_dist);
	cudaFree(d_max);
	cudaFree(d_min);
}
