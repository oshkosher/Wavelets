#include "histogram_gpu.h"
#include "cucheck.h"
#include "test_compress_common.h"
#include "cuda_timer.h"

#define HISTOGRAM_CACHE_SIZE 256

// #define DO_DETAILED_CHECKS

/**
   Histogram computation.
   
   Scan through the input data, incrementing the appropriate entry
   in a shared memory array for each input value. At the end, add
   the data from the shared memory array from each thread block to
   the global result array.

   There may be more bins than the shared array can hold. If so, do
   multiple passes over the data, each time ignoring values outside
   the current range of values being summarized. For example, if
   HISTOGRAM_CACHE_SIZE is 100, then on the first pass ignore anything
   that is outside the range 0..99. On the second pass ignore anything
   outside the range 100..199, and so on.

   The wavelet compression algorithm generates many elements with the
   value 0, so there will be much contention for that bin. The caller
   can pass in the value of that bin, and the counters for that bin
   will be stored in a local variable to eliminate contention until the
   final summation.

   Creating a special case for the zero bin is a huge win.
   In one test, it reduces the runtime from 848ms to 48.3 ms.
*/
__global__ void computeFrequenciesKernel
(int *freqCounts, int binCount, const int *data, int count, int zeroBin) {

  // compute frequencies for HISTOGRAM_CACHE_SIZE bins at a time.
  // If binCount > HISTOGRAM_CACHE_SIZE, do multiple passes over the 
  // data.
  __shared__ int sharedFreq[HISTOGRAM_CACHE_SIZE];
  __shared__ int sharedZeroCount;

  // clear the shared memory
  if (threadIdx.x == 0) sharedZeroCount = 0;
  for (int i=threadIdx.x; i < HISTOGRAM_CACHE_SIZE; i += blockDim.x) {
    sharedFreq[i] = 0;
  }

  __syncthreads();

  int zeroCount = 0;

  // There may be more bins than HISTOGRAM_CACHE_SIZE.
  // Do multiple passes over the data if necessary.
  for (int startBin = 0;
       startBin < binCount;
       startBin += HISTOGRAM_CACHE_SIZE) {

    for (int i = threadIdx.x + blockIdx.x * blockDim.x;
         i < count;
         i += blockDim.x * gridDim.x) {

      int b = data[i];

      // The zero bin will have many updates and experience high contention.
      // Handle it with a local variable.
      if (b == zeroBin) {

        // only process the zero bin when in the first pass over the data,
        // regardless of which pass range its value fits in.
        if (startBin == 0) {
          zeroCount++;
        }

      } else {
        unsigned bin = (unsigned) (b - startBin);
        if (bin < HISTOGRAM_CACHE_SIZE) {
          atomicAdd(&sharedFreq[bin], 1);
        }

      }
    }

    __syncthreads();

    // write shared memory cache to global memory
    for (int i=threadIdx.x;
         i < HISTOGRAM_CACHE_SIZE && startBin+i < binCount;
         i += blockDim.x) {
      if (sharedFreq[i] > 0) {
        atomicAdd(&freqCounts[startBin + i], sharedFreq[i]);
        sharedFreq[i] = 0;
      }
    }

    // sync so the zeros written to shared memory to clear it are seen
    __syncthreads();
  }

  // gather the results for zeroCount into shared memory
  atomicAdd(&sharedZeroCount, zeroCount);
  __syncthreads();

  // gather the results for zeroCount into global memory
  if (threadIdx.x == 0) {
    atomicAdd(&freqCounts[zeroBin], sharedZeroCount);
  }

}


void computeFrequenciesGPU(int *&freqCounts_dev, int binCount,
                           const int *data_dev, int count, int zeroBin) {
                        
  CudaTimer timer("Frequency count");

  // 16 thread blocks worked well with the cards I tried.
  // With the laptop GPU I tried, 512 threads/block gave the best performance,
  // but 1024 threads/block worked best with the desktop GPU I tried.
  int blockSize, gridSize = 16;

  int gpuId;
  cudaDeviceProp prop;
  CUCHECK(cudaGetDevice(&gpuId));
  CUCHECK(cudaGetDeviceProperties(&prop, gpuId));
  blockSize = prop.multiProcessorCount <= 2 ? 512 : 1024;

  timer.start();
  CUCHECK(cudaMalloc((void**)&freqCounts_dev, binCount * sizeof(int)));
  CUCHECK(cudaMemset(freqCounts_dev, 0, binCount * sizeof(int)));

  computeFrequenciesKernel<<<gridSize, blockSize>>>
    (freqCounts_dev, binCount, data_dev, count, zeroBin);

  timer.end();
  if (!QUIET) {
    timer.sync();
    timer.print();
    fflush(stdout);
  }

#ifdef DO_DETAILED_CHECKS
  // copy the results to the CPU, compute them on the CPU and compare
  int *qdata = new int[count];
  CUCHECK(cudaMemcpy(qdata, data_dev, sizeof(int) * count,
                     cudaMemcpyDeviceToHost));
  vector<int> counts;
  counts.resize(binCount);
  for (int i=0; i < count; i++) counts[qdata[i]]++;
  delete[] qdata;

  int *freqCounts = new int[binCount];
  CUCHECK(cudaMemcpy(freqCounts, freqCounts_dev, sizeof(int) * binCount,
                     cudaMemcpyDeviceToHost));
  for (int i=0; i < binCount; i++) {
    if (counts[i] != freqCounts[i]) {
      fprintf(stderr, "Mismatch for %d: %d cpu, %d gpu\n", i, counts[i], freqCounts[i]);
    }
  }
  delete[] freqCounts;
#endif
         
}
