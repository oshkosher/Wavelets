#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <algorithm>
#include "thrust/random.h"
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "thrust/sort.h"
#include "thrust/execution_policy.h"
#include "cucheck.h"
#include "cuda_timer.h"

using namespace std;
using namespace std::chrono;

#define DEFAULT_SIZE (1024*1024)

class RandFloat {
  thrust::default_random_engine uniformRand;
  thrust::random::uniform_real_distribution<float> randGen;
public:
  RandFloat() : uniformRand(time(NULL)) {
  }

  float operator()() {
    return randGen(uniformRand);
  }
};

class CompareFloats {
public:
  __host__ __device__ bool operator()(float a, float b) {
    return a < b;
  }
};
    

void printHelp() {
  cout << "\n  test_sort [size]\n\n";
  exit(1);
}


int main2(int argc, char **argv) {

  if (argc > 2) printHelp();

  int size = DEFAULT_SIZE;
  if (argc > 1) {
    if (1 != sscanf(argv[1], "%d", &size) || size < 1) {
      cout << "Invalid size: \"" << argv[1] << "\"\n";
      return 1;
    }
  }

  // create a vector on the host
  thrust::host_vector<float> input(size);
  thrust::host_vector<float>::iterator iter;

  // fill it with random values between 0 and 1
  // default_random_engine rand;
  // uniform_real_distribution<float> randFloats(0, 1);
  // float f= randFloats(rand);

  RandFloat randGen;

  cout << "Input\n";
  high_resolution_clock::time_point startTime, endTime;
  startTime = high_resolution_clock::now();
  // *iter++ = randGen(uniformRand);
  thrust::generate(input.begin(), input.end(), randGen);
  endTime = high_resolution_clock::now();
  double sec = duration_cast<duration<double>>(endTime - startTime).count();
  cout << "Time to initialize: " << sec << " seconds\n";

  // make a copy of it and sort the copy
  thrust::host_vector<float> sorted = input;
  startTime = high_resolution_clock::now();
  thrust::sort(thrust::host, sorted.begin(), sorted.end());
  endTime = high_resolution_clock::now();
  sec = duration_cast<duration<double>>(endTime - startTime).count();
  cout << "Time to sort: " << sec << " seconds\n";
  cout.flush();

  // make a copy on the GPU
  thrust::device_vector<float> dev_data = input;
  float *dev_data_p = dev_data.data().get();
  cout << "Data copied to GPU\n";
  cudaEvent_t startEvent, endEvent;
  CUCHECK(cudaEventCreate(&startEvent));
  CUCHECK(cudaEventCreate(&endEvent));
  CUCHECK(cudaEventRecord(startEvent));

  thrust::sort(dev_data.begin(), dev_data.end());

  CUCHECK(cudaEventRecord(endEvent));
  CUCHECK(cudaEventSynchronize(endEvent));
  float elapsed;
  CUCHECK(cudaEventElapsedTime(&elapsed, startEvent, endEvent));
  cout << "Data sorted on GPU in " << elapsed << "ms\n";

  thrust::host_vector<float> sorted2 = dev_data;
  for (int i=0; i < sorted.size(); i++) {
    if (sorted[i] != sorted2[i]) {
      cout << "Error at [" << i << "]: mismatch " 
           << sorted[i] << " != " << sorted2[i] << endl;
      break;
    }
  }


  return 0;
}


int main(int argc, char **argv) {

  if (argc > 2) printHelp();

  int size = DEFAULT_SIZE;
  if (argc > 1) {
    if (1 != sscanf(argv[1], "%d", &size) || size < 1) {
      cout << "Invalid size: \"" << argv[1] << "\"\n";
      return 1;
    }
  }

  RandFloat randGen;

  // create data on the host
  float *input = new float[size];
  for (int i=0; i < size; i++) input[i] = randGen();

  // make a copy of it and sort the copy
  float *sorted = new float[size];
  memcpy(sorted, input, sizeof(float) * size);

  high_resolution_clock::time_point startTime, endTime;
  startTime = high_resolution_clock::now();

  std::sort(sorted, sorted+size);
  // thrust::sort(sorted, sorted+size);

  endTime = high_resolution_clock::now();
  double sec = duration_cast<duration<double>>(endTime - startTime).count();
  cout << "Time to sort on the CPU: " << sec << " seconds\n";

  // send the data to the GPU
  CudaTimer timer;
  float *dev_data;
  CUCHECK(cudaMalloc(&dev_data, sizeof(float) * size));
  timer.start();
  CUCHECK(cudaMemcpy(dev_data, input, sizeof(float) * size,
                     cudaMemcpyHostToDevice));
  timer.end();

  // sort it on the GPU
  thrust::device_ptr<float> dev_start_p(dev_data), dev_end_p(dev_data + size);
  timer.start();

  // Even if the comparison function object is just "a < b", it will
  // be much slower than a native comparison. On my test machine, it
  // slows from 12ms to 43ms.
  // CompareFloats compareFloats;
  thrust::sort(dev_start_p, dev_end_p /* ,compareFloats */ );

  timer.end();

  // copy it back from the GPU
  float *gpu_result = new float[size];
  timer.start();
  CUCHECK(cudaMemcpy(gpu_result, dev_data, sizeof(float) * size,
                     cudaMemcpyDeviceToHost));
  timer.end();
  CUCHECK(cudaStreamSynchronize(0));
  cout << "Copy to GPU: " << timer.getEventTime(0) << " ms\n";
  cout << "Sort on GPU: " << timer.getEventTime(1) << " ms\n";
  cout << "Copy from GPU: " << timer.getEventTime(2) << " ms\n";

  for (int i=0; i < size; i++) {
    if (sorted[i] != gpu_result[i]) {
      cout << "Error at [" << i << "]: mismatch " 
           << sorted[i] << " != " << gpu_result[i] << endl;
      break;
    }
  }

  return 0;
}

