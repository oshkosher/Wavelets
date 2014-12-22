#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/transform.h>

struct TenX {
  __host__ __device__ int operator() (int x) const {
    return x*10;
  }
} myFunctor;

void hostVectors() {
  thrust::host_vector<int> vec1(4), vec2(4);

  printf("Host\n");
  for (size_t i=0; i < vec1.size(); i++) vec1[i] = i;
  thrust::transform(vec1.begin(), vec1.end(), vec2.begin(), myFunctor);

  for (size_t i=0; i < vec1.size(); i++) {
    printf("%d\t%d\n", vec1[i], vec2[i]);
  }
}

void deviceVectors() {
  thrust::host_vector<int> vec1(4), vec2(4);
  thrust::device_vector<int> vec1_dev(4), vec2_dev(4);

  printf("Device\n");
  for (size_t i=0; i < vec1.size(); i++) vec1[i] = i;
  vec1_dev = vec1;
  thrust::transform(vec1_dev.begin(), vec1_dev.end(), vec2_dev.begin(),
                    myFunctor);
  vec2 = vec2_dev;

  for (size_t i=0; i < vec1.size(); i++) {
    printf("%d\t%d\n", vec1[i], vec2[i]);
  }
}

void deviceVectorDevPointers() {
  thrust::host_vector<int> vec1(4), vec2(4);
  thrust::device_vector<int> vec1_dev(4), vec2_dev(4);
  thrust::device_ptr<int> vec1_dev_start(thrust::raw_pointer_cast(&vec1_dev[0]));
  thrust::device_ptr<int> vec2_dev_start(thrust::raw_pointer_cast(&vec2_dev[0]));

  printf("Device device pointers\n");
  for (size_t i=0; i < vec1.size(); i++) vec1[i] = i;
  vec1_dev = vec1;
  thrust::transform(thrust::device, vec1_dev_start, vec1_dev_start + 4,
                    vec2_dev_start, myFunctor);
  vec2 = vec2_dev;

  for (size_t i=0; i < vec1.size(); i++) {
    printf("%d\t%d\n", vec1[i], vec2[i]);
  }
}

void deviceVectorPointers() {
  thrust::host_vector<int> vec1(4), vec2(4);
  thrust::device_vector<int> vec1_dev(4), vec2_dev(4);
  int *vec1_dev_start = thrust::raw_pointer_cast(&vec1_dev[0]);
  int *vec2_dev_start = thrust::raw_pointer_cast(&vec2_dev[0]);

  printf("Device pointers\n");
  for (size_t i=0; i < vec1.size(); i++) vec1[i] = i;
  vec1_dev = vec1;
  thrust::transform(thrust::device, vec1_dev_start, vec1_dev_start + 4,
                    vec2_dev_start, myFunctor);
  vec2 = vec2_dev;

  for (size_t i=0; i < vec1.size(); i++) {
    printf("%d\t%d\n", vec1[i], vec2[i]);
  }
}

void deviceVecToPointerToDevPointer() {
  thrust::host_vector<int> vec1(4), vec2(4);
  thrust::device_vector<int> vec1_dev(4), vec2_dev(4);
  int *vec1_dev_start = thrust::raw_pointer_cast(&vec1_dev[0]);
  int *vec2_dev_start = thrust::raw_pointer_cast(&vec2_dev[0]);
  thrust::device_ptr<int> vec1_dev_ptr(vec1_dev_start);
  thrust::device_ptr<int> vec2_dev_ptr(vec2_dev_start);

  printf("Device vec->ptr->devptr\n");
  for (size_t i=0; i < vec1.size(); i++) vec1[i] = i;
  vec1_dev = vec1;
  thrust::transform(thrust::device, vec1_dev_ptr, vec1_dev_ptr + 4,
                    vec2_dev_ptr, myFunctor);
  vec2 = vec2_dev;

  for (size_t i=0; i < vec1.size(); i++) {
    printf("%d\t%d\n", vec1[i], vec2[i]);
  }
}
  

int main() {

  hostVectors();
  deviceVectors();
  deviceVectorDevPointers();
  deviceVectorPointers();
  deviceVecToPointerToDevPointer();

  return 0;
}
