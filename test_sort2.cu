#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <cstdlib>

#include <iostream>
#include <iterator>

int main(void)
{
  // generate random data on the host
  thrust::host_vector<int> h_vec(20);
  thrust::generate(h_vec.begin(), h_vec.end(), rand);
  std::cerr << "input..." << std::endl;
  std::copy( h_vec.begin(), h_vec.end(), std::ostream_iterator<int>(std::cerr, "\n") );
  std::cerr << "" << std::endl;

  // transfer to device and sort
  thrust::device_vector<int> d_vec = h_vec;
  std::cerr << "gpu[0] = " << d_vec[0] << ", gpu[1] = " << d_vec[1] << '\n';
  thrust::sort(d_vec.begin(), d_vec.end());
  std::cerr << "sort...\n";
  std::cerr << "gpu[0] = " << d_vec[0] << ", gpu[1] = " << d_vec[1] << '\n';

  // show result
  thrust::host_vector<int> h_result = d_vec;
  std::cerr << "output..." << std::endl;
  std::copy( h_result.begin(), h_result.end(), std::ostream_iterator<int>(std::cerr, "\n") );
  std::cerr << "" << std::endl;

  std::cerr << "third item in sorted data:" << d_vec[2] << std::endl;

  return 0;
}
