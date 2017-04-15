#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

//  g++ -std=c++11 -O2 thresh_sort.cc -o thresh_sort

int main(int argc, char **argv) {

  if (argc != 2) {
    printf("\n  thresh_sort <threshold value> < in.txt > out.txt\n\n");
    return 1;
  }

  float threshValue;
  if (1 != sscanf(argv[1], "%f", &threshValue) || threshValue < 0) {
    printf("Invalid threshold value: %s\n", argv[1]);
    return 1;
  }
  
  vector<float> values;
  float f;

  while (scanf("%f", &f) == 1) {
    f = fabsf(f);
    if (f >= threshValue) values.push_back(f);
  }

  sort(values.begin(), values.end());

  for (auto it = values.begin(); it != values.end(); it++) {
    printf("%f\n", *it);
  }

  return 0;
}

    
