#include <stdlib.h>
#include <stdio.h>
#include <cmath>

int main() {
  float f = -5.3f;
  f = std::abs(f);
  printf("abs(-5.3f) = %f\n", f);
  return 0;
}
