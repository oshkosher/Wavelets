#include <cstdio>
#include "dwt_cpu.h"

int main() {
  
  MirroredIterator iter;

  /*
    -5: 1
    -4: 0
    -3: 1
    -2: 2
    -1: 1
    0: 0
    1: 1
    2: 2
    3: 1
    4: 0
  */
  iter.init(3, -5);

  for (int i=-5; i < 10; i++) {
    printf("%d\n", iter++);
  }

  return 0;
}

  
