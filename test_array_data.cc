#include <cstdio>
#include "array_data.h"
#include "nixtimer.h"


void testTranspose2d() {
  Data3d src, dest;
  src.setSize(1000, 1000, 200);
  dest.setSize(1000, 1000, 200);
  src.allocate();
  dest.allocate();
  
  for (int z=0; z < src.depth(); z++)
    for (int y=0; y < src.height(); y++)
      for (int x=0; x < src.width(); x++)
        src.set(x, y, z, x + y*100 + z*10000);

  double startTime = NixTimer::time();
  src.transpose2dFast(dest);
  double elapsed = NixTimer::time() - startTime;

  printf("Transpose %dx%dx%d in %.3fms\n", src.width(), src.height(),
         src.depth(), elapsed*1000);

  for (int z=0; z < src.depth(); z++)
    for (int y=0; y < src.height(); y++)
      for (int x=0; x < src.width(); x++)
        if (src.get(x, y, z) != dest.get(y, x, z)) {
          printf("Mismatch at %d, %d, %d\n", x, y, z);
          return;
        }
  
}  


int main() {
  testTranspose2d();

}
