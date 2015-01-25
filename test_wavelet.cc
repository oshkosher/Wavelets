#include "wavelet.h"
#include "nixtimer.h"

using namespace scu_wavelet;

void randomize(CubeFloat &c, const int3 &size) {
  c.size = size;
  c.allocate();
  for (int i=0; i < c.count(); i++)
    c.set(i, 0, 0, rand() % 100);
}
  

bool verifyTranspose2D(const CubeFloat &before, const CubeFloat &after) {

  // before.print("before", 0);
  // after.print("after", 0);

  if (before.size.x != after.size.y ||
      before.size.y != after.size.x ||
      before.size.z != after.size.z) {
    printf("size mismatch\n");
    return false;
  }

  for (int z=0; z < before.size.z; z++) {
    for (int y=0; y < before.size.y; y++) {
      for (int x=0; x < before.size.x; x++) {
        if (after.get(y, x, z) != before.get(x, y, z)) {
          printf("Mismatch before(%d,%d,%d)=%f, after(%d,%d,%d)=%f\n",
                 x, y, z, before.get(x, y, z),
                 y, x, z, after.get(y, x, z));
          return false;
        }
      }
    }
  }
  return true;
}


bool verifyTranspose3D(const CubeFloat &before, const CubeFloat &after) {

  // before.print("before", 0);
  // after.print("after", 0);

  if (before.size.x != after.size.z ||
      before.size.y != after.size.x ||
      before.size.z != after.size.y) {
    printf("size mismatch\n");
    return false;
  }

  for (int z=0; z < before.size.z; z++) {
    for (int y=0; y < before.size.y; y++) {
      for (int x=0; x < before.size.x; x++) {
        if (after.get(y, z, x) != before.get(x, y, z)) {
          printf("Mismatch before(%d,%d,%d)=%f, after(%d,%d,%d)=%f\n",
                 x, y, z, before.get(x, y, z),
                 y, z, x, after.get(y, x, z));
          return false;
        }
      }
    }
  }
  return true;
}


int main_step_sizes() {
  // srand(time(NULL));

  CubeFloat orig, cube;
  orig.size = int3(512, 512, 512);
  orig.allocate();
  for (int i=0; i < orig.count(); i++)
    orig.set(i, 0, 0, rand() % 100);

  cube.size = orig.size;
  cube.allocate();


  printf("      ");
  for (int xstep = 1; xstep <= 16; xstep++)
    printf("%6d ", xstep);
  putchar('\n');

  for (int ystep = 17; ystep <= 32; ystep++) {
    printf("%5d ", ystep);
    for (int xstep = 1; xstep <= 16; xstep++) {

      cube.copyFrom(orig);

      // orig.print("original");

      double startTime = NixTimer::time();

      // cube.transpose2d();
      // cube.transpose2dFast(xstep, ystep);

      double elapsedMs = (NixTimer::time() - startTime) * 1000;

      printf("%6.1f ", elapsedMs); fflush(stdout);

      // printf("transpose %dx%dx%d step %d: %.3f ms\n", cube.size.x, cube.size.y, cube.size.z, step, elapsedMs);

      // if (!verifyTranspose2D(orig, cube)) return 1;
    }

    putchar('\n'); fflush(stdout);
  }

  return 0;
}


int main_2d_perf() {
  CubeFloat orig, cube;
  randomize(orig, int3(512, 512, 512));

  cube.size = orig.size;
  cube.allocate();
  float totalTime = 0;

  // orig.print("Original");

  for (int i=0; i < 5; i++) {
    cube.size = orig.size;
    cube.totalSize = orig.totalSize;
    cube.copyFrom(orig);
    
    double startTime = NixTimer::time();
    
    cube.transpose2d();

    // cube.transpose2dFast();
    
    double elapsedMs = (NixTimer::time() - startTime) * 1000;
    totalTime += elapsedMs;
    
    printf("%.1f ", elapsedMs);
    fflush(stdout);
    verifyTranspose2D(orig, cube);
  }
  printf("   total = %.1f\n", totalTime);

  return 0;
}


int main_3d_perf() {
  CubeFloat orig, cube;
  randomize(orig, int3(512, 512, 512));

  cube.size = orig.size;
  cube.allocate();
  cube.copyFrom(orig);
  float totalTime = 0;

  for (int i=0; i < 5; i++) {
    cube.size = orig.size;
    cube.totalSize = orig.totalSize;
    cube.copyFrom(orig);
    
    // orig.print("Original");

    double startTime = NixTimer::time();
    
    cube.transpose3dBack();

    // cube.transpose3dBackSlow();
    
    double elapsedMs = (NixTimer::time() - startTime) * 1000;
    totalTime += elapsedMs;

    // cube.print("Transposed");
    
    printf("%.1f ", elapsedMs);
    fflush(stdout);
    verifyTranspose3D(cube, orig);
  }

  printf("   total = %.1f\n", totalTime);

  return 0;
}


int main() {
  main_3d_perf();
}
