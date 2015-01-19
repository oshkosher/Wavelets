#include <cstdio>
#include "cubelet_file.h"

using namespace scu_wavelet;

void writeFile(const char *filename) {
  CubeletStreamWriter out;
  CubeFloat cube;
  // float *data;
  
  // cube.datatype = Cubelet::CUBELET_FLOAT32;
  
  cube.size = int3(100,100,100);
  cube.allocate();

  out.open(filename);

  cube.size = int3(4, 5, 6);
  cube.parentOffset = int3(100, 101, 102);
  for (int i=0; i < 4*5*6; i++)
    cube.set(i, 0, 0, i * .25);
  if (!out.addCubelet(&cube)) return;

  cube.size = int3(7, 8, 9);
  cube.parentOffset = int3(200, 201, 202);
  for (int i=0; i < 7*8*9; i++)
    cube.set(i, 0, 0, i * .5);
  if (!out.addCubelet(&cube)) return;

  out.close();
}

void readFile(const char *filename) {
  CubeletStreamReader in;

  if (!in.open(filename)) return;

  Cube cube;

  while (true) {
    if (!in.next(&cube)) break;
    printf("Cubelet %dx%dx%d, offset %d,%d,%d\n",
           cube.width(), cube.height(), cube.depth(),
           cube.parentOffset.x, cube.parentOffset.y, cube.parentOffset.z);

    if (cube.isWaveletCompressed) continue;

    in.getCubeData(&cube);

    if (cube.datatype == WAVELET_DATA_FLOAT32) {
      CubeFloat *data = (CubeFloat*) &cube;
      data->print();
    } else if (cube.datatype == WAVELET_DATA_UINT8) {
      CubeByte *data = (CubeByte*) &cube;
      data->print();
    } else if (cube.datatype == WAVELET_DATA_INT32) {
      CubeInt *data = (CubeInt*) &cube;
      data->print();
    }

  }

}

int main(int argc, char **argv) {
  const char *filename = argv[1];
  if (filename && !strcmp(filename, "-")) filename = NULL;

  // writeFile(filename);
  readFile(filename);

  return 0;
}
