#include <cstdio>
#include "cubelet_file.h"

void writeFile(const char *filename) {
  CubeletStreamWriter out;
  Cubelet cube;
  float *data;
  
  cube.datatype = Cubelet::CUBELET_FLOAT32;
  
  cube.data = data = (float*) malloc(sizeof(float) * 100 * 100 * 100);

  out.open(filename);

  cube.setSize(4, 5, 6);
  for (int i=0; i < 4*5*6; i++)
    data[i] = i * .1;
  if (!out.addCubelet(&cube)) return;

  cube.setSize(7, 8, 9);
  cube.x_offset = 1;
  for (int i=0; i < 7*8*9; i++)
    data[i] = i * .2;
  if (!out.addCubelet(&cube)) return;

  out.close();
}

void readFile(const char *filename) {
  CubeletStreamReader in;

  if (!in.open(filename)) return;

  Cubelet cube;

  while (true) {
    if (!in.next(&cube)) break;
    printf("Cubelet %dx%dx%d, offset %d,%d,%d\n",
           cube.width, cube.height, cube.depth,
           cube.x_offset, cube.y_offset, cube.z_offset);
  }

}

int main() {
  const char *filename = "test_cubelet_file.out";
  // const char *filename = NULL;
  writeFile(filename);
  readFile(filename);

  return 0;
}
