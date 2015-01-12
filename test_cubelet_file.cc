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
  cube.setOffset(100, 101, 102);
  for (int i=0; i < 4*5*6; i++)
    data[i] = i * .25;
  if (!out.addCubelet(&cube)) return;

  cube.setSize(7, 8, 9);
  cube.setOffset(200, 201, 202);
  for (int i=0; i < 7*8*9; i++)
    data[i] = i * .5;
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
           cube.xOffset, cube.yOffset, cube.zOffset);

    /*
    if (cube.datatype == Cubelet::CUBELET_FLOAT32) {

      float *data = (float*) in.getData();

      for (unsigned y=0; y < cube.height; y++) {
        for (unsigned x=0; x < cube.width; x++) {
          printf("%6.2f", data[y*cube.width + x]);
        }
        printf("\n");
      }

    } else if (cube.datatype == Cubelet::CUBELET_UINT8) {

      unsigned char *data = (unsigned char*) in.getData();

      for (unsigned y=0; y < cube.height; y++) {
        for (unsigned x=0; x < cube.width; x++) {
          printf("%4d", data[y*cube.width + x]);
        }
        printf("\n");
      }

    }

    printf("\n");
    */

  }

}

int main(int argc, char **argv) {
  const char *filename = argv[1];
  if (filename && !strcmp(filename, "-")) filename = NULL;

  // writeFile(filename);
  readFile(filename);

  return 0;
}
