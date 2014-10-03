#include <cstdio>
#include <cassert>
#include "dwt_cpu.h"
#include "data_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dwt_cpu.h"
#include "data_io.h"
#include "thresh_cpu.h"

#include "thresh_apply_cpu.h"
#include "thresh_apply_cpu.h"
 
int main_full(int argc, char **argv) {
  if (argc != 6) {
    fprintf(stderr, "\n  test_haar_thresh_quantUnif_cpu <steps> <input data> <output data> compRatio bits\n"
            "  Negative steps does inverse transform\n\n");
    return 1;
  }

  float *data, compRatio;
  int width, height, stepCount, bits;
  if (1 != sscanf(argv[1], "%d", &stepCount)) {
    printf("Invalid step count: \"%s\"\n", argv[1]);
    return 1;
  }
stepCount = abs(stepCount);  // We are using pos in the forward direction and neg in the reverse dir.
  
  const char *inputFile = argv[2], *outputFile = argv[3];
  
  if (1 != sscanf(argv[4], "%g", &compRatio)) {
	  printf("Invalid step count: \"%s\"\n", argv[4]);
	  return 1;
  }
  if (1 != sscanf(argv[5], "%d", &bits)) {
	  printf("Invalid bits: \"%s\"\n", argv[5]);
	  return 1;
  }
  printf("Reading %s...", inputFile);
  fflush(stdout);
  
  if (!readDataFile(inputFile, &data, &width, &height)) return 1;
  printf("%d x %d\n", width, height);
  fflush(stdout);

  // pad the size to a power of two
  int size = (width > height) ? width : height;
  int padSize = dwt_padded_length(size, 0, true);
  if (width != size || height != size) 
  {
    printf("Padding to %dx%d\n", padSize, padSize);
    fflush(stdout);
    float *padData = dwt_pad_2d(height, width, width, data,
				padSize, padSize, padSize, NULL, REPEAT);
    delete[] data;
    data = padData;
    size = width = height = padSize;
  }

  if (width != height) {
    printf("Only square matrices are supported.\n");
    return 1;
  }


   haar_2d(size, data, false, stepCount);

    float maxVal;
    float threshold = thresh_cpu(size, data, compRatio, &maxVal);  // Calculate the threshold
    thresh_apply_cpu(size, data, bits, threshold, maxVal);         // Apply threshold and uniform quantization
    //dquant_unif_cpu(size, data, bits, threshold, maxVal,1);      // reverse quantization
    haar_2d(size, data, true, stepCount);	           // Take the inverse transform
  // printMatrix(width, height, data);
  printf("Writing...\n");
  fflush(stdout);
  if (writeDataFile(outputFile, data, size, size, true))      // Write the reconstructed image!
    printf("%s written\n", outputFile);

  delete[] data;
  
  return 0;
}


int main(int argc, char **argv) {
  // testMisc();
  // testPad();

  main_full(argc, argv);

  return 0;
}
