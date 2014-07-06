#include "cuda.h"
#include "cucheck.h"
#include "data_io.h"
#include "nixtimer.h"
#include "dwt_cpu.h"
#include "dwt_gpu.h"
#include "transpose_gpu.h"


bool similar(float a, float b) {
  return fabs(a-b) < .00001;
}


void printHelp() {
  fprintf(stderr, 
          "\n"
          "  haar [-inverse] [-text] <input datafile> <output datafile> [steps]\n"
          "  Do a Haar discrete wavelet transform.\n"
          "    -inverse : invert the transform\n"
          "    -text : output data in text format rather than binary\n"
          "  By default, one transformation step will be done.\n"
          "\n");
  exit(1);
}


int main(int argc, char **argv) {
  if (argc < 3) printHelp();

  bool inverse = false, textOutput = false;
  int argNo = 1;

  while (argNo < argc && argv[argNo][0] == '-') {
    if (!strcmp(argv[argNo], "-inverse")) {
      inverse = true;
      argNo++;
    }

    else if (!strcmp(argv[argNo], "-text")) {
      textOutput = true;
      argNo++;
    }

    else printHelp();
  }

  // not enough arguments for the input file and output file
  if (argNo+2 > argc) printHelp();

  // read the input file
  const char *inputFilename = argv[argNo++];
  const char *outputFilename = argv[argNo++];
  int stepCount = 1;

  if (argNo < argc) {
    const char *stepsArg = argv[argNo++];
    if (1 != sscanf(stepsArg, "%d", &stepCount)) {
      printf("Invalid step count: \"%s\"\n", stepsArg);
      return 1;
    }
  }
  if (argNo < argc) printHelp();

  float *data_cpu, *data_gpu, elapsed;
  int width, height;
  printf("Reading %s...", inputFilename);
  fflush(stdout);
  if (!readDataFile(inputFilename, &data_cpu, &width, &height)) return 1;
  printf("%d x %d\n", width, height);
  fflush(stdout);

  if (width != height) {
    printf("Error: only square data is currently supported.\n");
    return 1;
  }

  int size = width;
  
  CUCHECK(cudaSetDevice(0));

  // make a copy of the data for the GPU to use
  data_gpu = new float[size*size];
  memcpy(data_gpu, data_cpu, sizeof(float)*height*width);

  // run the CPU version of the algorithm
  printf("CPU: "); fflush(stdout);
  elapsed = haar_not_lifting_2d(size, data_cpu, inverse, stepCount);
  printf("%.6f ms\n", elapsed);

  // run the GPU version of the algorithm
  elapsed = haar_not_lifting_2d_cuda(size, data_gpu, inverse, stepCount);
  printf("CUDA: %.6f ms\n", elapsed);

  /*
    // try a variety of thread block sizes
  float *data_gpu_copy = new float[size*size];
  memcpy(data_gpu_copy, data_gpu, sizeof(float)*height*width);
  for (int threadBlockSize = 32; threadBlockSize <= 1024; threadBlockSize*=2) {
    printf("Thread block size: %d\n", threadBlockSize);
    memcpy(data_gpu, data_gpu_copy, sizeof(float)*height*width);
    elapsed = haar_not_lifting_2d_cuda(size, data_gpu, inverse, stepCount,
                                       threadBlockSize);
    printf("CUDA: %.6f ms\n", elapsed);
  }
  delete[] data_gpu_copy
  */

  double totalErr = 0;
  for (int i=0; i < size*size; i++) {
    totalErr += fabs(data_cpu[i] - data_gpu[i]);
  }
  
  double averageErr = totalErr / (size*size);

  if (averageErr < 0.000001) {

    // if the CPU version and the GPU version produced similar results,
    // output the requested file.
    writeDataFile(outputFilename, data_gpu, size, size, !textOutput);
    printf("Wrote %s\n", outputFilename);

  } else {

    // if the results look bad, output two files; one with the CPU results
    // and one with the GPU results.
    printf("Average error = %.7f\n", averageErr);
    
    /*
      printf("CPU:\n");
      printMatrix(width, height, data_cpu);
      
      printf("GPU:\n");
      printMatrix(width, height, data_gpu);
    */
    
    writeDataFile("err_cpu.data", data_cpu, size, size, !textOutput);
    writeDataFile("err_gpu.data", data_gpu, size, size, !textOutput);
    printf("Wrote err_cpu.data and err_gpu.data\n");
  }

  return 0;
}
