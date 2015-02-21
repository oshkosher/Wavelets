#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>

#include "wavelet.h"
#include "cubelet_file.h"
using namespace scu_wavelet;

#include "mat.h"

/*
  On Windows, this must be built from a 64-bit Visual Studio command
  line (using Makefile.nmake), since it links with Matlab's 64 Windows
  libraries.

  Matlab's "mat.h" is found at MATLAB/extern/include/mat.h

  To set variables to enable 64-bit code:
    C:\mvs2012\VC\bin\x86_amd64\vcvarsx86_amd64.bat

  And add MATLAB\bin\win64 to the PATH:
    set PATH=%PATH%;"c:\Program Files\MATLAB\R2014a\bin\win64"

*/


using namespace std;


void printHelp() {
  printf("\n  mat2cube <output_cube_file> <input_mat_file> [<mat_variable_name>]\n\n");
  exit(1);
}


void getAllVariableNames(vector<string> &vars, MATFile *matfile) {
  vars.clear();

  while (true) {
    const char *name;
    mxArray *array = matGetNextVariableInfo(matfile, &name);
    if (array == NULL) break;
    vars.push_back(name);
    mxDestroyArray(array);
  }
}


size_t getCount(mwSize dimCount, const mwSize *dims) {
  size_t count = 1;
  for (int i=0; i < dimCount; i++) {
    count *= dims[i];
  }
  return count;
}


void scanUint8(unsigned char *data, mwSize dimCount, const mwSize *dims) {
}

void scanInt8(signed char *data, mwSize dimCount, const mwSize *dims) {
  size_t count = getCount(dimCount, dims);
  size_t zeroCount = 0, ltZeroCount = 0;
  for (size_t i=0; i < count; i++) {
    signed char c = data[i];
    if (c < 0) {
      ltZeroCount++;
    } else if (c == 0) {
      zeroCount++;
    }
    printf("%d\n", (int)c);
  }

  printf("%d zero, %d < zero, %d > zero\n",
         (int)zeroCount, (int)ltZeroCount,
         (int)(count - ltZeroCount - zeroCount));
}



int main(int argc, char **argv) {

  if (argc < 3 || argc > 4) printHelp();

  const char *outputCubeFilename = argv[1];
  const char *inputMatlabFilename = argv[2];
  const char *matVariableName = argv[3];

  MATFile *matfile = matOpen(inputMatlabFilename, "r");
  if (!matfile) {
    printf("ERROR: cannot open \"%s\".\n", inputMatlabFilename);
    return 1;
  }

  vector<string> variables;
  if (matVariableName == NULL) {
    getAllVariableNames(variables, matfile);

    if (variables.size() == 0) {
      printf("No variables found in \"%s\":\n", inputMatlabFilename);
      return 1;
    }

    /*    
    printf("Variables found in \"%s\":\n", inputMatlabFilename);
    for (size_t i=0; i < variables.size(); i++) {
      printf("  %s\n", variables[i].c_str());
    }
    */
    
    if (variables.size() > 1) {
      printf("Choose a variable and add it to the command line.\n");
      return 1;
    }

    matVariableName = variables[0].c_str();
  }

  mxArray *var = matGetVariable(matfile, matVariableName);

  if (var == NULL) {
    printf("\"%s\" variable not found in \"%s\"\n", matVariableName,
           inputMatlabFilename);
    return 1;
  }

  mwSize dimCount = mxGetNumberOfDimensions(var);
  const mwSize *dims = mxGetDimensions(var);
  printf("size = ");
  for (int i = 0; i < dimCount; i++) {
    if (i > 0) putchar('x');
    printf("%d", (int)dims[i]);
  }
  putchar('\n');
  if (dimCount > 3) {
    printf("Cubelet file cannot handle data with > 3 dimensions.\n");
    return 1;
  }

  printf("element size = %d bytes\n", (int)mxGetElementSize(var));

  Cube cube;
  cube.size = int3(1,1,1);
  if (dimCount > 0) {
    cube.size.x = dims[0];
    if (dimCount > 1) {
      cube.size.y = dims[1];
      if (dimCount > 2) {
        cube.size.z = dims[2];
      }
    }
  }
  cube.totalSize = cube.size;
  cube.data_ = mxGetData(var);

  mxClassID dataType = mxGetClassID(var);
  switch (dataType) {
  case mxUINT8_CLASS:
    printf("datatype: uint8\n");
    cube.datatype = WAVELET_DATA_UINT8;
    break;
  case mxINT8_CLASS:
    printf("datatype: int8\n");
    cube.datatype = WAVELET_DATA_UINT8;
    break;
  default:
    printf("datatype: %d\n", (int)dataType, dimCount, dims);
  }

  CubeletStreamWriter cubeWriter;

  if (!cubeWriter.open(outputCubeFilename)) return 1;
  if (!cubeWriter.addCubelet(&cube)) return 1;
  if (!cubeWriter.close()) return 1;

  mxDestroyArray(var);

  printf("Wrote \"%s\"\n", outputCubeFilename);

  return 0;
}

