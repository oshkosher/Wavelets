#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>

#include <octave/oct.h>

/*
  /usr/include/octave-3.8.2/octave/oct.h
  /usr/lib/octave/3.8.2/liboctave.dll.a
*/


using namespace std;


void printHelp() {
  printf("\n  oct2cube <output_cube_file> <input_oct_file> [<oct_variable_name>]\n\n");
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

  


int main(int argc, char **argv) {

  if (argc < 3 || argc > 4) printHelp();

  const char *outputCubeFilename = argv[1];
  const char *inputOctaveFilename = argv[2];
  const char *octVariableName = argv[3];

  MATFile *matfile = matOpen(inputMatlabFilename, "r");
  if (!matfile) {
    fprintf(stderr, "ERROR: cannot open \"%s\".\n", inputMatlabFilename);
    return 1;
  }

  vector<string> variables;
  getAllVariableNames(variables, matfile);

  printf("Variables found in \"%s\":\n", inputMatlabFilename);
  for (size_t i=0; i < variables.size(); i++) {
    printf("  %s\n", variables[i].c_str());
  }

  return 0;
}

