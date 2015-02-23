#include <string>
#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <cstdint>
#include <vector>
#include <algorithm>
#include "cubelet_file.h"

using namespace std;
using namespace scu_wavelet;


/*
struct Options {
  std::string inputFile;
  int3 inputShape;

  std::string outputFile;
};
*/
  

void printHelp();
bool isCubeletFilename(const char *filename);

void cube_list(const char *cubelet_file);
void cube_cat(const char *cubelet_file, int3 id);
void cube_frombytes(const char *output_cubelet_file, int3 size,
                    const char *raw_input_file);
bool cube_split(const char *output_cubelet_file, int3 size, int overlap,
                const char *input_cubelet_file);

template<class CubeType>
bool cube_split_typed(CubeType *allCubes, int allCubesCount,
		      CubeletStreamReader &streamIn,
		      const char *output_cubelet_file,
		      int3 output_cubelet_size, int overlap);
bool cube_extract(const char *output_file, bool do_append,
		  const char *input_file, const std::vector<int3> &id_list);
bool cube_partial(const char *output_file, const char *input_file,
		  int3 offset, int3 size, const std::vector<int3> &id_list);
		  

template<class T>
void copy2d(T *dest, int destWidth, const T *src, int srcWidth,
	    int copyWidth, int copyHeight);


int main(int argc, const char **argv) {
  if (argc < 3) printHelp();

  const char *command = argv[1];
  const char *cubelet_file = argv[2];

  if (!strcmp(command, "list")) {
    if (argc != 3) printHelp();
    cube_list(cubelet_file);
  }

  else if (!strcmp(command, "cat")) {
    if (argc > 4) printHelp();
    int3 id(-1,-1,-1);
    if (argc > 3) {
      if (id.parse(argv[3]) != 3) printHelp();
    }
    cube_cat(cubelet_file, id);
  }

  else if (!strcmp(command, "frombytes")) {
    int3 size;
    if (argc != 5) printHelp();
    if (size.parse(argv[3]) != 3) printHelp();
    const char *input_file = argv[4];
    cube_frombytes(cubelet_file, size, input_file);
  }

  else if (!strcmp(command, "split")) {
    int3 size;
    int overlap;
    if (argc != 6) printHelp();
    if (size.parse(argv[3]) != 3) printHelp();
    if (1 != sscanf(argv[4], "%d", &overlap)) printHelp();
    const char *input_file = argv[5];
    cube_split(cubelet_file, size, overlap, input_file);
  }

  else if (!strcmp(command, "extract")) {
    vector<int3> id_list;
    bool do_append = false;
    const char *output_file, *input_file;
    int argno = 2;
    if (!strcmp(argv[argno], "-append")) {
      do_append = true;
      argno++;
    }
    if (argc-argno < 3) printHelp();
    output_file = argv[argno++];
    input_file = argv[argno++];
    while (argno < argc) {
      int3 id;
      if (3 != id.parse(argv[argno])) {
	printf("Invalid cubelet id \"%s\" (should be in the form x,y,z)\n",
	       argv[argno]);
	return 1;
      }
      id_list.push_back(id);
      argno++;
    }

    cube_extract(output_file, do_append, input_file, id_list);
  }

  else if (!strcmp(command, "partial")) {
    if (argc < 6) printHelp();
    const char *output_file = argv[2];
    const char *input_file = argv[3];
    int3 offset, size;
    if (3 != offset.parse(argv[4]) || !(offset >= int3(0,0,0))) {
      printf("Invalid <x,y,z> offset: %s\n", argv[4]);
      return 1;
    }
    if (3 != size.parse(argv[5]) || !(size > int3(0,0,0))) {
      printf("Invalid <x,y,z> size: %s\n", argv[5]);
      return 1;
    }
    vector<int3> id_list;
    for (int argno = 6; argno < argc; argno++) {
      int3 id;
      if (3 != id.parse(argv[argno])) {
	printf("Invalid <x,y,z> cubelet id: %s\n", argv[argno]);
	return 1;
      }
      id_list.push_back(id);
    }
    cube_partial(output_file, input_file, offset, size, id_list);
  }
    

  else printHelp();

  // deallocate static protobuf data
  google::protobuf::ShutdownProtobufLibrary();
  
  return 0;
}


void printHelp() {
  printf("\n"
"  Operations on cubelet files.\n"
"\n"
"  Where <x,y,z> is found, it specifies the id of a cubelet within a\n"
"  cublet file. When omitted, most operations default to the first cubelet\n"
"  in the file.\n"
"\n"
"  cube list <cube>\n"
"    List all the cubelets in the given file.\n"
"\n"
"  cube cat <cube> [<x,y,z>]\n"
"    Print the value of each pixel, one per line, to stdout. If <x,y,z>\n"
"    is specified, use that cubelet, otherwise the first one in the\n"
"    file.\n"
"\n"
"  cube frombytes <outputCube> <x,y,z> <inputFile>\n"
"    Given an input file consisting of raw unsigned bytes, convert it\n"
"    to a single cubelet file.\n"
"\n"
"  cube split <outputCube> <W,H,D> <overlap> <inputCube>\n"
"    Given an existing cubelet file, split it into smaller cubelets of\n"
"    size WxHxD, with 'overlap' pixels of overlap wherever cubelets meet.\n"
"    This accepts input cubelet files in two formats:\n"
"      1. One large cubelet.\n"
"      2. Many 2-d cubelets (z dimension == 1), all the same size and type,\n"
"         and their Z coordinates range from 0..N.\n"
"\n"
"  cube extract [-append] <outputCube> <inputCube> <x,y,z> [<x,y,z> ...]\n"
"    Copy one or more cublets from one file to another.\n"
"    With the \"-append\" flag, add to <outputCube>, otherwise overwrite.\n"
"\n"
"  cube partial <outputCube> <inputCube> <offset> <size> [<x,y,z> ...]\n"
"    Extract part of a cubelet, starting at the given offset, of the given size.\n"
"    Both <offset> and <size> are x,y,z triplets.\n"
"    If no cubelet id is specified, extract the same part from every cubelet.\n"
"\n"
         );
  exit(1);
}


// return true iff the given filename has a ".cube" suffix.
bool isCubeletFilename(const char *f) {
  string filename = f;
  return filename.length() >= 5 &&
    filename.substr(filename.length() - 5, 5) == ".cube";
}


void cube_list(const char *filename) {

  CubeletStreamReader in;

  if (!in.open(filename)) return;

  Cube cube;
  
  while (true) {
    if (!in.next(&cube)) break;

    const char *compressed = cube.isWaveletCompressed ? " compressed" : "";

    printf("Cubelet %d,%d,%d: %dx%dx%d of %s, %d%s bytes at %" PRIu64 "\n",
           cube.parentOffset.x, cube.parentOffset.y, cube.parentOffset.z,
           cube.width(), cube.height(), cube.depth(),
           waveletDataTypeName(cube.datatype),
           cube.getSizeInBytes(), compressed, cube.dataFileOffset);
  }
}


template<class T>
void printAll(const T *cube) {
  int count = cube->count();
  const char *format = cube->printfFormatSimple();
  for (int i=0; i < count; i++) {
    printf(format, cube->get(i,0,0));
    putchar('\n');
  }
}


void cube_cat(const char *cubelet_file, int3 id) {
  CubeletStreamReader in;

  if (!in.open(cubelet_file)) return;

  // find the requested cubelet
  Cube cube;
  if (!in.find(&cube, id)) {
    fprintf(stderr, "Cubelet %d,%d,%d not found.\n", id.x, id.y, id.z);
    return;
  }

  if (!in.getCubeData(&cube)) return;

  in.close();

  if (cube.isWaveletCompressed) {
    fprintf(stderr, "Cannot print compressed cubelets.\n");
    return;
  }
  
  switch (cube.datatype) {
  case WAVELET_DATA_UINT8:
    printAll((CubeByte*)&cube);
    break;
  case WAVELET_DATA_INT32:
    printAll((CubeInt*)&cube);
    break;
  case WAVELET_DATA_FLOAT32:
    printAll((CubeFloat*)&cube);
    break;
  default:
    fprintf(stderr, "Unrecognized datatype: %d\n", cube.datatype);
  }
}


void cube_frombytes(const char *output_cubelet_file, int3 size,
                    const char *raw_input_file) {

  if (!(size > int3(0,0,0))) {
    fprintf(stderr, "Input shape is invalid.\n");
    return;
  }
  
  FILE *inf = fopen(raw_input_file, "rb");
  if (!inf) {
    fprintf(stderr, "Failed to read %s\n", raw_input_file);
    return;
  }

  CubeByte cube;
  cube.size = size;
  cube.allocate();
  uint64_t dataSize = (uint64_t)size.x * size.y * size.z;
  uint64_t bytesRead;
  bytesRead = fread(cube.data(), 1, dataSize, inf);
  if (bytesRead != dataSize) {
    fprintf(stderr, "Data size mismatch: expected %" PRIu64 " bytes, "
            "got %" PRIu64 "\n", dataSize, bytesRead);
  }

  CubeletStreamWriter out;
  if (!out.open(output_cubelet_file)) return;
  if (!out.addCubelet(&cube)) return;
  out.close();

  printf("Wrote %s\n", output_cubelet_file);
}


struct OrderCubesByZLevel {
  bool operator () (const Cube &a, const Cube &b) const {
    if (a.parentOffset.z == b.parentOffset.z) {
      if (a.parentOffset.y == b.parentOffset.y) {
        return a.parentOffset.x < b.parentOffset.x;
      } else {
        return a.parentOffset.y < b.parentOffset.y;
      }
    } else {
      return a.parentOffset.z < b.parentOffset.z;
    }
  }
};


// Check that all the cubelets are 2-d frames and that there are no
// frames missing.
bool cubeletsAre2DFrames(const Cube *allCubes, int count) {
  const Cube &first = allCubes[0];
  if (first.size.z != 1) {
    printf("Error: all cubes should have thickness of 1, but cube %s "
           "is %dx%dx%d\n", first.getId(), first.size.x, first.size.y,
           first.size.z);
    return false;
  }

  for (int i=0; i < count; i++) {
    const Cube &c = allCubes[i];
    if (c.size != first.size) {
      printf("Error: cube %s should have size %dx%dx%d, but is %dx%dx%d\n",
             c.getId(), first.size.x, first.size.y, first.size.z,
             c.size.x, c.size.y, c.size.z);
      return false;
    }

    if (c.parentOffset.z != i) {
      printf("Error: cubelet 0,0,%d missing.\n", i);
      return false;
    }

    if (c.parentOffset.x != 0 || c.parentOffset.y != 0) {
      printf("Error: cube %s should have id 0,0,%d\n",
             c.getId(), i);
      return false;
    }

    if (c.datatype != first.datatype) {
      printf("All cubes should have datatype %s, but cube %s has type %s\n",
             waveletDataTypeName(first.datatype), c.getId(),
             waveletDataTypeName(c.datatype));
      return false;
    }

    if (c.maxPossibleValue != first.maxPossibleValue) {
      printf("All cubes should have maxPossibleValue %d, but cube %s has "
             "max=%d\n",
             first.maxPossibleValue, c.getId(), c.maxPossibleValue);
      return false;
    }
  }

  return true;
}


bool cube_split(const char *output_cubelet_file, int3 size, int overlap,
                const char *input_cubelet_file) {

  vector<Cube> allCubes;
  CubeletStreamReader in;
  if (!in.open(input_cubelet_file)) return false;
  in.listCubelets(allCubes);

  OrderCubesByZLevel zOrder;
  sort(allCubes.begin(), allCubes.end(), zOrder);

  /*
  for (size_t i=0; i < allCubes.size(); i++) {
    Cube c = allCubes[i];
    printf("Cube %d: %d,%d,%d\n", (int)i, c.parentOffset.x, c.parentOffset.y,
           c.parentOffset.z);
  }
  */

  if (allCubes.size() == 0) {
    printf("Error: %s is empty\n", input_cubelet_file);
    return false;
  }

  if (!(size > int3(0,0,0))) {
    printf("Invalid output cubelet size: %d,%d,%d (must all be > 0)\n",
           size.x, size.y, size.z);
    return false;
  }

  WaveletDataType datatype = allCubes[0].datatype;
  int count = (int) allCubes.size();

  switch (datatype) {
  case WAVELET_DATA_UINT8:
    return cube_split_typed((CubeByte*)allCubes.data(), count, in,
			    output_cubelet_file, size, overlap);
			    
    break;
  case WAVELET_DATA_INT32:
    return cube_split_typed((CubeInt*)allCubes.data(), count, in,
			    output_cubelet_file, size, overlap);
    break;
  case WAVELET_DATA_FLOAT32:
    return cube_split_typed((CubeFloat*)allCubes.data(), count, in,
			    output_cubelet_file, size, overlap);
    break;
  default:
    printf("Unrecognized data type: %d\n", datatype);
    return false;
  }

}


template<class CubeType>
bool cube_split_typed(CubeType *allCubes, int allCubesCount,
		      CubeletStreamReader &streamIn,
		      const char *output_cubelet_file,
		      int3 size, int overlap) {

  CubeletStreamWriter streamOut;
  if (!streamOut.open(output_cubelet_file)) return false;
  
  bool isSingleCube;

  // total size of the input data, as if it was one big cubelet
  int3 inputSize;

  if (allCubesCount == 1) {
    // input is one big cube
    inputSize = allCubes[0].size;
    isSingleCube = true;
    // read the data for the one big cubelet
    streamIn.getCubeData(&allCubes[0]);
  } else {
    // input is many flat frames
    if (!cubeletsAre2DFrames(allCubes, allCubesCount)) return false;
    inputSize = allCubes[0].size;
    inputSize.z = allCubesCount;
    isSingleCube = false;
  }

  // number of cubelets in each direction
  int3 cubeCount;
  cubeCount.x = (inputSize.x - overlap - 1) / (size.x - overlap) + 1;
  cubeCount.y = (inputSize.y - overlap - 1) / (size.y - overlap) + 1;
  cubeCount.z = (inputSize.z - overlap - 1) / (size.z - overlap) + 1;

  // set up the temporary cubelets
  int tempCubeCount = cubeCount.x * cubeCount.y;
  CubeType *tempCubes = new CubeType[tempCubeCount];
  for (int y=0; y < cubeCount.y; y++) {
    for (int x=0; x < cubeCount.x; x++) {
      int cubeNo = y * cubeCount.x + x;
      CubeType *c = &tempCubes[cubeNo];
      c->setType();

      c->size = size;
      // set the x and y position of the cubes; the z position will
      // be set during processing
      c->parentOffset = int3(x * (size.x-overlap),
			     y * (size.y-overlap),
			     0);

      c->maxPossibleValue = allCubes[0].maxPossibleValue;

      // trim cubes that stick out over the end
      int end = c->size.x + c->parentOffset.x;
      if (end > inputSize.x)
	c->size.x = inputSize.x - c->parentOffset.x;

      end = c->size.y + c->parentOffset.y;
      if (end > inputSize.y)
	c->size.y = inputSize.y - c->parentOffset.y;

      c->allocate();
    }
  }

  // how many z-levels of data have been added to the temp cubelets
  int currentDepth = 0;

  for (int frameNo=0; frameNo < inputSize.z; frameNo++) {
    printf("\rFrame %d of %d", frameNo+1, inputSize.z);
    fflush(stdout);
    
    // get a pointer to the next frame of data

    typename CubeType::MyType *inputData, *outputData;

    // copy the data into subcubes

    for (CubeType *c = tempCubes; c < tempCubes + tempCubeCount; c++) {

      if (isSingleCube) {
	inputData = allCubes->pointer(c->parentOffset.x,
				      c->parentOffset.y, frameNo);
      } else {
	// read the data for this frame
	streamIn.getCubeData(&allCubes[frameNo]);
	inputData = allCubes[frameNo].pointer(c->parentOffset.x,
					      c->parentOffset.y, 0);
      }

      outputData = c->pointer(0, 0, currentDepth);
      
      copy2d(outputData, c->size.x,
	     inputData, inputSize.x,
	     c->size.x, c->size.y);

      // deallocate the data for this frame
      if (!isSingleCube) {
	allCubes[frameNo].deallocate();
      }

    }

    currentDepth++;
    if (currentDepth == size.z) {

      // flush temp cubelets

      for (CubeType *c = tempCubes; c < tempCubes + tempCubeCount; c++) {
	c->size.z = currentDepth;

	if (!streamOut.addCubelet(c)) return false;

	// set the z offset for the next copy of this cubelet
	c->parentOffset.z = frameNo + 1 - overlap;
      }

      // reset the depth
      currentDepth = 0;

      // if this isn't the end of the frames, back up a bit 
      if (frameNo < inputSize.z-1) frameNo -= overlap;
    }

  }
  putchar('\n');

  // write out leftover frames
  if (currentDepth > 0) {
    for (CubeType *c = tempCubes; c < tempCubes + tempCubeCount; c++) {
      c->size.z = currentDepth;
      if (!streamOut.addCubelet(c)) return false;
    }
  }

  delete[] tempCubes;
  streamOut.close();

  return true;
}


template<class T>
void copy2d(T *dest, int destWidth, const T *src, int srcWidth,
	    int copyWidth, int copyHeight) {

  for (int y = 0; y < copyHeight; y++) {
    memcpy(dest, src, sizeof(T) * copyWidth);
    dest += destWidth;
    src += srcWidth;
  }

}


bool cube_extract(const char *output_file, bool do_append,
		  const char *input_file, const std::vector<int3> &id_list) {
  
  CubeletStreamWriter writer;
  if (!writer.open(output_file, do_append)) return false;

  CubeletStreamReader reader;
  if (!reader.open(input_file)) return false;

  bool any_errors = false;

  for (size_t i=0; i < id_list.size(); i++) {
    Cube cube;
    int3 id = id_list[i];

    if (!reader.find(&cube, id)) {
      printf("Cannot find cubelet %d,%d,%d. Skipping.\n", id.x, id.y, id.z);
      continue;
    }
    
    if (!reader.getCubeData(&cube)) continue;
    
    writer.addCubelet(&cube);
    printf("Extracted %d,%d,%d: %d bytes\n", id.x, id.y, id.z,
	   cube.getSizeInBytes());
    fflush(stdout);
  }
  
  reader.close();
  writer.close();

  return !any_errors;
}


template<class CubeType>
void copy_partial_typed(const CubeType *cube, CubeletStreamWriter &writer,
			int3 offset, int3 size) {
  const typename CubeType::MyType *readp;
  typename CubeType::MyType *writep;

  if (!(offset < cube->size)) {
    printf("Cubelet %s too small for offset %d,%d,%d\n",
	   cube->getId(), size.x, size.y, size.z);
    return;
  }
  
  CubeType outCube = *cube;
  outCube.size = outCube.totalSize = size;
  outCube.data_ = NULL;
  outCube.allocate();

  int width;
  if (offset.x + size.x > cube->size.x) {
    width = cube->size.x - offset.x;
  } else {
    width = size.x;
  }

  for (int z=0; z < size.z && z + offset.z < cube->size.z; z++) {
    for (int y=0; y < size.y && y + offset.y < cube->size.y; y++) {
      readp = cube->pointer(offset.x, offset.y + y, offset.z + z);
      writep = outCube.pointer(0, y, z);
      memcpy(writep, readp, width * (sizeof *readp));
    }
  }

  writer.addCubelet(&outCube);
}


void copy_partial(const Cube *cube, CubeletStreamWriter &writer,
		  int3 offset, int3 size) {

  switch(cube->datatype) {
  case WAVELET_DATA_UINT8:
    copy_partial_typed((const CubeByte*)cube, writer, offset, size);
    break;
  case WAVELET_DATA_INT32:
    copy_partial_typed((const CubeInt*)cube, writer, offset, size);
    break;
  case WAVELET_DATA_FLOAT32:
    copy_partial_typed((const CubeFloat*)cube, writer, offset, size);
    break;
  default:
    printf("Error on cube %s: unrecognized datatype %d\n",
	   cube->getId(), cube->datatype);
  }
}


bool cube_partial(const char *output_file, const char *input_file,
		  int3 offset, int3 size, const std::vector<int3> &id_list) {
  
  CubeletStreamWriter writer;
  if (!writer.open(output_file)) return false;

  CubeletStreamReader reader;
  if (!reader.open(input_file)) return false;

  Cube cube;

  if (id_list.empty()) {
    // loop through all the input cubelets
    while (reader.next(&cube)) {
      cube.allocate();
      reader.getRawData(cube.data_);
      copy_partial(&cube, writer, offset, size);
      cube.deallocate();
    }
  } else {
    // extract just the listed cubelets
    for (size_t i=0; i < id_list.size(); i++) {
      int3 id = id_list[i];
      if (!reader.find(&cube, id)) {
	printf("Cubelet %d,%d,%d not found.\n", id.x, id.y, id.z);
      } else {
	cube.allocate();
	reader.getRawData(cube.data_);
	copy_partial(&cube, writer, offset, size);
	cube.deallocate();
      }
    }
  }

  reader.close();
  writer.close();

  return true;
}
		  
