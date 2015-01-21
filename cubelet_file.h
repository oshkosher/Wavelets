#ifndef __CUBELET_FILE__
#define __CUBELET_FILE__

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include "wavelet_compress.pb.h"


#include "wavelet.h"  // "Cube" class is defined here

/*

  File format:

  Header (64 bytes)
    24 bytes: text description ending with newline (so you can use 'head -1'
              to check it).

    40 bytes: unused

  For each cubelet:
    4 bytes: length of encoded CubeletBuffer
    <variable bytes> encoded CubeletBuffer object
    <variable bytes> raw data

  Cubelet index footer:
    4 bytes: all bits set (0xFFFFFFFF)
    encoded CubeletIndexBuffer
    4 bytes: length of encoded CubeletIndexBuffer
    4 bytes: "cube"

  The 4 byte "cube" suffix marks this file as a completed cubelet file.
  Without that, it is just a truncated cubelet file without an index.

  To quickly read the index, seek to 8 bytes from the end, read the length
  of the CubeletIndexBuffer and the "cube" suffix. Verify that the last 4
  bytes are "cube", and seek back to the start of the encoded
  CubeletIndexBuffer.
*/


class CubeletStreamWriter {

 public:

  // or NULL or "-" to write to stdout
  bool open(const char *filename);

  bool addCubelet(const Cube *cubelet);

  bool close();

  CubeletStreamWriter() {
    outf = NULL;
  }

  ~CubeletStreamWriter() {
    if (outf) close();
  }

 private:
  FILE *outf;

  // metadata for all the cubelets that have been written to this stream
  CubeletIndexBuffer index;
};


class CubeletStreamReader {
 public:
  typedef scu_wavelet::int3 int3;

  // or NULL or "-" to read from stdin
  bool open(const char *filename);

  bool isOpen() {return inf != NULL;}

  // Read the next cubelet in the stream and fill 'cube' with all
  // the metadata for it. Do not read the content data for the cubelet.
  // If we reach the end of the file (or empty cubelet that marks the end
  // of the file) return false, otherwise true;
  bool next(Cube *cube);

  // Return the content data for the current cubelet.
  // If data is NULL, memory for the data will be allocated with malloc().
  // Otherwise, it will be written to the given pointer.
  // Since this is designed to operate on streaming data, this can only
  // be called once for each cubelet (otherwise we would need to seek
  // backwards, which isn't possible with stdin, or buffer the data, which
  // would be a waste most of the time).
  // Return the address of where the data was written, or NULL on error.
  void *getRawData(void *data = NULL);

  // read the data directly into a cubelet, which may have padding between
  // data rows
  bool getCubeData(Cube *cube);

  void close();

  CubeletStreamReader() {
    inf = NULL;
    eofReached = false;
    dataSizeBytes = 0;
    dataHasBeenRead = false;
  }

  ~CubeletStreamReader() {
    close();
  }

 private:
  FILE *inf;

  bool eofReached;

  // size of the data for the curent cubelet, so we know how far ahead
  // to seek to read the next
  uint32_t dataSizeBytes;

  // Set this if the data for this cubelet has been read, so we know if
  // the file pointer is at the beginning or end of the data.
  bool dataHasBeenRead;

  int3 currentCubeSize;

  template<class NUM>
    void copyIntoTypedCube(CubeNum<NUM> *cube, void *rawData,
                         int3 inputDataSize) {
    if (cube->data() == NULL) cube->allocate();
    cube->copyFrom((NUM*)rawData, inputDataSize);
  }
};


#endif // __CUBELET_FILE__
