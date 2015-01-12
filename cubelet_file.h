#ifndef __CUBELET_FILE__
#define __CUBELET_FILE__

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include "wavelet_compress.pb.h"

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


struct Cubelet {
  unsigned width, height, depth;
  unsigned xOffset, yOffset, zOffset;

  enum DataType {
    CUBELET_UINT8,
    CUBELET_FLOAT32
  };

  DataType datatype;

  uint64_t dataFileOffset;

  // if data is not null, the destructor will call free() on it
  void *data;

  void setSize(int w=1, int h=1, int d=1) {
    width = w;
    height = h;
    depth = d;
  }

  void setOffset(int x=0, int y=0, int z=0) {
    xOffset = x;
    yOffset = y;
    zOffset = z;
  }

  Cubelet() {
    width = height = depth = 1;
    xOffset = yOffset = zOffset = 0;
    datatype = CUBELET_FLOAT32;
    data = NULL;
    dataFileOffset = 0;
  }

  ~Cubelet() {
    if (data) free(data);
  }

};  


class CubeletStreamWriter {
 public:

  // or NULL to write to stdout
  bool open(const char *filename);

  bool addCubelet(const Cubelet *cubelet);

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

  // or NULL to read from stdin
  bool open(const char *filename);

  // Read the next cubelet in the stream and fill 'cube' with all
  // the metadata for it. Do not read the content data for the cubelet.
  // If we reach the end of the file (or empty cubelet that marks the end
  // of the file) return false, otherwise true;
  bool next(Cubelet *cube);

  // Return the content data for the current cubelet.
  // If data is NULL, memory for the data will be allocated with malloc().
  // Otherwise, it will be written to the given pointer.
  // Since this is designed to operate on streaming data, this can only
  // be called once for each cubelet (otherwise we would need to seek
  // backwards, which isn't possible with stdin, or buffer the data, which
  // would be a waste most of the time).
  // Return the address of where the data was written, or NULL on error.
  void *getData(void *data = NULL);

  void close();

  CubeletStreamReader() {
    inf = NULL;
    eofReached = false;
    indexOffset = 0;
    dataSizeBytes = 0;
    dataHasBeenRead = false;
  }

  ~CubeletStreamReader() {
    close();
  }

 private:
  FILE *inf;

  bool eofReached;

  // Read from the header of the file, this is the index of the
  // CubeletIndexBuffer at the end of the file. If zero, then it is not present.
  uint64_t indexOffset;

  // size of the data for the curent cubelet, so we know how far ahead
  // to seek to read the next
  uint32_t dataSizeBytes;

  // Set this if the data for this cubelet has been read, so we know if
  // the file pointer is at the beginning or end of the data.
  bool dataHasBeenRead;

};


#endif // __CUBELET_FILE__
