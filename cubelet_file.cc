#include <cassert>
#include <cstdio>
#include <cstring>
#include "cubelet_file.h"
#include "wavelet_compress.pb.h"


// 20 character header line
static const char *CUBELET_HEADER = "SCU cubelets 1.0\n";
#define HEADER_SIZE 64

static bool writeProtobuf(FILE *outf, google::protobuf::Message *message);
static bool readProtobuf(FILE *inf, google::protobuf::Message *message,
                         uint32_t *length = NULL);


// This is designed to be 64 bytes
struct CubeletHeader {
  char headerLine[24];

  // Offset into this file for the index containing the metadata
  // for all the cubelets. This will initially be zero, and will be
  // backpatched when the file is complete.
  uint64_t cubeletIndexOffset;

  char unused[32];

  CubeletHeader() {
    assert(sizeof(CubeletHeader) == HEADER_SIZE);
    strncpy(headerLine, CUBELET_HEADER, 24);
    headerLine[23] = 0;  // null-terminate it just in case

    memset(unused, 0, sizeof unused);

    cubeletIndexOffset = 0;
  }
};


/* Serialize a protocol buffer message to the given file by writing the
   size of the encoded data (4 byte unsigned integer) followed by the
   encoded data.
*/
static bool writeProtobuf(FILE *outf, google::protobuf::Message *message) {
                          
  uint32_t codedLen = message->ByteSize();
  fwrite(&codedLen, sizeof codedLen, 1, outf);

  // encode the metadata
  char *codedBuf = new char[codedLen];
  assert(codedBuf);
  if (!message->SerializeToArray(codedBuf, codedLen)) {
    fprintf(stderr, "Failed to encode metadata\n");
    return false;
  }

  // write the metadata to the file
  if (fwrite(codedBuf, 1, codedLen, outf) != codedLen) {
    fprintf(stderr, "Failed to write encoded metadata to file\n");
    return false;
  }
  delete[] codedBuf;

  return true;
}  


static bool readProtobuf(FILE *inf, google::protobuf::Message *message,
                         uint32_t *length) {

  // read the length
  uint32_t codedLen;
  if (fread(&codedLen, sizeof codedLen, 1, inf) != 1) {
    // this will happen naturally at EOF
    // fprintf(stderr, "Failed to read message buffer size\n");
    return false;
  }

  // return the length, if the caller is interested
  if (length) *length = codedLen;

  // this is a bit of a hack, but if the length is 0xFFFFFFFF,
  // consider it invalid and return false. This marks EOF.
  if (codedLen == 0xFFFFFFFF) return false;

  // allocate a buffer
  // XXX might want to have a small (32 bytes or so) buffer array on the
  // stack as a possible shortcut
  char *codedBuf = new char[codedLen];
  assert(codedBuf);

  // read the encoded data from the file
  if (fread(codedBuf, 1, codedLen, inf) != codedLen) {
    fprintf(stderr, "Failed to read metadata\n");
    return false;
  }

  // decode the metadata
  if (!message->ParseFromArray(codedBuf, codedLen)) {
    fprintf(stderr, "Failed to decode metadata\n");
    return false;
  }

  delete[] codedBuf;
  
  return true;
}


bool CubeletStreamWriter::open(const char *filename) {
  // basic assumption in a 32-bit world
  assert(sizeof(unsigned) == 4);

  if (outf) close();

  if (filename) {
    outf = fopen(filename, "wb");
    if (!outf) {
      fprintf(stderr, "Error: cannot open \"%s\" for writing.\n", filename);
      return false;
    }
  } else {
    outf = stdout;
  }

  CubeletHeader header;

  size_t bytesWritten = fwrite(&header, 1, sizeof(header), outf);
  if (bytesWritten != sizeof(header)) {
    fprintf(stderr, "Failed to write cubelet file header\n");
    if (outf != stdout) fclose(outf);
    outf = NULL;
    return false;
  }
  
  return true;
}


bool CubeletStreamWriter::addCubelet(const Cubelet *cubelet) {
  if (!outf) {
    fprintf(stderr, "Cannot add cubelets. Stream is not open yet.\n");
    return false;
  }
  
  CubeletBuffer *buffer = index.add_cubelets();

  assert(cubelet->width != 0);
  assert(cubelet->height != 0);
  assert(cubelet->depth != 0);

  buffer->set_width(cubelet->width);
  buffer->set_height(cubelet->height);
  buffer->set_depth(cubelet->depth);

  if (cubelet->x_offset) buffer->set_x_offset(cubelet->x_offset);
  if (cubelet->y_offset) buffer->set_y_offset(cubelet->y_offset);
  if (cubelet->z_offset) buffer->set_z_offset(cubelet->z_offset);

  assert(cubelet->datatype == Cubelet::CUBELET_UINT8 ||
         cubelet->datatype == Cubelet::CUBELET_FLOAT32);

  unsigned pixelSize = 0;

  switch (cubelet->datatype) {
  case Cubelet::CUBELET_UINT8:
    buffer->set_data_type(CubeletBuffer_DataType_UINT8);
    pixelSize = 1;
    break;
  case Cubelet::CUBELET_FLOAT32:
    buffer->set_data_type(CubeletBuffer_DataType_FLOAT32);
    pixelSize = sizeof(float);
    break;
  }

  unsigned dataByteCount;

  // if the data hasn't been set, don't output any data
  if (cubelet->data)
    dataByteCount = cubelet->width * cubelet->height * cubelet->depth
      * pixelSize;
  else
    dataByteCount = 0;
    
  // when the data is compressed, byte_count will be harder to compute
  buffer->set_byte_count(dataByteCount);

  writeProtobuf(outf, buffer);

  // save the position of the data so it can be included in the index
  buffer->set_data_file_offset(ftell(outf));
  
  // write the data
  if (cubelet->data) {
    if (fwrite(cubelet->data, 1, dataByteCount, outf) != dataByteCount) {
      fprintf(stderr, "Failed to write cubelet data to file\n");
      return false;
    }
  }
  
  return true;
}


bool CubeletStreamWriter::close() {
  // If the file isn't open do nothing.
  if (!outf) return true;

  // If we're writing to stdout, we can't seek back to the beginning of the
  // file to write the index offset, so there's no reason to write the index,
  // so just close the file.
  if (outf != stdout) {

    int32_t eofMarker = -1; // set all bits
    fwrite(&eofMarker, sizeof eofMarker, 1, outf);

    // get the current file position
    uint64_t indexOffset = ftell(outf);

    // write the cubelet index
    writeProtobuf(outf, &index);

    size_t headerOffset = offsetof(CubeletHeader, cubeletIndexOffset);

    if (fseek(outf, headerOffset, SEEK_SET)) {
      fprintf(stderr, "Failed to seek to file header.\n");
      return false;
    }

    if (fwrite(&indexOffset, sizeof indexOffset, 1, outf) != 1) {
      fprintf(stderr, "Failed to backpatch index offset in the file header.\n");
      return false;
    }

    fclose(outf);
  }

  outf = NULL;

  return true;
}


bool CubeletStreamReader::open(const char *filename) {
  if (inf) close();

  if (filename) {
    inf = fopen(filename, "rb");
    if (!inf) {
      fprintf(stderr, "Cannot open \"%s\" for reading.\n", filename);
      return false;
    }
  } else {
    inf = stdin;
  }

  CubeletHeader header;

  // read the header
  if (fread(&header, sizeof header, 1, inf) != 1) {
    fprintf(stderr, "Failed to read header in \"%s\".\n", filename);
    if (inf != stdout) fclose(inf);
    inf = NULL;
    return false;
  }

  if (strcmp(header.headerLine, CUBELET_HEADER)) {
    fprintf(stderr, "Input file \"%s\" is not a cubelet file.\n", filename);
    if (inf != stdout) fclose(inf);
    inf = NULL;
    return false;
  }

  // Save the offset of the index. this currently unused, but it might
  // be useful later.
  indexOffset = header.cubeletIndexOffset;

  dataSizeBytes = 0;
  dataHasBeenRead = true;

  return true;
}


bool CubeletStreamReader::next(Cubelet *cube) {
  // if already EOF, don't read more
  if (eofReached) return false;

  // if the data for this cubelet has not been read, seek past it
  if (!dataHasBeenRead) {
    if (fseek(inf, dataSizeBytes, SEEK_CUR)) {
      eofReached = true;
      return false;
    }
  }

  CubeletBuffer buf;

  // if whatever comes next is not a parseable cubelet, we're done
  if (!readProtobuf(inf, &buf)) {
    eofReached = true;
    return false;
  }

  cube->width = buf.width();
  cube->height = buf.height();
  cube->depth = buf.depth();

  cube->x_offset = buf.x_offset();
  cube->y_offset = buf.y_offset();
  cube->z_offset = buf.z_offset();

  switch (buf.data_type()) {
  case CubeletBuffer_DataType_UINT8:
    cube->datatype = Cubelet::CUBELET_UINT8;
    break;
  case CubeletBuffer_DataType_FLOAT32:
    cube->datatype = Cubelet::CUBELET_FLOAT32;
    break;
  }

  cube->data_file_offset = 0;
  cube->data = NULL;

  dataSizeBytes = buf.byte_count();
  dataHasBeenRead = false;

  return true;
}


void *CubeletStreamReader::getData(void *data) {

  if (dataSizeBytes == 0) return NULL;
  
  if (!data)
    data = malloc(dataSizeBytes);

  // the data can only be read once for each cubelet to avoid buffering
  // or seeking backwards
  if (dataHasBeenRead) {
    fprintf(stderr, "Error: cannot read the data for one cubelet twice.\n");
    return NULL;
  }

  if (fread(data, 1, dataSizeBytes, inf) != dataSizeBytes) {
    fprintf(stderr, "Failed to read cubelet content data.\n");
    return NULL;
  }

  dataHasBeenRead = true;

  return data;
}

    
void CubeletStreamReader::close() {
  if (inf != stdout) fclose(inf);
  inf = NULL;
}
