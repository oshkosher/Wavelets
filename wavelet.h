#ifndef __WAVELET_H__
#define __WAVELET_H__

#include <climits>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>
#include "wavelet_compress.pb.h"

/*
  This file defines:
    QuantizeAlgorithm - enum of quantization algorithms
    int3 - wrapper for three ints x, y, z
    WaveletCompressionParam - all the parameters used when compressing and
                              decompressing data via our wavelet-based
                              algorithm
    Cube - 3-d data storage
*/

typedef enum {
  QUANT_ALG_UNKNOWN = -1,
  QUANT_ALG_UNIFORM = 1,
  QUANT_ALG_LOG,
  QUANT_ALG_LLOYD
} QuantizeAlgorithm;


typedef enum {
  WAVELET_UNKNOWN = -1,
  WAVELET_HAAR = 1,
  WAVELET_CDF97,
  WAVELET_CDF53
} WaveletAlgorithm;


typedef enum {
  WAVELET_DATA_UNKNOWN = -1,
  WAVELET_DATA_UINT8 = 1,
  WAVELET_DATA_INT32,
  WAVELET_DATA_FLOAT32
} WaveletDataType;

#define DEFAULT_WAVELET_STEPS -1
#define DEFAULT_THRESHOLD_FRACTION .5
#define DEFAULT_QUANTIZE_BINS 256
#define DEFAULT_WAVELET WAVELET_HAAR
#define DEFAULT_QUANTIZE_ALGORITHM QUANT_ALG_LOG
#define DEFAULT_WAVELET_TRANSPOSE_STANDARD true
#define DEFAULT_PADDING_METHOD REFLECT

// convert between quantization algorithms and their names as string
QuantizeAlgorithm quantAlgName2Id(const char *name);
const char *quantAlgId2Name(QuantizeAlgorithm id);

// convert between my quantization alg enum and the protobuf enum
CompressionParametersBuffer_QuantizationAlgorithm quantAlgId2ProtoId
  (QuantizeAlgorithm id);
QuantizeAlgorithm quantProtoId2AlgId
  (CompressionParametersBuffer_QuantizationAlgorithm protoId);

// convert between my wavelet algorithm enum and the protobuf enum
CompressionParametersBuffer_WaveletAlgorithm waveletAlgToProtoId
  (WaveletAlgorithm id);
WaveletAlgorithm protoIdToWaveletAlg
  (CompressionParametersBuffer_WaveletAlgorithm protoId);
WaveletAlgorithm waveletAlgNameToId(const char *name);
const char *waveletAlgToName(WaveletAlgorithm id);

// convert between my wavelet data type enum and the protobuf enum
CubeletBuffer_DataType datatypeToProtoId(WaveletDataType id);
WaveletDataType protoIdToWaveletDatatype(CubeletBuffer_DataType id);
const char *waveletDataTypeName(WaveletDataType id);
int waveletDataTypeSize(WaveletDataType id);


// to prevent collision with int3 class defined in CUDA

namespace scu_wavelet {  

class int3 {
 public:
  int x, y, z;
  
  int3() : x(0), y(0), z(0) {}
  int3(int x_, int y_, int z_) : x(x_), y(y_), z(z_) {}

  int count() const {return x*y*z;}

  bool operator == (const int3 &other) const {
    return x == other.x && y == other.y && z == other.z;
  }

  bool operator != (const int3 &other) const {
    return !(*this == other);
  }

  bool operator < (const int3 &other) const {
    return x < other.x && y < other.y && z < other.z;
  }

  bool operator > (const int3 &other) const {
    return other < *this;
  }

  bool operator >= (const int3 &other) const {
    return !(*this < other);
  }

  bool operator <= (const int3 &other) const {
    return !(other < *this);
  }

  int3 operator + (const int3 &other) const {
    return int3(x + other.x, y + other.y, z + other.z);
  }

  int parse(const char *str) {
    x = y = z = 1;
    return sscanf(str, "%d,%d,%d", &x, &y, &z);
  }
};

} // namespace scu_wavelet


class WaveletCompressionParam {
 public:
  typedef scu_wavelet::int3 int3;

  // the cube might have data added; this is the original size of the data
  int3 originalSize;

  WaveletDataType originalDatatype;

  // Number of wavelet steps done along each axis. A zero value represents
  // no wavelet steps along that axis. A negative value as input means
  // to perform the maximum number. That number will be filled in before
  // these parameters are saved.
  int3 transformSteps;

  WaveletAlgorithm waveletAlg;

  // if true, then do all wavelet steps in the x direction before doing all
  //          steps in the y direction
  // if false, do one x step, then one y step, and repeat
  bool isWaveletTransposeStandard;

  // Fraction of the data that was rounded down to zero before
  // quantizing.  This is not needed to decompress the data, but it
  // was used to determine threshValue, and may be interesting to
  // refer back to.
  float thresholdFraction;
  
  // smallest absolute value in the data; everthing below this was
  // rounded down to zero.
  float thresholdValue;

  // number of quantization bins
  int binCount;

  // which quantization algorithm was used
  QuantizeAlgorithm quantAlg;

  // greatest absolute value in the data; used in UNIFORM and LOG
  // quantization algorithms
  float maxValue;
  
  // size of the huffman-encoded data in bytes
  int compressedSize;

  // codebook and codebook boundaries; used in LLOYD
  // quantization algorithms
  std::vector<float> binBoundaries;
  std::vector<float> binValues;

  // table used to decode Huffman-encoded bits. See huffman.h (HuffmanDecoder)
  // for an explanation
  std::vector<int> huffDecode;

  // set default values
  void init() {
    originalDatatype = WAVELET_DATA_UNKNOWN;
    transformSteps = int3(DEFAULT_WAVELET_STEPS, DEFAULT_WAVELET_STEPS,
                          DEFAULT_WAVELET_STEPS);
    waveletAlg = DEFAULT_WAVELET;
    isWaveletTransposeStandard = DEFAULT_WAVELET_TRANSPOSE_STANDARD;
    thresholdFraction = DEFAULT_THRESHOLD_FRACTION;
    thresholdValue = 0;
    maxValue = 0;
    binCount = DEFAULT_QUANTIZE_BINS;
    quantAlg = DEFAULT_QUANTIZE_ALGORITHM;
    binBoundaries.clear();
    binValues.clear();
    huffDecode.clear();
  }
};


/**
   This class defines all the metadata about a cubelet. It contains
   a (void*) pointer to the data, which may or may not be set.
   To access the data, cast a pointer to a Cube to a pointer to
   a CubeNum<>, or use one of the typedefs CubeFloat, CubeInt, or CubeByte.
   CubeNum<> defines routines for accessing the actual data.
*/
class Cube {
 public:
  typedef scu_wavelet::int3 int3;

  // the apparent size (does not include padding)
  // int width, height, depth;
  int3 size;

  // Apply this as an offset when accessing data.
  // For example, if inset = (3,4,5), then accessing (1,1,1) will return (4,5,6)
  // This can be used to define a Cube that is a sub-cube of another.
  int3 inset;

  // the actual size, including padding
  // this is needed when calculating data locations when there are offsets.
  // for example, if the actual data contains 10 elements per row, but this
  // view of it has an inset of 3 on either side, then the apparent width
  // is 4, but to get to the next row, we need to add 10 to the pointer.
  int3 totalSize;

  // Offset of this cubelet in a containing dataset.
  int3 parentOffset;

 private:
  char idString[50];
 public:
  const char *getId() {
    sprintf(idString, "(%d,%d,%d)", 
            parentOffset.x, parentOffset.y, parentOffset.z);
    return idString;
  }
  
  // If this references data in a file, this is the offset of the data
  // in that file.
  uint64_t dataFileOffset;

  WaveletDataType datatype;

  // If true, then this object "owns" the data, so if 'data' is non-NULL
  // when the destructor is called, delete[] it.
  // False by default, calling 'allocate' sets this to true.
  bool ownsData;

  bool isWaveletCompressed;

  // If this cubelet represents data that has been wavelet-compressed, this
  // struture will contain all the parameters needed to decompress it.
  WaveletCompressionParam param;

  void *data_;

  bool allocate() {
    if (ownsData) deallocate();
    assert(size > int3(0,0,0));
    assert(inset >= int3(0,0,0));
    assert(datatype != WAVELET_DATA_UNKNOWN);

    size_t datumSize = waveletDataTypeSize(datatype);

    // if total size is unset, fill it in
    if (totalSize == int3(0,0,0)) totalSize = size + inset;

    int64_t length = (int64_t) totalSize.x * totalSize.y * totalSize.z;

    // everything is easier when you can index with ints, and this is designed
    // to handle reasonably-sized in-memory objects, so check that the indexes
    // are small enough to fit in an int
    assert(length <= INT_MAX);

    size_t byteCount = (size_t)datumSize * length;
    data_ = malloc(byteCount);
    assert(data_ != NULL);

    ownsData = true;

    memset(data_, byteCount, 0);

    return data_ != NULL;
  }

  /*
  bool allocate_untyped() {
    int datumSize = waveletDataTypeSize(datatype);
    if (datumSize < 0) {
      fprintf(stderr, "Cannot allocate memory for unknown data type %d\n",
              datatype);
      assert(false);
    }
    return allocate(datumSize);
  }
  */

  void deallocate() {
    if (data_) free(data_);
    data_ = NULL;
  }

  int width() const {return size.x;}
  int height() const {return size.y;}
  int depth() const {return size.z;}

  // allocate memory with the current settings

  int offset(int x, int y, int z) const {
    return (((z+inset.z) * totalSize.y)
            + (y+inset.y)) * totalSize.x
            + (x+inset.x);
  }

  int count() const {return size.count();}
  int totalCount() const {return totalSize.count();}
  bool hasInsets() const {return !(totalSize == size);}

  int getSizeInBytes() const {
    if (isWaveletCompressed) {
      return param.compressedSize;
    } else {
      int pixelSize = (datatype == WAVELET_DATA_UINT8) ? 1 : 4;
      return pixelSize * count();
    }
  }

  void copyFromCubeletBuffer(const CubeletBuffer *buf);
  void copyToCubeletBuffer(CubeletBuffer *buf) const;

  Cube() {
    init();
  }

  void init() {
    size = inset = totalSize = parentOffset = int3(0,0,0);
    dataFileOffset = 0;
    datatype = WAVELET_DATA_UNKNOWN;
    ownsData = false;
    isWaveletCompressed = false;
    param.init();
    data_ = NULL;
  }

  // default - do nothing

  void setType() {
    printf("warning: default setType() called\n");
  }

  const char *printfFormat() {
    printf("warning: default printfFormat() called\n");
    return "%f ";
  }
    
};


/**
   Define routines for accessing cube data.
   To access a single element, use get(), set() ref(), or pointer().
   To efficiently scan through all the data, use visit().
*/
template<class NUM>
class CubeNum : public Cube {
 public:

  // set the datatype field appropriately
  void setType();
  const char *printfFormat();

  const NUM *data() const {return (NUM*) data_;}
  NUM *data() {return (NUM*) data_;}

  CubeNum() {
    setType();
  }

  ~CubeNum() {
    if (ownsData && data_) free(data_);
  }

  // get one value
  NUM get(int x, int y, int z) const {
    return *pointer(x, y, z);
  }

  // set one value
  void set(int x, int y, int z, NUM value) {
    *pointer(x, y, z) = value;
  }

  // get a pointer into the data
  NUM* pointer(int x, int y, int z) {
    return data() + offset(x, y, z);
  }

  // get a pointer into the data
  const NUM* pointer(int x, int y, int z) const {
    return data() + offset(x, y, z);
  }

  // return an assignable reference to an entry
  NUM& ref(int x, int y, int z) {
    return data()[offset(x, y, z)];
  }

  const NUM& ref(int x, int y, int z) const {
    return data()[offset(x, y, z)];
  }


  // copy raw data into this one, automatically padding if necessary
  void copyFrom(const NUM *srcData, int3 srcSize) {

    // copy one row at a time, in case the sizes are different or
    // or there are insets in my data

    const NUM *readp = srcData;

    for (int z = 0; z < srcSize.z; z++) {
      for (int y = 0; y < srcSize.y; y++) {
        NUM *writep = pointer(0, y, z);
        memcpy(writep, readp, sizeof(NUM) * srcSize.x);
        readp += srcSize.x;
      }
    }

  }


  // copy raw data into this one, automatically padding if necessary
  void copyFrom(const CubeNum<NUM> &other) {

    assert(other.size <= size);

    // copy one row at a time

    for (int z = 0; z < other.depth(); z++) {
      for (int y = 0; y < other.height(); y++) {
        NUM *writep = pointer(0, y, z);
        const NUM *readp = other.pointer(0, y, z);
        memcpy(writep, readp, sizeof(NUM) * other.width());
      }
    }

  }


  void copyTo(NUM *data, int width, int height, int depth);

  // Visit all the normal elements (not the pad data), and call the
  // given function on each element.
  template <class Sink>
  void visit(Sink &sink) const {

    for (int z = 0; z < depth(); z++) {
      for (int y = 0; y < height(); y++) {
        const NUM *p = pointer(0, y, z);
        const NUM *e = p + width();
        while (p < e) sink.visit(*p++);
      }
    }

    /*
    for (int z = 0; z < depth(); z++) {
      for (int y = 0; y < height(); y++) {
        for (int x = 0; x < width(); x++) {
          sink.visit(get(x, y, z));
        }
      }
    }
    */
  }


  /* Visit every row of data, calling
       sink.visitRow(NUM *data, int length, int y, int z)
     on each one.
  */
  template <class Sink>
  void visitRows(Sink &sink) {
    for (int z = 0; z < depth(); z++) {
      for (int y = 0; y < height(); y++) {
        sink.visitRow(pointer(0, y, z), width(), y, z);
      }
    }
  }    


  void transpose2d(CubeNum<NUM> &dest) const {

    // transpose x&y of each level, but do not change levels
    assert(width() == dest.height());
    assert(height() == dest.width());
    assert(depth() == dest.depth());

    for (int z=0; z < depth(); z++)
      for (int y=0; y < height(); y++)
        for (int x=0; x < width(); x++)
          dest.set(y, x, z, get(x, y, z));

  }


  void transpose2dFast(CubeNum<NUM> &dest) const {

    const int subX = 10, subY = 10;

    for (int z=0; z < depth(); z++)
      for (int majorY=0; majorY < height(); majorY += subY)
        for (int majorX=0; majorX < width(); majorX += subX)

          for (int y = 0; y < subY; y++)
            for (int x = 0; x < subX; x++)
              // dest.set(majorY+y, majorX+x, z, get(majorX+x, majorY+y, z));
              dest.set(majorX+x, majorY+y, z, get(majorY+y, majorX+x, z));
  }


  // transpose xyz -> yzx
  void transpose3dFwd() {
    assert(size == totalSize);
    CubeNum<NUM> tmp;
    tmp.totalSize = tmp.size = int3(size.y, size.z, size.x);
    tmp.allocate();

    for (int z=0; z < depth(); z++)
      for (int y=0; y < height(); y++)
        for (int x=0; x < width(); x++)
          tmp.set(y, z, x, get(x, y, z));

    deallocate();
    data_ = tmp.data_;
    tmp.data_ = NULL;
    size = tmp.size;
    totalSize = tmp.totalSize;
  }

  // transpose xyz -> zxy
  void transpose3dBack() {
    assert(size == totalSize);
    CubeNum<NUM> tmp;
    tmp.size = tmp.totalSize = int3(size.z, size.x, size.y);
    tmp.allocate();

    for (int z=0; z < depth(); z++)
      for (int y=0; y < height(); y++)
        for (int x=0; x < width(); x++)
          tmp.set(z, x, y, get(x, y, z));

    deallocate();
    data_ = tmp.data_;
    tmp.data_ = NULL;
    size = tmp.size;
    totalSize = tmp.totalSize;
  }    

  void print(const char *title = NULL) {
    if (title) printf("%s\n", title);
    const char *format = printfFormat();

    
    for (int x=0; x<width(); x++) {
      printf(format, get(x,0,0));
    }
    putchar('\n');

    /*    
    for (int z=0; z<depth(); z++) {
      printf("z=%d\n", z);
      for (int y=0; y<height(); y++) {
        for (int x=0; x<width(); x++) {
          printf(format, get(x,y,z));
        }
        putchar('\n');
      }
      putchar('\n');
    }
    */
  }

};

typedef CubeNum<float> CubeFloat;
template<> void CubeNum<float>::setType();
template<> const char *CubeNum<float>::printfFormat();

typedef CubeNum<int> CubeInt;
template<> void CubeNum<int>::setType();
template<> const char *CubeNum<int>::printfFormat();

typedef CubeNum<unsigned char> CubeByte;
template<> void CubeNum<unsigned char>::setType();
template<> const char *CubeNum<unsigned char>::printfFormat();


#endif // __WAVELET_H__
