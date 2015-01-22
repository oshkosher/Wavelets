/**
   Test the full path of tools:
     Read data file
     Perform wavelet transform
     Apply cutoff threshold
     Quantize
     Run-length encode
     Save to file

  And the inverse.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "nixtimer.h"
#include "rle.h"
#include "param_string.h"
#include "test_compress_common.h"
#include "huffman.h"
#include "quant.h"

using namespace std;
using namespace scu_wavelet;


// read&write just the data part of the file to&from f.intData
// using Huffman encoding
static void initHuffman(Huffman &huff, const CubeInt *cube, bool quiet);

static bool writeQuantDataHuffman(Huffman &huff, vector<uint32_t> *outData,
                                  const CubeInt *data, bool quiet);


void printHelp() {
  printf("\n"
         "  test_compress [opts] <input file> <output file>\n"
         "  Options:\n"
	 "    -d : decompress (default is to compress)\n"
         "    -steps <n> : # of wavelet transform steps (default = %d)\n"
         "    -thresh <n> : proportion of data to be eliminated in threshold cutoff\n"
         "                  Must be between [0..1]. (default = %.3f)\n"
         "    -qcount <n> : # of bins into which data is quantized.\n"
         "                 Must be >= 1 (default=%d)\n"
         "    -wave <wavelet> : Wavelet to use: haar or cdf97. (default = %s)\n"
         "    -qalg <alorithm> : quantization algorithm: uniform, log, count, or lloyd\n"
         "                       (default = %s)\n"
         "    -bq <filename> : before quantizing, save a copy of the data this file\n"
         "    -enc : print the bit encoding of each value\n"
         "    -nonstd : use nonstandard wavelet transpose order\n"
         "    -err : compute error metrics (slow, disabled by default)\n"
         "    -q : be quiet; suppess all output\n"
         "    -v : be verbose; print the data after each step\n"
         "\n",
         DEFAULT_WAVELET_STEPS,
         DEFAULT_THRESHOLD_FRACTION,
         DEFAULT_QUANTIZE_BINS,
         waveletAlgToName(DEFAULT_WAVELET),
         quantAlgId2Name(DEFAULT_QUANTIZE_ALGORITHM)
         );

  // defaults are defined in wavelet.h, just after enum declarations

  exit(1);
}


bool parseOptions(int argc, char **argv, Options &opt, int &nextArg) {
  opt.init();

  for (nextArg = 1; nextArg < argc; nextArg++) {
    const char *arg = argv[nextArg];
    if (arg[0] != '-') break;

    if (!strcmp(arg, "-d")) {
      opt.doCompress = false;
    }

    else if (!strcmp(arg, "-steps")) {
      if (++nextArg >= argc) printHelp();
      arg = argv[nextArg];
      int steps;
      if (1 != sscanf(arg, "%d", &steps) || steps < 0) {
        fprintf(stderr, "Invalid # of wavelet transform steps: \"%s\"\n", arg);
        return false;
      } else {
        opt.param.transformSteps = int3(steps, steps, steps);
      }
    }

    else if (!strcmp(arg, "-thresh")) {
      if (++nextArg >= argc) printHelp();
      arg = argv[nextArg];
      if (1 != sscanf(arg, "%f", &opt.param.thresholdFraction) ||
          opt.param.thresholdFraction < 0 ||
          opt.param.thresholdFraction > 1) {
        fprintf(stderr, "Invalid threshold proportion: \"%s\"\n", arg);
        return false;
      }
    }

    else if (!strcmp(arg, "-qcount")) {
      if (++nextArg >= argc) printHelp();
      arg = argv[nextArg];
      if (1 != sscanf(arg, "%d", &opt.param.binCount) ||
          opt.param.binCount < 1) {
        fprintf(stderr, "Invalid # quantize bins: \"%s\"\n", arg);
        return false;
      }
    }

    else if (!strcmp(arg, "-bq")) {
      if (++nextArg >= argc) printHelp();
      opt.saveBeforeQuantizingFilename = argv[nextArg];
    }

    else if (!strcmp(arg, "-qalg")) {
      if (++nextArg >= argc) printHelp();
      opt.param.quantAlg = quantAlgName2Id(argv[nextArg]);
      if (opt.param.quantAlg == QUANT_ALG_UNKNOWN) {
        fprintf(stderr, "Invalid quantize algorithm: \"%s\"\n", argv[nextArg]);
        return false;
      }
    }

    else if (!strcmp(arg, "-wave")) {
      if (++nextArg >= argc) printHelp();
      opt.param.waveletAlg = waveletAlgNameToId(argv[nextArg]);
      if (opt.param.waveletAlg == WAVELET_UNKNOWN) {
        fprintf(stderr, "Invalid wavelet: \"%s\"\n", argv[nextArg]);
        return false;
      }
    }

    else if (!strcmp(arg, "-enc")) {
      opt.printHuffmanEncoding = true;
    }

    else if (!strcmp(arg, "-nonstd")) {
      opt.param.isWaveletTransposeStandard = false;
    }

    else if (!strcmp(arg, "-q")) {
      opt.quiet = true;
    }

    else if (!strcmp(arg, "-v")) {
      opt.verbose = true;
    }

    else if (!strcmp(arg, "-err")) {
      opt.doComputeError = true;
    }

    else if (!strcmp(arg, "-experiment")) {
      opt.runQuantizationExperiments = true;
      opt.quiet = true;
    }

    else {
      fprintf(stderr, "Unrecognized option: \"%s\"\n", arg);
      return false;
    }
  }

  return true;
}


/**
   Format:
   16 bytes: File type identification string "SCU wavelet 1.0\n"
   4 bytes: The length of the header in bytes, not including the 20 bytes
            for this and the ID string, encoded as a little-endian int.
   variable length: Header, encoded as a Google Protocol Buffer
   variable length: Data. Layout can be determined by reading the header.
*/
#if 0
bool writeQuantData(const char *filename, const Data2d &data, Options &opt,
                    int *fileSizeBytes) {
  
  if (fileSizeBytes) *fileSizeBytes = 0;

  // initialize the huffman encoding
  Huffman huff;
  initHuffman(huff, data, opt.quantizeBits, opt.quiet);
  if (opt.printHuffmanEncoding) huff.printEncoding();
  
  FILE *outf = fopen(filename, "wb");
  if (!outf) {
    printf("Error writing to \"%s\"\n", filename);
    return false;
  }

  fputs(FILE_ID_STRING, outf);

  // build the protobuf
  WaveletCompressedImage buf;

  // copy parameters into protocol buffer
  buf.set_width(data.width);
  buf.set_height(data.height);
  buf.set_wavelet_transform_step_count(opt.waveletSteps.x);
  buf.set_standard_transpose(opt.isWaveletTransposeStandard);
  buf.set_quantize_bits(opt.quantizeBits);
  buf.set_threshold_value(opt.thresholdValue);
  buf.set_wavelet_algorithm(WaveletCompressedImage_WaveletAlgorithm_HAAR);
  buf.set_quantization_algorithm(quantAlgId2ProtoId(opt.quantizeAlgorithm));

  int sizeBeforeCodebook = buf.ByteSize();

  if (opt.quantizeAlgorithm == QUANT_ALG_UNIFORM ||
      opt.quantizeAlgorithm == QUANT_ALG_LOG) {
    buf.set_quant_max_value(opt.maxAbsVal);
  } else {
    for (size_t i=0; i < opt.quantBinBoundaries.size(); i++)
      buf.add_quant_bin_boundaries(opt.quantBinBoundaries[i]);
    for (size_t i=0; i < opt.quantBinValues.size(); i++)
      buf.add_quant_bin_values(opt.quantBinValues[i]);
  }

  int sizeBeforeHufftable = buf.ByteSize();
  int codebookSize = sizeBeforeHufftable - sizeBeforeCodebook;

  // add the Huffman decode table to the protocol buffer
  vector<int> huffDecodeTable;
  huff.getDecoderTable(huffDecodeTable);
  for (size_t i=0; i < huffDecodeTable.size(); i++) {
    buf.add_huffman_encode_table(huffDecodeTable[i]);
  }

  int huffDecodeTableSize = buf.ByteSize() - sizeBeforeHufftable;
  /*
  printf("Decode table, %d entries\n", (int)huffDecodeTable.size());
  for (size_t i=0; i < huffDecodeTable.size(); i += 2) {
    printf("%4d: %d %d\n", (int)i, huffDecodeTable[i], huffDecodeTable[i+1]);
  }
  */

  // get the size of the header and write the size to the file
  assert(sizeof(unsigned) == 4);
  unsigned codedLen = (unsigned) buf.ByteSize();
  fwrite(&codedLen, sizeof codedLen, 1, outf);
  if (!opt.quiet)
    printf("Header %u bytes (codebook %d bytes, %d bytes huff decode[%d])\n",
           codedLen, codebookSize, huffDecodeTableSize,
           (int)huffDecodeTable.size());;

  // encode the header
  char *codedBuf = new char[codedLen];
  assert(codedBuf);
  if (!buf.SerializeToArray(codedBuf, codedLen)) {
    fprintf(stderr, "Failed to encode parameters\n");
    return false;
  }

  // write the header to the file
  if (fwrite(codedBuf, 1, codedLen, outf) != codedLen) {
    fprintf(stderr, "Failed to write encoded parameters to file\n");
    return false;
  }
  delete[] codedBuf;
    
  bool success = writeQuantDataHuffman(huff, outf, data, opt.quantizeBits,
                                       opt.quiet);

  if (fileSizeBytes) *fileSizeBytes = (int) ftell(outf);
  
  fclose(outf);

  return success;
}
#endif


/**
   Write the data as one cubelet to the cubelet stream.
*/
bool writeQuantData(CubeletStreamWriter &cubeletStream,
                    CubeInt *cube, Options &opt,
                    int *sizeBytes) {
  
  if (sizeBytes) *sizeBytes = 0;

  // initialize the huffman encoding
  Huffman huff;
  initHuffman(huff, cube, opt.quiet);
  if (opt.printHuffmanEncoding) huff.printEncoding();

  vector<uint32_t> encodedData;

  if (!writeQuantDataHuffman(huff, &encodedData, cube, opt.quiet))
    return false;

  huff.getDecoderTable(cube->param.huffDecode);
  cube->param.compressedSize = encodedData.size() * sizeof(uint32_t);
  cube->isWaveletCompressed = true;

  // temporarily swap in the compressed data
  void *saveData = cube->data_;
  cube->data_ = encodedData.data();
  
  cubeletStream.addCubelet(cube);

  cube->data_ = saveData;

  // return success;
  return true;
}


bool readQuantData(CubeletStreamReader &cubeletStream, CubeInt *cube) {

  while (true) {
    if (!cubeletStream.next(cube)) {
      fprintf(stderr, "No compressed cubelet found\n");
      return false;
    }

    // look for a compressed cubelet in the right format
    if (cube->datatype == WAVELET_DATA_INT32 &&
        cube->isWaveletCompressed)
      break;
  }

  int wordCount = cube->param.compressedSize / sizeof(uint32_t);
  vector<uint32_t> bitData(wordCount, 0);

  // read the compressed data into an array of uint32_t's
  if (!cubeletStream.getRawData(bitData.data()))
    return false;

  // allocate storage for the decoded integer values
  cube->allocate();

  HuffmanDecoder huffDecoder;
  huffDecoder.init(cube->param.huffDecode);

  BitStreamMemorySource mem(&bitData);
  BitStreamReader<BitStreamMemorySource> bitReader(&mem);

  assert(cube->count());
  assert(cube->data());

  int readCount = huffDecoder.decodeFromStream
    (cube->data(), cube->count(), &bitReader);

  if (cube->count() != readCount) {
    printf("ERROR: read only %d of %d values\n", readCount, cube->count());
    return false;
  }

  return true;
}


#if 0
static void initHuffman(Huffman &huff, const Data2d &data, int quantizeBits,
                        bool quiet) {

  double startTime = NixTimer::time();

  int count = data.count();

  // number of possible values
  int valueCount = 1 << quantizeBits;
  huff.init(valueCount);

  assert(data.floatData == NULL && data.intData != NULL);

  // train the huffman encoder
  for (int i=0; i < count; i++) huff.increment(data.intData[i]);

  huff.computeHuffmanCoding();
  double elapsed = NixTimer::time() - startTime;
  if (!quiet)
    printf("Huffman build table %.3f ms\n", elapsed*1000);
}
#endif

class TraverseForHuffman {
public:
  Huffman &huff;
  TraverseForHuffman(Huffman &h) : huff(h) {}
  
  void visit(int value) {
    huff.increment(value);
  }
};

static void initHuffman(Huffman &huff, const CubeInt *cube, bool quiet) {

  double startTime = NixTimer::time();

  // number of possible values
  huff.init(cube->param.binCount);

  // train the huffman encoder
  TraverseForHuffman init(huff);
  cube->visit<TraverseForHuffman>(init);

  huff.computeHuffmanCoding();
  double elapsed = NixTimer::time() - startTime;
  if (!quiet)
    printf("Huffman build table %.3f ms\n", elapsed*1000);
}

#if 0
static bool writeQuantDataHuffman(Huffman &huff, FILE *outf,
                                  const Data2d &data, int quantizeBits,
                                  bool quiet) {

  // write the data
  BitStreamFileSink fileSink(outf);
  BitStreamWriter<BitStreamFileSink> bitWriter(&fileSink);
  int count = data.count();

  // number of possible values
  int valueCount = 1 << quantizeBits;
  
  huff.encodeToStream(&bitWriter, data.intData, count);
  bitWriter.flush();

  // printf("%llu bits written\n", (long long unsigned) bitWriter.size());
  size_t bitsWritten = bitWriter.size();
  int bytesWritten = (bitsWritten + 31) / 32 * 4;

  long long unsigned totalBits = 0;
  for (int i=0; i < valueCount; i++)
    totalBits += huff.encodedLength(i) * huff.getCount(i);

  if (!quiet)
    printf("Huffman encoding: %d bytes, %.2f bits/pixel, "
           "longest encoding = %d bits\n",
           bytesWritten, (double)totalBits / count,
           huff.getLongestEncodingLength());

  return true;
}
#endif

static bool writeQuantDataHuffman(Huffman &huff, vector<uint32_t> *outData,
                                  const CubeInt *data, bool quiet) {

  BitStreamMemorySink memorySink(outData);
  BitStreamWriter<BitStreamMemorySink> bitWriter(&memorySink);

  // this routine does not accept padded data
  assert(data->inset == int3(0,0,0));
  huff.encodeToStream(&bitWriter, data->pointer(0,0,0), data->count());

  // printf("%llu bits written\n", (long long unsigned) bitWriter.size());
  size_t bitsWritten = bitWriter.size();
  int bytesWritten = (bitsWritten + 31) / 32 * 4;

  long long unsigned totalBits = 0;
  for (int i=0; i < data->param.binCount; i++)
    totalBits += huff.encodedLength(i) * huff.getCount(i);

  if (!quiet)
    printf("Huffman encoding: %d bytes, %.2f bits/pixel, "
           "longest encoding = %d bits\n",
           bytesWritten, (double)totalBits / data->count(),
           huff.getLongestEncodingLength());

  return true;
}


Quantizer *createQuantizer(const WaveletCompressionParam &param) {

  int bits = ceilLog2(param.binCount);

  switch (param.quantAlg) {
  case QUANT_ALG_LOG:
    return new QuantLog(bits, param.thresholdValue, param.maxValue);
  case QUANT_ALG_UNIFORM:
    return new QuantUniform(bits, param.thresholdValue, param.maxValue);
  case QUANT_ALG_LLOYD:
    return new QuantCodebook(param.binBoundaries, param.binValues);
  default:
    fprintf(stderr, "Unknown quantization algorithm id: %d\n",
            (int)param.quantAlg);
    return NULL;
  }
}
