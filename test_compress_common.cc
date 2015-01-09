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

using namespace std;


// read&write just the data part of the file to&from f.intData
// using Huffman encoding
static void initHuffman(Huffman &huff, const Data2d &data, int quantizeBits,
                        bool quiet);
static bool writeQuantDataHuffman(Huffman &huff, FILE *outf,
                                  const Data2d &data, int quantizeBits,
                                  bool quiet);
static bool readQuantDataHuffman(HuffmanDecoder &huff, FILE *inf,
                                 Data2d &data);


void printHelp() {
  printf("\n"
         "  test_compress [opts] <input file> <output file>\n"
         "  Options:\n"
	 "    -d : decompress (default is to compress)\n"
         "    -steps <n> : # of wavelet transform steps (default = %d)\n"
         "    -thresh <n> : proportion of data to be eliminated in threshold cutoff\n"
         "                  Must be between [0..1]. (default = %.3f)\n"
         "    -qbits <n> : # of bits into which data is quantized. Must be between \n"
         "                 1 and 32.  (default=%d)\n"
         "    -qalg <alorithm> : quantization algorithm: uniform, log, count, or lloyd\n"
         "                       (default = %s)\n"
         "    -bq <filename> : before quantizing, save a copy of the data this file\n"
         "    -enc : print the bit encoding of each value\n"
         "    -std : use standard wavelet transpose order\n"
         "    -q : be quiet; suppess all output\n"
         "\n",
         DEFAULT_WAVELET_STEPS,
         DEFAULT_THRESHOLD_FRACTION,
         DEFAULT_QUANTIZE_BITS,
         quantAlgId2Name(DEFAULT_QUANTIZE_ALGORITHM)
         );
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
      if (1 != sscanf(arg, "%d", &opt.waveletSteps) ||
          opt.waveletSteps < 0) {
        fprintf(stderr, "Invalid # of wavelet transform steps: \"%s\"\n", arg);
        return false;
      }
    }

    else if (!strcmp(arg, "-thresh")) {
      if (++nextArg >= argc) printHelp();
      arg = argv[nextArg];
      if (1 != sscanf(arg, "%f", &opt.thresholdFraction) ||
          opt.thresholdFraction < 0 ||
          opt.thresholdFraction > 1) {
        fprintf(stderr, "Invalid threshold proportion: \"%s\"\n", arg);
        return false;
      }
    }

    else if (!strcmp(arg, "-qbits")) {
      if (++nextArg >= argc) printHelp();
      arg = argv[nextArg];
      if (1 != sscanf(arg, "%d", &opt.quantizeBits) ||
          opt.quantizeBits < 1 ||
          opt.quantizeBits > 32) {
        fprintf(stderr, "Invalid # quantize bits: \"%s\"\n", arg);
        return false;
      }
    }

    else if (!strcmp(arg, "-bq")) {
      if (++nextArg >= argc) printHelp();
      opt.saveBeforeQuantizingFilename = argv[nextArg];
    }

    else if (!strcmp(arg, "-qalg")) {
      if (++nextArg >= argc) printHelp();
      opt.quantizeAlgorithm = quantAlgName2Id(argv[nextArg]);
      if (opt.quantizeAlgorithm < 0) {
        fprintf(stderr, "Invalid quantize algorithm: \"%s\"\n", argv[nextArg]);
        return false;
      }
    }

    else if (!strcmp(arg, "-enc")) {
      opt.printHuffmanEncoding = true;
    }

    else if (!strcmp(arg, "-std")) {
      opt.isWaveletTransposeStandard = true;
    }

    else if (!strcmp(arg, "-q")) {
      opt.quiet = true;
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
  buf.set_wavelet_transform_step_count(opt.waveletSteps);
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


bool readQuantData(const char *filename, Data2d &data, Options &opt) {

  FILE *inf = fopen(filename, "rb");
  if (!inf) {
    printf("Cannot read \"%s\"\n", filename);
    return false;
  }

  char idString[17] = {0};
  if (16 != fread(idString, 1, 16, inf)) return false;
  if (strcmp(idString, FILE_ID_STRING)) {
    fprintf(stderr, "Invalid file format. Expected protocol buffer header.\n");
    return false;
  }

  // read the header data
  unsigned codedLen;
  if (1 != fread(&codedLen, sizeof codedLen, 1, inf)) return false;
  
  char *codedBuf = new char[codedLen];
  assert(codedBuf);
  if (codedLen != (unsigned)fread(codedBuf, 1, codedLen, inf)) {
    fprintf(stderr, "Failed to read header.\n");
    return false;
  }

  // decode the protobuf
  WaveletCompressedImage buf;
  if (!buf.ParseFromArray(codedBuf, codedLen)) {
    fprintf(stderr, "Failed to decode header data.\n");
    return false;
  }
  delete[] codedBuf;

  // copy parameters from protocol buffer
  data.initInts(buf.width(), buf.height());
  opt.waveletSteps = buf.wavelet_transform_step_count();
  opt.isWaveletTransposeStandard = buf.standard_transpose();
  opt.quantizeBits = buf.quantize_bits();
  opt.thresholdValue = buf.threshold_value();
  opt.maxAbsVal = buf.quant_max_value();
  opt.quantizeAlgorithm = quantProtoId2AlgId(buf.quantization_algorithm());

  opt.quantBinBoundaries.resize(buf.quant_bin_boundaries_size());
  for (int i=0; i < buf.quant_bin_boundaries_size(); i++)
    opt.quantBinBoundaries[i] = buf.quant_bin_boundaries(i);

  opt.quantBinValues.resize(buf.quant_bin_values_size());
  for (int i=0; i < buf.quant_bin_values_size(); i++)
    opt.quantBinValues[i] = buf.quant_bin_values(i);

  // get the Huffman decode table from the protocol buffer
  vector<int> huffDecodeTable;
  for (int i=0; i < buf.huffman_encode_table_size(); i++)
    huffDecodeTable.push_back(buf.huffman_encode_table(i));

  HuffmanDecoder huffDecoder;
  huffDecoder.init(huffDecodeTable);

  bool success = readQuantDataHuffman(huffDecoder, inf, data);
  
  fclose(inf);

  return success;
}


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


static bool writeQuantDataHuffman(Huffman &huff, FILE *outf,
                                  const Data2d &data, int quantizeBits,
                                  bool quiet) {

  // write the data
  BitStreamWriter bitWriter(outf);
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


static bool readQuantDataHuffman(HuffmanDecoder &huff, FILE *inf,
                                 Data2d &data) {

  BitStreamReader bitReader(inf);

  int count = data.count();
  assert(data.intData && !data.floatData);

  int readCount = huff.decodeFromStream(data.intData, count, &bitReader);

  if (count != readCount) {
    printf("ERROR: read only %d of %d values\n", count, readCount);
    return false;
  }

  return true;
}


QuantizeAlgorithm quantAlgName2Id(const char *name) {
  if (!strcmp(name, "uniform")) return QUANT_ALG_UNIFORM;
  if (!strcmp(name, "log")) return QUANT_ALG_LOG;
  if (!strcmp(name, "lloyd")) return QUANT_ALG_LLOYD;
  if (!strcmp(name, "count")) return QUANT_ALG_COUNT;
  return QUANT_ALG_UNKNOWN;
}
  

const char *quantAlgId2Name(QuantizeAlgorithm id) {
  switch (id) {
  case QUANT_ALG_UNIFORM: return "uniform";
  case QUANT_ALG_LOG: return "log";
  case QUANT_ALG_COUNT: return "count";
  case QUANT_ALG_LLOYD: return "lloyd";
  default: return NULL;
  }
}
  

WaveletCompressedImage_QuantizationAlgorithm quantAlgId2ProtoId
  (QuantizeAlgorithm id) {

  switch (id) {
  case QUANT_ALG_UNIFORM:
    return WaveletCompressedImage_QuantizationAlgorithm_UNIFORM;
  case QUANT_ALG_LOG: 
    return WaveletCompressedImage_QuantizationAlgorithm_LOG;
  case QUANT_ALG_COUNT:
    return WaveletCompressedImage_QuantizationAlgorithm_COUNT;
  case QUANT_ALG_LLOYD:
    return WaveletCompressedImage_QuantizationAlgorithm_LLOYD;
  default:
    return WaveletCompressedImage_QuantizationAlgorithm_UNIFORM;
  }
}
  

QuantizeAlgorithm quantProtoId2AlgId
  (WaveletCompressedImage_QuantizationAlgorithm protoId) {

  switch (protoId) {
  case WaveletCompressedImage_QuantizationAlgorithm_UNIFORM:
    return QUANT_ALG_UNIFORM;
  case WaveletCompressedImage_QuantizationAlgorithm_LOG:
    return QUANT_ALG_LOG;
  case WaveletCompressedImage_QuantizationAlgorithm_COUNT:
    return QUANT_ALG_COUNT;
  case WaveletCompressedImage_QuantizationAlgorithm_LLOYD:
    return QUANT_ALG_LLOYD;
  default:
    return QUANT_ALG_UNIFORM;
  }
}

