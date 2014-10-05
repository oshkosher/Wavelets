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
#include "rle.h"
#include "param_string.h"
#include "test_compress_common.h"

using namespace std;


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
         "\n",
         DEFAULT_WAVELET_STEPS,
         DEFAULT_THRESHOLD_FRACTION,
         DEFAULT_QUANTIZE_BITS,
         quantAlgId2Name(DEFAULT_QUANTIZE_ALGORITHM)
         );
  exit(1);
}


bool parseOptions(int argc, char **argv, Options &opt, int &nextArg) {
  opt.doCompress = true;
  opt.waveletSteps = DEFAULT_WAVELET_STEPS;
  opt.thresholdFraction = DEFAULT_THRESHOLD_FRACTION;
  opt.quantizeBits = DEFAULT_QUANTIZE_BITS;
  opt.quantizeAlgorithm = DEFAULT_QUANTIZE_ALGORITHM;

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

    else if (!strcmp(arg, "-qalg")) {
      if (++nextArg >= argc) printHelp();
      opt.quantizeAlgorithm = quantAlgName2Id(argv[nextArg]);
      if (opt.quantizeAlgorithm < 0) {
        fprintf(stderr, "Invalid quantize algorithm: \"%s\"\n", argv[nextArg]);
        return false;
      }
    }

    else {
      fprintf(stderr, "Unrecognized option: \"%s\"\n", arg);
      return false;
    }
  }

  return true;
}


/**
   Format: See the source.
   width (int, 4 bytes)
   height (int, 4 bytes)
   wavelet steps (int, 4 bytes)
   quantize bits (int, 4 bytes)
   threshold value (float, 4 bytes)
   quantization algorithm enum value (int, 4 bytes) 
   maximum value before quantization (float, 4 bytes)
   all the quantized data (width*height*(# of bits) bits, rounded up to a word)
*/
bool writeQuantDataSimple(const char *filename, FileData &f) {

  FILE *outf = fopen(filename, "wb");
  if (!outf) {
    printf("Error writing to \"%s\"\n", filename);
    return false;
  }

  assert(sizeof f.quantizeAlgorithm == sizeof(int));

  // write parameters to the header
  fwrite(&f.width, sizeof(int), 1, outf);
  fwrite(&f.height, sizeof(int), 1, outf);
  fwrite(&f.waveletSteps, sizeof(int), 1, outf);
  fwrite(&f.quantizeBits, sizeof(int), 1, outf);
  fwrite(&f.threshold, sizeof(float), 1, outf);
  fwrite(&f.quantizeAlgorithm, sizeof(int), 1, outf);
  // XXX we will need something different for Lloyd's algorithm codebook
  fwrite(&f.quantMaxVal, sizeof(float), 1, outf);
    
  bool success = writeQuantData(filename, outf, f);

  fclose(outf);

  return success;
}

/**
   Read the data written by writeQuantData.
*/
bool readQuantDataSimple(const char *filename, FileData &f) {

  FILE *inf = fopen(filename, "rb");
  if (!inf) {
    printf("Cannot read \"%s\"\n", filename);
    return false;
  }

  assert(sizeof f.quantizeAlgorithm == sizeof(int));

  fread(&f.width, sizeof(int), 1, inf);
  fread(&f.height, sizeof(int), 1, inf);
  fread(&f.waveletSteps, sizeof(int), 1, inf);
  fread(&f.quantizeBits, sizeof(int), 1, inf);
  fread(&f.threshold, sizeof(float), 1, inf);
  fread(&f.quantizeAlgorithm, sizeof(int), 1, inf);
  fread(&f.quantMaxVal, sizeof(float), 1, inf);

  bool success = readQuantData(filename, inf, f);

  fclose(inf);

  return success;
}


bool writeQuantDataParamStrings(const char *filename, FileData &f) {

  FILE *outf = fopen(filename, "wb");
  if (!outf) {
    printf("Error writing to \"%s\"\n", filename);
    return false;
  }

  // write parameters to the header
  ParamString p;
  p.setInt("w", f.width);
  p.setInt("h", f.height);
  p.setInt("ws", f.waveletSteps);
  p.setInt("qb", f.quantizeBits);
  p.setFloat("th", f.threshold);
  p.set("qa", quantAlgId2Name(f.quantizeAlgorithm));
  if (f.quantizeAlgorithm == QUANT_ALG_UNIFORM ||
      f.quantizeAlgorithm == QUANT_ALG_LOG) {
    p.setFloat("max", f.quantMaxVal);
  } else {
    p.setFloatList("qbound", f.quantBinBoundaries);
    p.setFloatList("cb", f.quantBinValues);
  }

  p.writeParameters(outf);
    
  bool success = writeQuantData(filename, outf, f);

  fclose(outf);

  return success;
}


/**
   Format:
   16 bytes: File type identification string "SCU wavelet 1.0\n"
   4 bytes: The length of the header in bytes, not including the 20 bytes
            for this and the ID string, encoded as a little-endian int.
   variable length: Header, encoded as a Google Protocol Buffer
   variable length: Data. Layout can be determined by reading the header.
*/
    
bool writeQuantDataProtoBuf(const char *filename, FileData &f) {

  FILE *outf = fopen(filename, "wb");
  if (!outf) {
    printf("Error writing to \"%s\"\n", filename);
    return false;
  }

  fputs(FILE_ID_STRING, outf);

  // build the protobuf
  WaveletCompressedImage buf;

  buf.set_width(f.width);
  buf.set_height(f.height);
  buf.set_wavelet_transform_step_count(f.waveletSteps);
  buf.set_quantize_bits(f.quantizeBits);
  buf.set_threshold_value(f.threshold);
  buf.set_wavelet_algorithm(WaveletCompressedImage_WaveletAlgorithm_HAAR);
  buf.set_quantization_algorithm(quantAlgId2ProtoId(f.quantizeAlgorithm));

  if (f.quantizeAlgorithm == QUANT_ALG_UNIFORM ||
      f.quantizeAlgorithm == QUANT_ALG_LOG) {
    buf.set_quant_max_value(f.quantMaxVal);
  } else {
    for (size_t i=0; i < f.quantBinBoundaries.size(); i++)
      buf.add_quant_bin_boundaries(f.quantBinBoundaries[i]);
    for (size_t i=0; i < f.quantBinValues.size(); i++)
      buf.add_quant_bin_values(f.quantBinValues[i]);
  }

  assert(sizeof(unsigned) == 4);
  unsigned codedLen = (unsigned) buf.ByteSize();
  fwrite(&codedLen, sizeof codedLen, 1, outf);
  printf("Header protobuf = %u bytes\n", codedLen);

  char *codedBuf = new char[codedLen];
  assert(codedBuf);
  if (!buf.SerializeToArray(codedBuf, codedLen)) {
    fprintf(stderr, "Failed to encode parameters\n");
    return false;
  }

  if (fwrite(codedBuf, 1, codedLen, outf) != codedLen) {
    fprintf(stderr, "Failed to write encoded parameters to file\n");
    return false;
  }
  delete[] codedBuf;
    
  bool success = writeQuantData(filename, outf, f);

  fclose(outf);

  return success;
}


bool readQuantDataParamStrings(const char *filename, FileData &f) {

  FILE *inf = fopen(filename, "rb");
  if (!inf) {
    printf("Cannot read \"%s\"\n", filename);
    return false;
  }

  ParamString p;
  p.readParameters(inf);
  if (!p.getInt("w", f.width)) printf("width not found\n");
  if (!p.getInt("h", f.height)) printf("height not found\n");
  if (!p.getInt("ws", f.waveletSteps)) printf("wavelet steps not found\n");
  if (!p.getInt("qb", f.quantizeBits)) printf("quant bits not found\n");
  if (!p.getFloat("th", f.threshold)) printf("threshold not found\n");
  string qa;
  if (!p.get("qa", qa)) printf("quant alg not found\n");
  f.quantizeAlgorithm = quantAlgName2Id(qa.c_str());
  if (f.quantizeAlgorithm == QUANT_ALG_UNIFORM ||
      f.quantizeAlgorithm == QUANT_ALG_LOG) {
    if (!p.getFloat("max", f.quantMaxVal)) printf("max value not found\n");
  } else {
    if (!p.getFloatList("cb", f.quantBinValues)) {
      printf("codebook not found\n");
      return false;
    }

    // check the codebook size
    if (f.quantBinValues.size() != (size_t)(1 << f.quantizeBits)) {
      fprintf(stderr, "Error: quantization of %d bits, but codebook has %d "
	      "entries (should be %d)\n",
	      f.quantizeBits, (int)f.quantBinValues.size(), 
	      (1 << f.quantizeBits));
      return false;
    }

  }

  bool success = readQuantData(filename, inf, f);

  fclose(inf);

  return success;
}


bool readQuantDataProtoBuf(const char *filename, FileData &f) {

  FILE *inf = fopen(filename, "rb");
  if (!inf) {
    printf("Cannot read \"%s\"\n", filename);
    return false;
  }

  char idString[17] = {0};
  fread(idString, 1, 16, inf);
  if (strcmp(idString, FILE_ID_STRING)) {
    fprintf(stderr, "Invalid file format. Expected protocol buffer header.\n");
    return false;
  }

  // read the header data
  unsigned codedLen;
  fread(&codedLen, sizeof codedLen, 1, inf);
  
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

  f.width = buf.width();
  f.height = buf.height();
  f.waveletSteps = buf.wavelet_transform_step_count();
  f.quantizeBits = buf.quantize_bits();
  f.threshold = buf.threshold_value();
  f.quantMaxVal = buf.quant_max_value();
  // buf.wavelet_algorithm();
  f.quantizeAlgorithm = quantProtoId2AlgId(buf.quantization_algorithm());

  f.quantBinBoundaries.resize(buf.quant_bin_boundaries_size());
  for (int i=0; i < buf.quant_bin_boundaries_size(); i++)
    f.quantBinBoundaries[i] = buf.quant_bin_boundaries(i);

  f.quantBinValues.resize(buf.quant_bin_values_size());
  for (int i=0; i < buf.quant_bin_values_size(); i++)
    f.quantBinValues[i] = buf.quant_bin_values(i);

  bool success = readQuantData(filename, inf, f);
  
  fclose(inf);

  return success;
}


bool readQuantData(const char *filename, FILE *inf, FileData &f) {
  int count = f.width * f.height;
  f.data = new float[count];

  BitStreamReader bits(inf);
  float *writePos = f.data;
  float *endPos = f.data + count;

  while (writePos < endPos) {
    if (bits.isEmpty()) {
      printf("Ran out of data reading %s, expected %d entries, got %d\n",
             filename, f.width * f.height, (int)(writePos - f.data));
      return false;
    }

    // read sign bit, value, length
    int sign = bits.read(1);
    // map 1,0 -> 2,0 -> -2,0 -> -1,1
    sign = 1 - (sign * 2);

    int quantized = bits.read(f.quantizeBits);
    int value = sign * quantized;

    int length = bits.read(8);
    
    for (int i=0; i < length; i++)
      *writePos++ = (float)value;
  }
  return true;
}


bool writeQuantData(const char *filename, FILE *outf, FileData &f) {

  // write the data
  BitStreamWriter bits(outf);
  int count = f.width*f.height;

  // this object takes (length,value) run-length pairs and writes
  // them in binary to the given bit stream
  WriteRLEPairsToBitStream rleToBits(&bits, f.quantizeBits);

  // this object takes data values as input, and passes (length,value)
  // run-length pairs to the rleToBits object
  EncodeRunLength<WriteRLEPairsToBitStream> rleEncoder(&rleToBits);

  // XXX see if it's faster to use pointer to traverse f.data[]

  for (int i=0; i < count; i++) {
    int x = (int) f.data[i];
    rleEncoder.data(x);
  }
  rleEncoder.end();

  printf("%d input values, %d RLE output pairs\n", count,
  rleEncoder.getOutputCount());
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
