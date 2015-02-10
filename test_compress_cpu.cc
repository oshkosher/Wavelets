#include <google/protobuf/stubs/common.h>
#include <algorithm>
#include "test_compress_cpu.h"
#include "test_compress_common.h"
#include "dwt_cpu.h"
#include "nixtimer.h"
#include "thresh_cpu.h"
#include "quant.h"
#include "lloyds.h"
#include "bit_stream.h"
#include "huffman.h"
#include "wavelet.h"
#include "cubelet_file.h"

using namespace scu_wavelet;
using namespace std;

bool compressFile(const char *inputFile, const char *outputFile, Options &opt);
bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt);


void computeLloydQuantization(const float *inputData, int count,
                              int binCount, float minVal, float maxVal,
                              float thresholdValue,
			      std::vector<float> &quantBinBoundaries,
			      std::vector<float> &quantBinValues);

void quantizationExperiments(CubeFloat &data, Options &opt,
                             const float *sortedAbsData,
                             float maxAbsVal);

  

int main(int argc, char **argv) {

  Options opt;
  int nextArg;

  if (!parseOptions(argc, argv, opt, nextArg)) return 1;
  
  // set global variable to enable/disable status output
  QUIET = opt.quiet;

  if (argc - nextArg != 2) printHelp();

  const char *inputFile = argv[nextArg++];
  const char *outputFile = argv[nextArg++];
  bool result;

  if (opt.doCompress) {
    result = compressFile(inputFile, outputFile, opt);
  } else {
    result = decompressFile(inputFile, outputFile, opt);
  }

  // deallocate static protobuf data
  google::protobuf::ShutdownProtobufLibrary();

  if (result == false) return 1;

  return 0;
}


bool compressFile(const char *inputFile, const char *outputFile,
                  Options &opt) {

  WaveletCompressionParam param = opt.param;
  Cube inputData;
  CubeFloat data;
  CubeInt quantizedData;
  double firstStartTime, startTime, elapsed;
  CubeletStreamReader cubeletStream;

  firstStartTime = NixTimer::time();

  // read and pad the data
  if (!readData(inputFile, &inputData, &data, &cubeletStream, &param))
    return false;

  // perform the wavelet transformation
  if (!waveletTransform(data, param, false, opt.verbose)) return false;

  // save the intermediate data to a file before quantizing
  if (opt.saveBeforeQuantizingFilename != "") {
    const char *filename = opt.saveBeforeQuantizingFilename.c_str();
    CubeletStreamWriter writer;
    data.param = param;

    // if the data is 1xHxD, transpose it to HxDx1
    bool isTransposed = false;
    if (data.width() == 1) {
      data.transpose3dFwd();
      isTransposed = true;
    }
    if (!writer.open(filename) ||
        !writer.addCubelet(&data) ||
        !writer.close()) {
      fprintf(stderr, "Failed to write intermediate data file \"%s\".\n",
              filename);
    } else {
      if (!QUIET)
        printf("Write intermediate data file \"%s\"\n", filename);
    }
    if (isTransposed)
      data.transpose3dBack();
  }


  // find the threshold value by sorting
  // XXX quickselect will speed up this selection, but then we'll lose the
  // benefits of having the sorted data
  const int count = data.count();
  float maxVal, minVal, maxAbsVal, *sortedAbsData = NULL;
  int nonzeroCount;

  startTime = NixTimer::time();

  param.thresholdValue
    = thresh_cpu(&data, param.thresholdFraction, &nonzeroCount, &maxVal,
                 &minVal, &sortedAbsData);

  maxAbsVal = sortedAbsData[count-1];
  
  // call the parameter optimization routine
  if (opt.doOptimize) {
    assert(inputData.datatype == WAVELET_DATA_UINT8);
    OptimizationData optData((CubeByte*)&inputData, &data, 
                             count, sortedAbsData,
                             nonzeroCount,
                             minVal, maxVal, maxAbsVal, param.transformSteps,
                             param.waveletAlg);
    if (!optimizeParameters(&optData, &param.thresholdValue, &param.binCount))
      return true;
  }

  // nonzeroCount is the number of nonzero entries (those whose absolute
  // values are greater than thresholdValue) in sortedAbsData[].
  // Since sortedAbsData is sorted, those are all at the end of the array.
  // Set nonzeroData to point to the beginning of that nonzero data.
  // NOTE: nonzeroData points to data in sortedAbsData, don't deallocate both
  float *nonzeroData = sortedAbsData + count - nonzeroCount;

  elapsed = NixTimer::time() - startTime;
  if (!QUIET)
    printf("min = %f, max = %f, threshold = %.10g: %.2f ms\n", minVal, maxVal,
           param.thresholdValue, elapsed*1000);

  // don't write a data file; just run some experiments testing different
  // quantization settings
  /*
  if (opt.runQuantizationExperiments) {
    quantizationExperiments(data, opt, sortedAbsData, maxAbsVal);
    return true;
  }
  */

  if (!quantize(data, quantizedData, maxAbsVal, param, nonzeroData,
                nonzeroCount, minVal, maxVal)) return false;

  // compute error rates
  if (opt.doComputeError) {
    float meanSqErr, peakSNR;
    // quantizedData.print("quantized before computing error rates");
    computeErrorRates(&quantizedData, param, &inputData, &meanSqErr, &peakSNR);
    printf("Mean squared error: %.3f, peak SNR: %.3f\n", meanSqErr, peakSNR);
  }

  if (opt.verbose) quantizedData.print("After quantization");
  
  // deallocate sortedAbsData
  delete[] sortedAbsData;
  nonzeroData = NULL;
  sortedAbsData = NULL;

  // write the quantized data to a file
  // for now, make a new file just for this cubelet
  startTime = NixTimer::time();
  quantizedData.param = param;
  CubeletStreamWriter cubeletWriter;
  if (!cubeletWriter.open(outputFile)) return false;
  if (!writeQuantData(cubeletWriter, &quantizedData, opt))
    return false;
  cubeletWriter.close();

    
  if (!QUIET) {
    printf("Write data file: %.2f ms\n", (NixTimer::time() - startTime)*1000);
    printf("Total: %.2f ms\n", (NixTimer::time() - firstStartTime)*1000);
  }

  return true;
}


bool decompressFile(const char *inputFile, const char *outputFile,
		    Options &opt) {

  CubeInt quantizedData;
  CubeFloat data;
  Cube outputData;
  double firstStartTime, startTime, elapsed;

  firstStartTime = startTime = NixTimer::time();

  // read the compressed file
  CubeletStreamReader cubeletStream;
  if (!cubeletStream.open(inputFile)) return false;
  if (!readQuantData(cubeletStream, &quantizedData)) return false;

  WaveletCompressionParam &param = quantizedData.param;

  elapsed = NixTimer::time() - startTime;
  if (!QUIET)
    printf("Read %dx%dx%d data file: %.2f ms\n",
           quantizedData.size.x, quantizedData.size.y, quantizedData.size.z,
           elapsed * 1000);

  data.size = quantizedData.size;
  data.param = param;
  data.allocate();

  if (opt.verbose) quantizedData.print("Before dequantize");
  
  // de-quantize the data
  dequantize(quantizedData, data, param);

  // perform inverse wavelet transform
  if (!waveletTransform(data, param, true, opt.verbose)) return false;

  // change the data back to the original datatype, if necessary
  outputData.size = param.originalSize;
  outputData.datatype = param.originalDatatype;
  translateCubeDataToOriginal(&data, &outputData, opt.verbose);

  // write the reconstructed data
  startTime = NixTimer::time();
  CubeletStreamWriter outStream;
  if (!outStream.open(outputFile) ||
      !outStream.addCubelet(&outputData) ||
      !outStream.close())
    return false;

  if (!QUIET) {
    printf("Write file: %.2f ms\n", (NixTimer::time() - startTime) * 1000);
    printf("Total: %.2f ms\n", (NixTimer::time() - firstStartTime) * 1000);
  }

  return true;
}


bool quantize(const CubeFloat &data, CubeInt &quantizedData,
              float maxAbsVal, WaveletCompressionParam &param,
              const float *nonzeroData, int nonzeroCount,
              float minVal, float maxVal,
              float *quantErrorOut) {

  // for now, to make the processing easier, don't accept padded data
  assert(data.size == data.totalSize);
  
  // allocate memory for the quantized data if it isn't already
  quantizedData.size = quantizedData.totalSize = data.size;
  quantizedData.inset = int3(0,0,0);
  quantizedData.allocate();

  double startTime, elapsed;

  switch (param.quantAlg) {

  case QUANT_ALG_UNIFORM:
  case QUANT_ALG_LOG:
    param.maxValue = maxAbsVal;
    break;

  case QUANT_ALG_LLOYD: 
    computeLloydQuantization(nonzeroData, nonzeroCount, param.binCount,
                             minVal, maxVal, param.thresholdValue,
                             param.binBoundaries, param.binValues);
    break;

  default: {}
  }

  Quantizer *quant = createQuantizer(param);
  if (!quant) return false;

  startTime = NixTimer::time();
  quant->quantizeRow(data.pointer(0,0,0), quantizedData.pointer(0,0,0),
                     data.count());

  if (!QUIET) {
    elapsed = NixTimer::time() - startTime;
    printf("Quantize %s: %.2f ms\n", quantAlgId2Name(param.quantAlg),
           elapsed*1000);
  }

  delete quant;
  
  return true;
}


/*
  If N is odd, distribute N codebook entries like this:

                                      
                                      
      (N-1)/2       1       (N-1)/2
  ---negatives--+-------+--positives---
                ^   0   ^             
                |       |             
          -thresh       +thresh      

  Except for the maximum positive value, which has its own codebook
  entry, negative and positive entries will be mirrored.

  For example, if N = 8, thresh=1, and max=10, one possible
  set of positive codebook entries is: 1, 3, 7
  All codebook entries:
   -7, -3, 1, 0, 1, 3, 7, 10

  If N is even, add the larger of minVal or maxVal to the
  beginning or end, respectively.

*/
void computeLloydQuantization
(const float *inputData, int count, 
 int binCount, float minVal, float maxVal,
 float thresholdValue,
 std::vector<float> &quantBinBoundaries,
 std::vector<float> &quantBinValues) {

  // Make 'codebookSize' entries on either size of 0
  int codebookSize = (binCount-1) / 2;

  quantBinBoundaries.clear();
  quantBinValues.clear();

  // float *binValues = new float[codebookSize];

  // inputData is sorted, use it to get minVal and maxVal
  // Skip the last entry in inputData[] because it is often much larger than
  // the rest and skews the results
  float maxAbsVal = inputData[count-2];
  assert(maxAbsVal > 0);

  // if the smallest input value is zero, set minAbsVal to the smallest
  // nonzero value
  /*
  float minAbsVal = inputData[0];
  if (minAbsVal <= 0) {
    minAbsVal = *std::upper_bound(inputData, inputData+count, 0);
  }
  */

  vector<float> codebook;
  initialLloydCodebook(codebook, codebookSize, thresholdValue, maxAbsVal);
  
  /*
  printf("Before Lloyd\n");
  for (int i=0; i < codebookSize; i++) printf("%d) %f\n", i, codebook[i]);
  */

  // fine-tune the codebook and bin boundaries using Lloyd's algorithm.
  // This also applies the quantization to each value, writing the values
  // to quantizedData[].
  float dist, reldist;
  unsigned *quantizedData = new unsigned[count];
  float *tmpBinBoundaries = new float[codebookSize-1];
  double startTime = NixTimer::time();
  int lloydIters = 0;
  lloyd(inputData, count-1, codebook.data(), (int)codebook.size(),
        tmpBinBoundaries, dist, reldist, quantizedData,
        DEFAULT_LLOYD_STOP_CRITERIA, &lloydIters);
  double elapsed = NixTimer::time() - startTime;
  if (!QUIET)
    printf("CPU Lloyd %d iters, %.2f ms\n", lloydIters, elapsed*1000);
  
  delete[] quantizedData;
  delete[] tmpBinBoundaries;

  /*
  printf("After Lloyd\n");
  for (int i=0; i < codebookSize; i++) printf("%d) %f\n", i, codebook[i]);
  */

  setBinsFromCodebook(quantBinValues, quantBinBoundaries, binCount,
                      codebook, thresholdValue, minVal, maxVal);
  
}


/*
void quantizationExperiments(const Data2d &data, Options &opt,
                             const float *sortedAbsData,
                             float maxAbsVal) {
  
  Data2d quantizedData;
  Data3d data3d;
  Data3dInt quantizedData3d;

  // allocate memory once for the quantized data
  quantizedData.initInts(data.width, data.height);
  const int count = data.width * data.height;

  float quantErr;
  int fileSize;

  // print header row
  printf("Bits\t" "ThresholdPct\t" "Log10MeanSqErr\t" "CompressionRatio\n");

  // try different values for quantizeBits and thresholdFraction

  for (opt.quantizeBits = 4; opt.quantizeBits <= 12; opt.quantizeBits++) {

    for (float thresholdFraction = 0.15f; thresholdFraction < 0.95f;
         thresholdFraction += 0.05f) {
      
      int threshIdx = (int)(thresholdFraction * count);
      opt.thresholdValue = sortedAbsData[threshIdx];

      int nonzeroCount = count - threshIdx - 1;
      const float *nonzeroData = sortedAbsData + count - nonzeroCount;

      // quantize, saving the error rate
      if (!quantize(data, quantizedData, data3d, quantizedData3d,
                    maxAbsVal, opt,
                    nonzeroData, nonzeroCount, &quantErr)) break;

      // write to a dummy file, saving the file size
      if (!writeQuantData("test_compress_cpu.experiment.wz",
                          quantizedData, opt, &fileSize)) break;

      // check quantError just in case, to avoid calling log(0)
      if (quantErr <= 0) {
        quantErr = 0;
      } else {
        quantErr = log(quantErr) / log(10);
      }

      double compressionRatio = 100.0 * (count - fileSize) / count;

      // write a row of data
      printf("%d\t%.1f\t%.7g\t%f\n", opt.quantizeBits, thresholdFraction*100,
             quantErr, compressionRatio);
      fflush(stdout);
    }
  }
}
*/


/*

//  test ErrorAccumulator
int main() {
  ErrorAccumulator err;
  err.add(10, 12);
  err.add(20, 17);
  err.add(30, 35);
  err.add(40, 50);
  err.setMaxPossible(255);
  printf("L1 %f, L2 %f, MSE %f, PSNR %f\n", err.getL1Error(), err.getL2Error(),
         err.getMeanSquaredError(), err.getPeakSignalToNoiseRatio());
  return 0;
}

*/
