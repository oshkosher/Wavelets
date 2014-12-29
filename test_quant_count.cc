#include <cstdio>
#include <cstring>
#include <vector>
#include <cmath>
#include "nixtimer.h"
#include "quant_count.h"
#include "quant_unif_cpu.h"
#include "quant_log_cpu.h"
#include "dquant_unif_cpu.h"
#include "dquant_log_cpu.h"
#include "quant.h"
#include "data_io.h"
#include "thresh_cpu.h"


int main2() {
  float data[100];
  int len = 100;
  for (int i=0; i < 100; i++) data[i] = i / 100.0f;

  std::vector<float> boundaries, codebook;

  quant_count_init_cpu(len, data, 2, .4, boundaries, codebook);
  for (int i=0; i < 3; i++)
    printf("  %f", boundaries[i]);
  printf("\n");
  for (int i=0; i < 4; i++)
    printf("  %f", codebook[i]);
  printf("\n");

  return 0;
}


int main3() {
  float *data;
  int width, height;
  
  if (!readDataFile("lenna.data", &data, &width, &height)) {
    printf("Failed to read lenna.data\n");
    return 1;
  }

  std::vector<float> boundaries, codebook;
  int bits = 3;

  quant_count_init_cpu(width * height, data, bits, .1, boundaries, codebook);
  int binCount = 1 << bits;
  for (int i=0; i < binCount-1; i++)
    printf("  %f", boundaries[i]);
  printf("\n");
  for (int i=0; i < binCount; i++)
    printf("  %f", codebook[i]);
  printf("\n");

  return 0;
}


int main() {
  float *data, *origData, *dequantizedData;
  int width, height, *quantizedData;
  
  if (!readDataFile("lenna.data", &origData, &width, &height)) {
    printf("Failed to read lenna.data\n");
    return 1;
  }

  // the old algorithms will quantize in-place, overwriting data[]
  // with floats that are integer values
  data = new float[width*height];

  // the new algorithms will quantize to this array
  quantizedData = new int[width*height];

  // and dequantize to this array
  dequantizedData = new float[width*height];


  std::vector<float> boundaries, codebook;
  int bits = 4, nonzeroCount;
  float maxVal, minVal;
  double startTime, elapsed1, elapsed2, elapsed3, err;
  float threshold = thresh_cpu(width*height, origData, .1, &nonzeroCount,
                               &maxVal, &minVal);

  bool doComputeErr = false;

  for (int algorithm=1; algorithm <= 3; algorithm++) {
    printf("Algorithm %d\n", algorithm);

    memcpy(data, origData, sizeof(float)*width*height);

    if (algorithm == 1) {
      startTime = NixTimer::time();
      quant_unif_cpu(width, data, bits, threshold, maxVal);
      elapsed1 = NixTimer::time() - startTime;
      QuantUniform qu;
      qu.init(bits, threshold, maxVal);
      QuantizationLooper<QuantUniform> qlu;
      qlu.init(&qu);
      qlu.quantize(width*height, origData, quantizedData, doComputeErr);
      elapsed2 = qlu.getExecuteTime();
      err = qlu.getError();
      qlu.dequantize(width*height, quantizedData, dequantizedData);
      elapsed3 = qlu.getExecuteTime();
    }

    else if (algorithm == 2) {
      startTime = NixTimer::time();
      quant_log_cpu(width, data, bits, threshold, maxVal);
      elapsed1 = NixTimer::time() - startTime;
      QuantLog qu;
      qu.init(bits, threshold, maxVal);
      QuantizationLooper<QuantLog> qlu;
      qlu.init(&qu);
      qlu.quantize(width*height, origData, quantizedData, doComputeErr);
      elapsed2 = qlu.getExecuteTime();
      err = qlu.getError();
      qlu.dequantize(width*height, quantizedData, dequantizedData);
      elapsed3 = qlu.getExecuteTime();
    }

    else {
      quant_count_init_cpu(width*height, origData, bits, threshold, boundaries,
			   codebook);
      startTime = NixTimer::time();
      quant_boundaries_array(boundaries, width*height, data);
      elapsed1 = NixTimer::time() - startTime;
      QuantCodebook qu(boundaries, codebook);
      QuantizationLooper<QuantCodebook> qlu;
      qlu.init(&qu);
      qlu.quantize(width*height, origData, quantizedData, doComputeErr);
      elapsed2 = qlu.getExecuteTime();
      err = qlu.getError();
      qlu.dequantize(width*height, quantizedData, dequantizedData);
      elapsed3 = qlu.getExecuteTime();
    }

    int mismatchCount = 0, matchCount = 0;
    for (int i=0; i < width*height; i++) {
      if ((int)data[i] != quantizedData[i]) {
	if (mismatchCount++ < 10) {
	  printf("Quant mismatch at %d: %.10g -> %.0f vs %d\n",
		 i, origData[i], data[i], quantizedData[i]);
	}
      } else {
	matchCount++;
      }
    }

    printf("  quant: %d matches, err=%g, %.3f ms old, %.3f ms new\n",
	   matchCount, err, elapsed1*1000, elapsed2*1000);

    if (algorithm <= 2) {

      startTime = NixTimer::time();
      if (algorithm == 1)
	dquant_unif_cpu(width, data, bits, threshold, maxVal);
      else if (algorithm == 2)
	dquant_log_cpu(width, data, bits, threshold, maxVal);
      elapsed1 = NixTimer::time() - startTime;

      int mismatchCount = 0, matchCount = 0;
      for (int i=0; i < width*height; i++) {
	if (fabsf(data[i] - dequantizedData[i]) > 0.000001) {
	  if (mismatchCount++ < 10) {
	    printf("Dequant mismatch at %d: %d -> %f vs %f\n",
		   i, quantizedData[i], data[i], dequantizedData[i]);
	  }
	} else {
	  matchCount++;
	}
      }

      printf("  dequant: %d matches, %.3f ms old, %.3f ms new\n",
	     matchCount, elapsed1*1000, elapsed3*1000);
    }
	     
  }
  
  delete[] data;
  delete[] origData;
  delete[] quantizedData;

  return 0;
}

#define LEN 15

int main4() {
  float data[LEN*LEN], orig[LEN], oldq[LEN*LEN], olddq[LEN], dq[LEN];
  int iq[LEN];
  int bits = 4, half = LEN/2;
  float threshold = 1.22f;


  orig[half] = 0;
  for (int i=1; i <= half; i++) {
    orig[half-i] = -pow(1.2, i);
    orig[half+i] = pow(1.2, i);
  }

  float maxVal = orig[LEN-1];

  memcpy(data, orig, sizeof orig);

  quant_log_cpu(LEN, data, bits, threshold, maxVal);
  memcpy(oldq, data, sizeof oldq);

  QuantLog q;
  q.init(bits, threshold, maxVal);
  QuantizationLooper<QuantLog> ql;
  ql.init(&q);
  ql.quantize(LEN, orig, iq);

  dquant_log_cpu(LEN, data, bits, threshold, maxVal);
  memcpy(olddq, data, sizeof olddq);
  ql.dequantize(LEN, iq, dq);

  for (int i=0; i < LEN; i++) {
    printf("%2d  %6.3f  %2d %2d  %6.3f %6.3f\n",
	   i, orig[i], (int)oldq[i], iq[i], olddq[i], dq[i]);
  }

  return 0;
}

  
