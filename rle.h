#ifndef __RLE_H__
#define __RLE_H__

#include <iostream>
#include <map>
// #include <unordered_map>

// Convert a 2-d grid of float values in the range 0..1 into integers
// in the range 0..255. Using an instance of the template class "Sink",
// call sink.data(x) on each value, and sink.end() after the last one.
template <class Sink>
class FloatMatrixToGrayscale {
  int width, height;
  float *data;
  Sink *sink;

 public:
 FloatMatrixToGrayscale(Sink *sink, int width_, int height_, float *data_)
    : width(width_), height(height_), data(data_) {
  }
  
  int convert(float f) {
    /*
    if (f <= 0) return 0;
    if (f >= 1) return 255;
    */
    return (int)((f * 255) + 0.5f);
  }

  int getSize() {return width * height;}
  int get(int offset) {return convert(data[offset]);}

  void scan() {

    float *p = data, *end = data + (width*height);
    for (; p != end; p++) {
      sink->data(convert(*p));
    }
    sink->end();
    /*
    int n = width*height;
    for (int i=0; i < n; i++) 
      sink.data(convert(data[i]));
    sink.end();
    */  
  }
};


class PrintValues {
 public:
  void data(int value) {
    std::cout << value << std::endl;
  }

  void end() {}
};

/**
   PairSink
     void data(int value, int length);
     void end();
*/
template <class PairSink>
class EncodeRunLength {
  int prev, length, outputCount;
  PairSink *sink;

 public:
  EncodeRunLength(PairSink *sink_) : prev(-1), length(0), outputCount(0),
    sink(sink_) {}

  void data(int value) {
    // printf("encode value %d\n", value);
    if (value == prev && length < getMaxLen()) {
      length++;
    } else {
      if (length > 0) {
        // printf("output %d, len %d\n", prev, length);
	// printf("%d\n", length);
        sink->data(prev, length);
	outputCount++;
      }
      prev = value;
      length = 1;
    }
  }

  void end() {
    sink->data(prev, length);
    outputCount++;
    sink->end();
  }

  static int getMaxLen() {return (1 << getBitSize()) - 1;}
  static int getBitSize() {return 8;}

  int getOutputCount() {return outputCount;}
};

    
class PrintRunLength {
 public:
  void data(int value, int length) {
    std::cout << value << "  *" << length << std::endl;
  }
};


// PairSink
class FrequencyCountPair {
  void add(int value) {
    /*
    // std::unordered_map<int,int>::iterator iter = counts.find(value);
    MapType::iterator iter = counts.find(value);
    if (iter == counts.end()) {
      counts[value] = 1;
    } else {
      iter->second++;
    }
    */

    if (value >= 0 && value <= 255) counts[value]++;
  }

 public:

  FrequencyCountPair() {
    for (int i=0; i < 256; i++) counts[i] = 0;
  }

  // value, frequency
  // typedef std::unordered_map<int,int> MapType;
  // typedef std::map<int,int> MapType;
  // MapType counts;

  int counts[256];

  void data(int value, int length) {
    add(value);
    add(length);
  }

  void end() {
  }
};
  


#endif // __RLE_H__
