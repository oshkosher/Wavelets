#ifndef __RLE_H__
#define __RLE_H__

#include <iostream>
#include <map>
#include <unordered_map>

template <class Sink>
class FloatMatrixToGrayscale {
  int width, height;
  float *data;
  
  int convert(float f) {
    /*
    if (f <= 0) return 0;
    if (f >= 1) return 255;
    */
    return (int)((f * 255) + 0.5f);
  }

 public:
  Sink sink;
  FloatMatrixToGrayscale(int width_, int height_, float *data_)
    : width(width_), height(height_), data(data_) {
  }

  int getSize() {return width * height;}
  int get(int offset) {return convert(data[offset]);}

  void scan() {

    float *p = data, *end = data + (width*height);
    for (; p != end; p++) {
      sink.data(convert(*p));
    }
    sink.end();
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


template <class PairSink>
class EncodeRunLength {
  int prev, length;

 public:
  PairSink sink;
  EncodeRunLength() : prev(-1), length(0) {}

  void data(int value) {
    if (value == prev && length < 255) {
      length++;
    } else {
      if (prev != -1) {
        sink.data(prev, length);
      }
      prev = value;
      length = 1;
    }
  }

  void end() {
    sink.data(prev, length);
  }
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

};
  


#endif // __RLE_H__
