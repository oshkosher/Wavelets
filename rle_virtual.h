#ifndef __RLE_H__
#define __RLE_H__

#include <iostream>
#include <map>
#include <unordered_map>

class GrayscaleSink {
 public:
  virtual void data(int value) {}
  virtual void end() {}
};


class FloatMatrixToGrayscale {
  int width, height;
  float *data;
  GrayscaleSink &sink;
  
  int convert(float f) {
    if (f <= 0) return 0;
    if (f >= 1) return 255;
    return (int)((f * 255) + 0.5f);
  }

 public:
  FloatMatrixToGrayscale(int width_, int height_, float *data_,
                         GrayscaleSink &sink_)
    : width(width_), height(height_), data(data_), sink(sink_) {
  }

  int getSize() {return width * height;}
  int get(int offset) {return convert(data[offset]);}

  void scan() {
    float *p = data, *end = data + (width*height);
    for (; p != end; p++) {
      sink.data(convert(*p));
    }
    sink.end();
  }
};


class PrintValues : public GrayscaleSink {
 public:
  void data(int value) {
    std::cout << value << std::endl;
  }

  void end() {}
};


class PairSink {
 public:
  virtual void data(int value, int length) {}
};


class EncodeRunLength : public GrayscaleSink {
  int prev, length;
  PairSink &sink;

 public:
  EncodeRunLength(PairSink &sink_) : prev(-1), length(0), sink(sink_) {}

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
class FrequencyCountPair : public PairSink {
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
