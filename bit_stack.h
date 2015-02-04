#ifndef __BIT_STACK_H__
#define __BIT_STACK_H__

#include <vector>


class BitStack {
 private:
  std::vector<unsigned> vec;
  int bitCount, bitOffset;
  
 public:
  BitStack() : bitCount(0), bitOffset(0) {
    vec.push_back(0);
  }

  void push(int bit) {
    if (bitOffset == 32) {
      vec.push_back(0);
      bitOffset = 0;
    }

    unsigned onebit = 1 << bitOffset;

    // set or unset the given bit
    if (bit) {
      vec[vec.size()-1] |= onebit;
    } else {
      vec[vec.size()-1] &= ~onebit;
    }

    bitOffset++;
  }

  int pop() {
    if (bitOffset == 0) {
      if (vec.empty()) return 0;
      bitOffset = 32;
      vec.pop_back();
    }
    
    bitOffset--;

    unsigned mask = 1u << bitOffset;
    int bit = vec[vec.size()-1] & mask;

    // clear the bit
    vec[vec.size()-1] &= ~mask;

    // make sure it's just 1 or 0
    return bit != 0;
  }

  // include partial words
  int getWordCount() const {
    return (int) vec.size();
  }

  int getBitCount() const {
    return ((int)vec.size() - 1) * 32 + bitOffset;
  }

  unsigned getWord(int i) const {
    return vec[i];
  }
};
    
  

#endif // __BIT_STACK_H__
