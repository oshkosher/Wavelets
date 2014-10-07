#ifndef __HUFFMAN_H__
#define __HUFFMAN_H__

#include <cstring>
#include <cassert>


class Huffman {
 public:

  // Initialize a histogram with 'size' entries
  Huffman(int size_) : size(size_) {
    counts = new int[size];
    memset(counts, 0, sizeof(int) * size);
    computed = false;
    nodes = 0;
  }

  ~Huffman() {
    delete[] counts;
    if (nodes) delete[] nodes;
  }

  // increase the count for one entry by one
  void increment(int value) {
    assert(value >= 0 && value < size);
    counts[value]++;
  }

  // set the count for one entry
  void update(int value, int offset) {
    assert(value >= 0 && value < size);
    counts[value] += offset;
  }

  // get the count for one entry
  int getCount(int value) {
    assert(value >= 0 && value < size);
    return counts[value];
  }

  // computing the encoding
  void computeHuffmanCoding();

  // return the length of the encoded form of this value
  int encodedLength(int value);

  // 
  void encode(int value, int &length, int &bits);
  void decode(int bits, int &value);


  struct Node {
    int count, value, tableOffset;
    Node *left, *right, *parent;

    Node(unsigned value_, unsigned count_) {
      value = value_;
      count = count_;
      left = right = parent = NULL;
    }

    Node(Node *a, Node *b) {
      value = -1;
      count = 0;
      if (a) count += a->count;
      if (b) count += b->count;
      left = a;
      right = b;
      parent = NULL;
      left->parent = right->parent = this;
    }
  };

 private:
  // number of possible input values. Range will be 0..(size-1)
  int size;

  // array of 'size' entries, each one representing a frequency count
  int *counts;

  // Node in the tree for each value
  Node **nodes;

  // computeHuffmanCoding has been called
  bool computed;
};


#endif // __HUFFMAN_H__
