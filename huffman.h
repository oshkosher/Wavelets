#ifndef __HUFFMAN_H__
#define __HUFFMAN_H__

#include <cstring>
#include <cassert>
#include <vector>
#include <string>
#include "bit_stack.h"
#include "bit_stream.h"


class HuffmanDecoder {

  // Table used to decode bits
  // This is an array of pairs of integers:
  // Offset 0:  [left  ] [right ]
  // Offset 2:  [left  ] [right ]
  // Offset 4:  [left  ] [right ]
  //   ...
  // If the left entry is 0, then this is a leaf node from the tree, and
  // the right entry contains the value.
  // Otherwise, if the current bit is 0, use the left entry to get the
  // offset of the next pair of entries, otherwise the bit is 1, and
  // use the right entry to get the next offset.
  std::vector<int> table;

 public:
  void init(const std::vector<int> &table_) {
    table = table_;
  }

  // Read up to 'count' integers from bitStream. Return the number read.
  template<class WordSource>
  int decodeFromStream(int *values, int count,
                       BitStreamReader<WordSource> *bitStream) {
    // printf("decode ");
    for (int i=0; i < count; i++) {
      int pos = 0;
      if (bitStream->isEmpty()) return i;
      while (true) {
        unsigned bit = bitStream->readBit();
        // putchar(bit ? '1' : '0');
        pos = table[pos + bit] * 2;
        if (table[pos] == 0) {
          values[i] = table[pos+1];
          // printf("%d %d\n", i, values[i]);
          break;
        }
      }
    }    
  
    return count;
  }

  void getDecoderTable(std::vector<int> &table_) {
    table_ = table;
  }

};


class Huffman {
 public:

  Huffman() {
    init(0);
  }

  // Initialize a histogram with 'size' entries
  Huffman(int size_) {
    init(size_);
  }

  void init(int size_) {
    size = size_;
    assert(counts.size() == 0);
    counts.resize(size, 0);
    computed = false;
    nodeStorage = NULL;
  }

  ~Huffman() {
    if (nodeStorage) delete[] nodeStorage;
  }

  int getSize() {return size;}

  // increase the count for one entry by one
  void increment(int value) {
    assert(value >= 0 && value < size);
    counts[value]++;
  }

  // change the count for one entry by adding 'change'
  void update(int value, int change) {
    assert(value >= 0 && value < size);
    counts[value] += change;
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

  int getLongestEncodingLength() {return longestBitString;}

  void encode(int value, int &length, unsigned bits[]);
  void decode(unsigned bits[], int &value);

  // Write values[0..count-1] encoded to bitStream.
  template<class Sink>
  void encodeToStream(BitStreamWriter<Sink> *bitStream,
                      const int *values, int count) {
    assert(computed);
    unsigned *encodedData = encodeTable.data();

    for (int i = 0; i < count; i++) {
      // printf("%d %d\n", i, values[i]);
      unsigned value = values[i];
      unsigned bitCount = encodedData[value*2];
      unsigned bitDataOffset = encodedData[value*2+1];
      bitStream->write(encodedData+bitDataOffset, bitCount);
    }
    bitStream->flush();
  }

  // Read up to 'count' integers from bitStream. Return the number read.
  template<class WordSource>
  int decodeFromStream(int *values, int count,
                       BitStreamReader<WordSource> *bitStream) {
    return decoder.decodeFromStream(values, count, bitStream);
  }

  // print the frequency and encoding for each value
  void printEncoding();

  // fill in the a table that can be used to initialize a decoder table
  void getDecoderTable(std::vector<int> &table) {
    decoder.getDecoderTable(table);
  }

  struct Node {
    int value, count;
    int nodeId;   // order in a breadth-first traversal of all nodes;

    Node *left, *right, *parent;

    // only set for leaf nodes
    int bitLength;
    std::vector<unsigned> encoding;

    Node() {
      count = value = nodeId = bitLength = -1;
      left = right = parent = NULL;
    }

    void initLeaf(unsigned value_, unsigned count_) {
      count = count_;
      value = value_;
    }

    void initInternalNode(Node *left_, Node *right_) {
      assert(left_);
      assert(right_);
      left = left_;
      right = right_;
      count = left->count + right->count;
      left->parent = right->parent = this;
    }

    bool isLeaf() {return value != -1;}

    // returns a string representation of the bits for this node
    std::string bitString();
  };

 private:
  // number of possible input values. Range will be 0..(size-1)
  int size;

  // array of 'size' entries, each one representing a frequency count
  std::vector<int> counts;

  // computeHuffmanCoding has been called
  bool computed;

  // Array of all the nodes. May contain unused entries.
  // The first 'size' nodes are leaf nodes, the rest are interior nodes.
  // Only computeHuffmanCoding() should use this; everything else should
  // use 'nodes'.
  Node *nodeStorage;

  // All used nodes
  // nodes.size() = 2 * values.size() - 1
  std::vector<Node*> nodes;

  Node *root;

  // All nodes with nonzero counts, in decreasing count order.
  // Consider this the 'owner' of all the Node* objects.
  // When this is destroyed, destroy the nodes.
  std::vector<Node*> values;

  HuffmanDecoder decoder;

  /*
    Table used to encode data.
    The first size*2 entries are pairs of entries for each value.
       encodeTable[i*2] = # of bits used to encode value i
       encodeTable[i*2+1] = index in encodeTable of the first word
                            of bits enocoding value i

       If the value's count is 0 these will be 0.

    For example, if the encoding for 3 takes 70 bits:
       encodeTable[6] = 70
       encodeTable[7] = 41
       ...
       encodeTable[41] = bits 0..31
       encodeTable[42] = bits 32..63
       encodeTable[43] = bits 64..69
  */
  std::vector<unsigned> encodeTable;

  int longestBitString;  // length of the longest encoding

  // Put the nodes (leaf and internal) in BFS order
  void orderNodes(std::vector<Node*> &orderedNodes);

  // fill in encodeTable
  void computeEncodeTable();

  // fill in decodeTable
  void computeDecodeTable();

  // traverse the tree, filling in Node::encoding for each node.
  void computeEncodedBits(Node *node, BitStack &bits);
};

#endif // __HUFFMAN_H__
