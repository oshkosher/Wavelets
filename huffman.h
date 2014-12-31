#ifndef __HUFFMAN_H__
#define __HUFFMAN_H__

#include <cstring>
#include <cassert>
#include <vector>
#include <string>
#include "bit_stack.h"
#include "bit_stream.h"


class Huffman {
 public:

  // Initialize a histogram with 'size' entries
  Huffman(int size_) : size(size_) {
    assert(counts.size() == 0);
    counts.resize(size, 0);
    computed = false;
  }

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
  int getLongestUsedEncodingLength() {return longestUsedBitString;}

  void encode(int value, int &length, unsigned bits[]);
  void decode(unsigned bits[], int &value);

  void encodeToStream(int value, BitStreamWriter *bitStream);
  int decodeFromStream(BitStreamReader *bitStream);

  struct Node {
    int count, value;
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

  // Table used to decode bits
  // This is an array of pairs of integers:
  // Offset 0:  [left  ] [right ]
  // Offset 2:  [left  ] [right ]
  // Offset 4:  [left  ] [right ]
  //   ...
  // If the left entry is -1, then this is a leaf node from the tree, and
  // the right entry contains the value.
  // Otherwise, if the current bit is 0, use the left entry to get the
  // offset of the next pair of entries, otherwise the bit is 1, and
  // use the right entry to get the next offset.
  int *decodeTable;

  /*
    Table used to encode data.
    The first size*2 entries are pairs of entries for each value:
       encodeTable[i*2] = # of bits used to encode value i
       encodeTable[i*2+1] = index in encodeTable of the first word
                            of bits enocoding value i

    For example, if the encoding for 3 takes 70 bits:
       encodeTable[6] = 70
       encodeTable[7] = 41
       ...
       encodeTable[41] = bits 0..31
       encodeTable[42] = bits 32..63
       encodeTable[43] = bits 64..69
  */
  std::vector<int> encodeTable;

  int longestBitString;  // length of the longest encoding
  int longestUsedBitString;  // length of the longest encoding of a value
                             // that appeared in the input data

  void buildTables(Node *root);

  // Put the nodes (leaf and internal) in order
  void orderNodes(Node *root, std::vector<Node*> &nodeTable,
		  std::vector<Node*> &valueTable);

  // fill in encodeTable
  void computeEncodeTable(Node *root,
			  const std::vector<Node*> &valueTable);

  // traverse the tree, filling in Node::encoding for each node.
  void computeEncodedBits(Node *node, BitStack &bits);
};


#endif // __HUFFMAN_H__
