#ifndef __HUFFMAN_H__
#define __HUFFMAN_H__

#include <cstring>
#include <cassert>
#include <vector>
#include <string>
#include "bit_stack.h"
#include "bit_stream.h"


// Logic for encoding long strings of duplicate values efficiently
class HuffmanDup {

 public:
  // use this many bits to encode the duplicate string length
  static const int DUP_ENCODE_BITS = 3;

  // offset the duplicate string length by this much
  static const int DUP_OFFSET = 3;

  static int decodeLength(unsigned v) {return 1 << (v + DUP_OFFSET);}
  static int dupMaximumValue() {return (1 << DUP_ENCODE_BITS) - 1;}
  static int dupMinimumLength() {return decodeLength(0);}
  static int dupMaximumLength() {return decodeLength(dupMaximumValue());}

  static const int DUP_MIN_ENCODED_LEN = 1 << DUP_OFFSET;

  // Returns encoding for the longest possible string less than or equal
  // to the given length. If length < dupMinimumLength(), returns 0.
  static unsigned dupFindLongest(int length) {

    if (length >= decodeLength(dupMaximumValue()))
      return dupMaximumValue();

    unsigned encoding = 0, encodedLength = dupMinimumLength();
    while (length >= (int)(encodedLength << 1)) {
      encoding++;
      encodedLength <<= 1;
    }
    return encoding;
  }

  // Given a length, count how many times the duplicate-string encoding
  // will be used, and how many time the regular single-encoding will be used
  static void dupCount(int length, int &singleEncodingCount,
		       int &dupEncodingCount) {
    dupEncodingCount = 0;

    while (length > dupMinimumLength()) {
      unsigned encoding = dupFindLongest(length);
      length -= decodeLength(encoding);
      dupEncodingCount++;
    }

    singleEncodingCount = length;
  }
};


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

  // When the value 'dupKey' is found in the input stream the next
  // DUP_ENCODE_BITS bits represent the length of a duplicates
  // of the value 'dupValue'. Use Huffman::decodeLength() to decode the bits
  // into a length.
  int dupKey, dupValue;

 public:
  void init(const std::vector<int> &table_) {
    table = table_;
    dupKey = -1;
    dupValue = 0;
  }

  void setDuplicate(int dupKey_, int dupValue_) {
    dupKey = dupKey_;
    dupValue = dupValue_;
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

	  int v = table[pos+1];

	  // this value marks a duplicate string, decode the length
	  if (v == dupKey) {
	    unsigned encodedLength =
	      bitStream->read(HuffmanDup::DUP_ENCODE_BITS);
	    int decodedLength = HuffmanDup::decodeLength(encodedLength);
	    // printf("%d. %d*[%d]\n", i, dupValue, decodedLength);
	    while (decodedLength > 0) {
	      values[i++] = dupValue;
	      decodedLength--;
	    }
	    i--;
	  } else {
	    values[i] = v;
	    // printf("%d %d\n", i, values[i]);
	  }
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
    dupValue = -1;
    dupValueCount = 0;
  }

  void init(const int *counts_, int size_) {
    size = size_;
    counts.assign(counts_, counts_ + size);
    computed = false;
    nodeStorage = NULL;
    dupValue = -1;
    dupValueCount = 0;
  }
  
  ~Huffman() {
    if (nodeStorage) delete[] nodeStorage;
  }

  // Mark a value as a common one that will have many strings of duplicates.
  // Must be in the range [0..getSize()-1].  Any value outside that range
  // will disable the 'common value' feature.
  void setDuplicateValue(int v) {
    if (v < 0 || v >= size) v = -1;
    dupValue = v;
  }

  int getDuplicateValue() {
    return dupValue;
  }

  // the value used internally to represent the start of the dup string
  int getDuplicateKey() {
    return size;
  }

  bool isDuplicateValueSet() {return dupValue >= 0;}

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

  // a string of 'length' copies of 'dupValue' was found.
  // update counters appropriately
  void addDupString(int length) {
    int single, dup;
    HuffmanDup::dupCount(length, single, dup);
    // add one for each instance of a single copy of 'dupValue'
    update(dupValue, single);
    // add one for each time the special dup encoding would be used
    dupValueCount += dup;
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

  // with the given frequencies for each value, return the total number
  // number of bytes needed to encode all of them.
  int totalEncodedLengthBytes();

  int getLongestEncodingLength() {return longestBitString;}

  void encode(int value, int &length, unsigned bits[]);
  void decode(unsigned bits[], int &value);

  // Write values[0..count-1] encoded to bitStream.
  template<class Sink>
  void encodeToStream(BitStreamWriter<Sink> *bitStream,
                      const int *values, int count) {
    assert(computed);
    // unsigned *encodedData = encodeTable.data();

    for (int i = 0; i < count; i++) {
      // printf("%d %d\n", i, values[i]);
      int value = values[i];

      // if this is the designated common value, see if it is duplicated
      if (value == dupValue) {
	int duplicateLength = 1;
	while (duplicateLength+i < count && 
	       values[i+duplicateLength] == dupValue)
	  duplicateLength++;
	dupEncodeString(duplicateLength, bitStream);
	i += duplicateLength-1;
      } else {
	// otherwise just write the encoding for this value
	encodeValue(value, bitStream);
      }
    }
    bitStream->flush();
  }

  // encode one value and write it to the string
  template<class Sink>
  void encodeValue(unsigned value, BitStreamWriter<Sink> *bitStream) {
    unsigned *encodedData = encodeTable.data();
    unsigned bitCount = encodedData[value*2];
    unsigned bitDataOffset = encodedData[value*2+1];
    bitStream->write(encodedData+bitDataOffset, bitCount);
  }

  // Encode 'length' copies of 'dupValue' as efficiently as possible.
  template<class Sink>
  void dupEncodeString(int length, BitStreamWriter<Sink> *bitStream) {

    // use the efficient encoding to reduce the length as much as possible
    // For example, if length==205: 128, 64, 8, leaving a remainder of 5

    while (length > HuffmanDup::dupMinimumLength()) {
      unsigned encoding = HuffmanDup::dupFindLongest(length);
      length -= HuffmanDup::decodeLength(encoding);
      encodeValue(getDuplicateKey(), bitStream);
      bitStream->write(encoding, HuffmanDup::DUP_ENCODE_BITS);
    }

    // write single instances; the remainder
    while (length > 0) {
      encodeValue(dupValue, bitStream);
      length--;
    }

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

  // If -1, this is ignored. Otherwise, the caller can set this to
  // mark a value that will be very common in the data. Repeated
  // instances of it use an encoding optimized for long strings of
  // duplicates.
  int dupValue, dupValueCount;

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
