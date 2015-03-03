#include <queue>
#include <algorithm>
#include <cstdio>
#include "huffman.h"
#include "bit_stack.h"

using namespace std;

class CompareHuffmanNode {
public:
  bool operator() (const Huffman::Node *a, const Huffman::Node *b) const {
    if (a->count == b->count) {
      return a < b;
    } else {
      return a->count > b->count;
    }
  }
};

static void indent(int depth) {
  for (int i=0; i < depth*2; i++) putchar(' ');
}

static void printTree(Huffman::Node *node, int depth = 0) {
  if (node->value != -1) {
    printf("%d", node->value);
  } else {
    putchar('\n');
    indent(depth);
    printf("0: ");
    printTree(node->left, depth+1);
    putchar('\n');
    indent(depth);
    printf("1: ");
    printTree(node->right, depth+1);
  }

  if (depth==0) putchar('\n');
}
    

void Huffman::computeHuffmanCoding() {
  if (nodeStorage) delete[] nodeStorage;
  nodes.clear();
  values.clear();
  longestBitString = 0;
  assert((int)counts.size() == size);

  // temporarily add dupValueCount to counts[]
  if (isDuplicateValueSet())
    counts.push_back(dupValueCount);

  // allocate space for all the nodes
  if (nodeStorage) delete[] nodeStorage;
  nodeStorage = new Node[counts.size() * 2 - 1];

  // Initialize all the leaf nodes,
  // add all the values that were used to encodedNodes and to queue
  priority_queue<Node*, vector<Node*>, CompareHuffmanNode> queue;
  int nodeCount = counts.size();
  for (int i=0; i < nodeCount; i++) {
    Node *node = nodeStorage + i;
    node->initLeaf(i, counts[i]);
    if (counts[i] > 0) {
      nodes.push_back(node);
      values.push_back(node);
      queue.push(node);
    }
  }

  // remove dupValueCount
  if (isDuplicateValueSet()) counts.pop_back();

  // sort by decreasing frequency
  CompareHuffmanNode sortComparison;
  sort(values.begin(), values.end(), sortComparison);

  // for N elements, do N-1 iterations, merging the two smallest nodes
  int storageIdx = nodeCount;
  while (queue.size() > 1) {
    Node *a = queue.top(); queue.pop();
    Node *b = queue.top(); queue.pop();
    Node *mergeNode = nodeStorage + storageIdx++;
    mergeNode->initInternalNode(a, b);
    queue.push(mergeNode);
    nodes.push_back(mergeNode);
  }

  root = queue.top(); queue.pop();
  // printTree(root);

  computeEncodeTable();

  computeDecodeTable();

  computed = true;
}


void Huffman::printEncoding() {
  for (size_t i=0; i < values.size(); i++) {
    Node *node = values[i];
    string bits = node->bitString();
    printf("%3d (%8d): %d %s\n", node->value, node->count, node->bitLength,
	   bits.c_str());
  }
}


void Huffman::computeEncodeTable() {

  // traverse the tree, saving the encoded form of each leaf node
  BitStack bits;
  computeEncodedBits(root, bits);

  // build encodeTable
  encodeTable.clear();
  int valueCount = size + (isDuplicateValueSet() ? 1 : 0);
  encodeTable.resize(valueCount * 2, 0);
  for (size_t i=0; i < values.size(); i++) {
    Node *node = values[i];

    // string encoding = node->bitString();
    // printf("Encode %d (count=%d) as %s\n", node->value, node->count, encoding.c_str());

    // save the bit length and the offset of where the bits are stored
    encodeTable[node->value*2] = node->bitLength;
    encodeTable[node->value*2+1] = encodeTable.size();
    
    for (size_t wordNo=0; wordNo < node->encoding.size(); wordNo++)
      encodeTable.push_back(node->encoding[wordNo]);

  }

  /*
  printf("Encode table\n");
  for (int i=0; i < size*2; i += 2) {
    printf("%3d: %d bits at %d\n", i, encodeTable[i], encodeTable[i+1]);
  }
  for (size_t i=size*2; i < encodeTable.size(); i++) {
    printf("%3d: %u\n", (int)i, encodeTable[i]);
  }
  */
}
  

void Huffman::computeEncodedBits(Node *node, BitStack &bits) {
  if (node->isLeaf()) {

    // save the bits the encode this value
    node->bitLength = bits.getBitCount();
    node->encoding.resize(bits.getWordCount());
    for (int i=0; i < bits.getWordCount(); i++) {
      node->encoding[i] = bits.getWord(i);
    }

    // track the longest bit strings
    if (node->bitLength > longestBitString)
      longestBitString = node->bitLength;

  } else {

    assert(node->left);
    assert(node->right);

    bits.push(0);
    computeEncodedBits(node->left, bits);
    bits.pop();
    bits.push(1);
    computeEncodedBits(node->right, bits);
    bits.pop();

  }
}


void Huffman::computeDecodeTable() {
  vector<int> decodeTable;
  
  // Do a breadth-first traversal of the tree to order the nodes by
  // increasing encoding length, and set their nodeId to their position
  // in the ordered list.
  vector<Node*> treeOrderedNodes;
  orderNodes(treeOrderedNodes);


  // now build a decode table by converting the nodes and references
  // to their children into an array of pairs of indices.
  decodeTable.resize(2 * treeOrderedNodes.size());
  for (size_t i=0; i < treeOrderedNodes.size(); i++) {
    Node *node = treeOrderedNodes[i];
    if (node->isLeaf()) {
      decodeTable[2*i] = 0;
      decodeTable[2*i+1] = node->value;
    } else {
      decodeTable[2*i] = node->left->nodeId;
      decodeTable[2*i+1] = node->right->nodeId;
    }
  }

  /*
  printf("Decode table\n");
  for (size_t i=0; i < treeOrderedNodes.size(); i++) {
    printf("%3d: %d  %d\n", (int)i, decodeTable[2*i], decodeTable[2*i+1]);
  }
  */

  decoder.init(decodeTable);
  if (isDuplicateValueSet())
    decoder.setDuplicate(size, dupValue);
}


void Huffman::orderNodes(std::vector<Node*> &orderedNodes) {

  queue<Node*> searchQueue;
  searchQueue.push(root);
  while (!searchQueue.empty()) {
    Node *node = searchQueue.front();
    searchQueue.pop();
    node->nodeId = orderedNodes.size();
    orderedNodes.push_back(node);
    if (!node->isLeaf()) {
      searchQueue.push(node->left);
      searchQueue.push(node->right);
    }
  }
}  


int Huffman::encodedLength(int value) {
  assert(computed);
  // special case for dup value
  if (isDuplicateValueSet() && value == size)
    return encodeTable[value*2] + HuffmanDup::DUP_ENCODE_BITS;
  assert(value >= 0 && value < size);

  return encodeTable[value*2];
}

int Huffman::totalEncodedLengthBytes() {
  if (!computed) return -1;
  int total = 0;
  for (int v = 0; v < size; v++)
    total += encodedLength(v) * getCount(v);

  // add the cost of encoding duplicate strings
  if (isDuplicateValueSet())
    total += dupValueCount * encodedLength(size);

  return (total + 31) / 32 * 4;
}
  

string Huffman::Node::bitString() {
  int wordIdx = 0;
  unsigned mask = 1;
  string str(bitLength, '0');
  for (int i=0; i < bitLength; i++) {
    unsigned bit = encoding[wordIdx] & mask;
    if (bit) str[i] = '1';
    mask <<= 1;
    if (mask == 0) {
      mask = 1;
      wordIdx++;
    }
  }

  return str;
}

