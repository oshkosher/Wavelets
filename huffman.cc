#include <queue>
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
  priority_queue<Node*, vector<Node*>, CompareHuffmanNode> queue;

  Node *nodes = new Node[size * 2 - 1];
  longestBitString = longestUsedBitString = 0;

  // the first 'size' nodes are leaf nodes, the rest are interior nodes
  for (int i=0; i < size; i++) {
    nodes[i].initLeaf(i, counts[i]);
    queue.push(&nodes[i]);
  }

  // for N elements, do N-1 iterations
  for (int i=size; i < size*2-1; i++) {
    Node *a = queue.top(); queue.pop();
    Node *b = queue.top(); queue.pop();
    nodes[i].initInternalNode(a, b);
    queue.push(&nodes[i]);
  }

  Node *root = queue.top(); queue.pop();
  printTree(root);

  buildTables(root);

  delete[] nodes;

  computed = true;
}


void Huffman::buildTables(Node *root) {

  // order the nodes and leaf nodes
  vector<Node *> nodeTable, valueTable;
  
  // Do a breadth-first traversal of the tree to order the nodes by
  // increasing encoding length, and set their nodeId to their position
  // in the ordered list.
  orderNodes(root, nodeTable, valueTable);

  // compute the encodeTable
  computeEncodeTable(root, valueTable);

  printf("Decode table\n");

  // now build a decode table by converting the nodes and references
  // to their children into an array of pairs of indices.
  for (int i=0; i < size*2 - 1; i++) {
    Node *node = nodeTable[i];
    if (node->isLeaf()) {
      decodeTable[2*i] = -1;
      decodeTable[2*i+1] = node->value;
    } else {
      decodeTable[2*i] = node->left->nodeId;
      decodeTable[2*i+1] = node->right->nodeId;
    }

    printf("%3d: %d  %d\n", i, decodeTable[2*i], decodeTable[2*i+1]);
  }

  // print the encoding for each string
  for (size_t i=0; i < valueTable.size(); i++) {
    Node *node = valueTable[i];
    string bits = node->bitString();
    printf("%3d: %s\n", node->value, bits.c_str());
  }
}



void Huffman::orderNodes(Node *root, std::vector<Node*> &nodeTable,
			 std::vector<Node*> &valueTable) {
  nodeTable.clear();
  queue<Node*> searchQueue;
  searchQueue.push(root);
  while (!searchQueue.empty()) {
    Node *node = searchQueue.front();
    searchQueue.pop();
    node->nodeId = nodeTable.size();
    nodeTable.push_back(node);
    if (node->isLeaf()) {
      valueTable.push_back(node);
    } else {
      searchQueue.push(node->left);
      searchQueue.push(node->right);
    }
  }
}  

void Huffman::computeEncodeTable(Node *root,
				 const std::vector<Node*> &valueTable) {

  // traverse the tree, saving the encoded form of each leaf node
  BitStack bits;
  computeEncodedBits(root, bits);

  // build encodeTable
  encodeTable.resize(size*2);
  for (size_t i=0; i < valueTable.size(); i++) {
    Node *node = valueTable[i];
    encodeTable[node->value*2] = node->bitLength;
    encodeTable[node->value*2+1] = encodeTable.size();

    for (size_t wordNo=0; wordNo < node->encoding.size(); wordNo++)
      encodeTable.push_back(node->encoding[wordNo]);
  }

  printf("Encode table\n");
  for (unsigned i=0; i < encodeTable.size(); i++) {
    printf("%3u: %d\n", i, encodeTable[i]);
  }
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

    if (node->count > 0 && node->bitLength > longestUsedBitString)
      longestUsedBitString = node->bitLength;

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


int Huffman::encodedLength(int value) {
  if (!computed) return -1;
  assert(value >= 0 && value < size);

  return encodeTable[value*2];
}
  
/*
bool Huffman::encodeWord(int value, int &length, unsigned &bits) {
  assert(computed);

  int dataOffset = encodeTable[value];

  length = encodeData[dataOffset];
  if (length < 0 || length > 32) return false;

  bits = encodeData[dataOffset+1];
}
*/

/*
bool Huffman::decode(unsigned bits, int &value, int &length) {
  assert(computed);
  pos = 0;
  
  for (int i=1; ; i++) {
    // 1 means go right, 0 means go left
    if (bits & 1) {
      pos = decodeTable[pos+1];
    } else {
      pos = decodeTable[pos];
    }

    // if the left index is -1, we found the value
    if (decodeTable[pos] == -1) {
      value = decodeTable[pos+1];
      length = i;
      return;
    }
  }
}
*/

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
