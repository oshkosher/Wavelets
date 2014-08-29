#include <queue>
#include <cstdio>
#include "huffman.h"

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

  nodes = new Node*[size];
  for (int i=0; i < size; i++) {
    nodes[i] = new Node(i, counts[i]);
    queue.push(nodes[i]);
  }

  /*
  // list the queue
  while (!queue.empty()) {
    Node *node = queue.top();
    queue.pop();
    printf("%d: %d\n", node->value, node->count);
  }
  */

  // for N elements, do N-1 iterations
  for (int i=1; i < size; i++) {
    Node *a = queue.top(); queue.pop();
    Node *b = queue.top(); queue.pop();
    Node *combine = new Node(a, b);
    queue.push(combine);
  }

  // printf("queue size: %d\n", (int)queue.size());
  Node *tree = queue.top();
  // printTree(tree);

  computed = true;
}

int Huffman::encodedLength(int value) {
  if (!computed) return -1;
  assert(value >= 0 && value < size);

  Node *n = nodes[value];
  int depth = 0;
  while (n->parent) {
    depth++;
    n = n->parent;
  }
  
  return depth;
}
  
  

void Huffman::encode(int value, int &length, int &bits) {
}

void Huffman::decode(int bits, int &value) {
}


