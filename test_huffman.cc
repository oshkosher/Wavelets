#include "huffman.h"


int main(int argc, char **argv) {
  
  Huffman huff(8);
  huff.update(0, 5);
  huff.update(1, 2);
  huff.update(2, 10);
  huff.update(3, 12);
  huff.update(4, 2);
  huff.update(5, 3);
  huff.update(6, 1);
  huff.update(7, 9);

  huff.computeHuffmanCoding();

  return 0;
}
