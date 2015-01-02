#ifndef __BIT_STREAM_H__
#define __BIT_STREAM_H__

#include <cstdio>
#include <cstring>
#include <cassert>


/**
  Little-endian bit stream.

  Each word of data is filled starting with the least significant bits.
  For example, if the string [1, 0, 1, 1] is written, there will be one word
  of data in the output; an integer with a value of 13.
*/
class BitStreamWriter {
  FILE *outf;
  unsigned *buffer;

  /** Size of the buffer */
  int capacity;

  /** Number of full words used in the buffer (the next bit to be added will
      go into buffer[used]. */
  int wordsBuffered;

  /** Number of bits full in the current word. */
  int bitsBuffered;

  /** Total number of bits passed to either of the write() functions. */
  size_t bitsWritten;

 public:
  /** capacity_ is the size of the buffer in words (not bytes) */
  BitStreamWriter(FILE *outf_, int capacity_=1024)
    : outf(outf_), capacity(capacity_), wordsBuffered(0), bitsBuffered(0),
      bitsWritten(0) {

    buffer = new unsigned[capacity];
    clearBuffer();

  }

  ~BitStreamWriter() {
    flush();
    delete[] buffer;
  }

  /** Add 32 or less bits to the stream. The bits are contains in 'bits',
      and 'count' is the number of bits being written. Least significant
      bits go first. */
  void write(unsigned bits, int count) {
    if (count <= 0) return;

    bitsWritten += count;

    // too much to fit in the current word
    if (bitsBuffered + count > 32) {
      // add the bits that will fit
      int bitsAdded = 32 - bitsBuffered;
      addBits(bits, bitsAdded);
      count -= bitsAdded;
      bits >>= bitsAdded;
    }
    
    addBits(bits, count);
  }      
  
  /** Add any number of bits to the stream. */
  void write(unsigned *bitArray, int count) {

    if (count <= 0) return;

    bitsWritten += count;

    // shortcut the quick case--just adding a few bits
    if (bitsBuffered + count <= 32) {
      addBits(bitArray[0], count);
      return;
    }
    
    // no offset; can copy whole words
    if (bitsBuffered == 0) {
      int wordCount = count / 32;  // # of whole words

      // if the incoming bits will fill the buffer, send it directly to outf
      if (wordCount + wordsBuffered >= capacity) {
	// flush the buffer
	flushWords();

	// write directory from the input data
	fwrite(bitArray, sizeof(unsigned), wordCount, outf);
      }

      // the incoming bits will not fill the buffer, so just add to it
      else {
	memcpy(&buffer[wordsBuffered], bitArray, sizeof(unsigned) * wordCount);
	wordsBuffered += wordCount;
      }
      count -= 32*wordCount;

      // add the last word
      write(bitArray[wordCount], count);
      return;
    }

    int inputBitOffset = 0;
    
    // fill the current word in the buffer
    if (bitsBuffered + count >= 32) {
      // add the bits that will fit
      int bitsAdded = 32 - bitsBuffered;
      addBits(bitArray[0], bitsAdded);
      count -= bitsAdded;
      inputBitOffset = bitsAdded;
    }

    // the current word should have been filled, leaving bitsBuffered==0.
    // if it wasn't, then all remaining bits must have been used up
    assert(bitsBuffered == 0 || count == 0);

    // keep adding words
    while (count >= 32) {
      buffer[wordsBuffered++] = 
	(bitArray[0] >> inputBitOffset) |
	(bitArray[1] << (32-inputBitOffset));

      count -= 32;
      bitArray++;

      if (wordsBuffered == capacity) flushWords();
    }

    // there are less than 32 bits remaining, and they may spread across
    // bitArray[0] and bitArray[1]
    bitsBuffered = 0;
    int bitsInWord0 = 32 - inputBitOffset;
    if (count <= bitsInWord0) {
      addBits(bitArray[0] >> inputBitOffset, count);
    } else {
      addBits(bitArray[0] >> inputBitOffset, bitsInWord0);
      count -= bitsInWord0;
      addBits(bitArray[1], count);
    }
  }      

  void flush() {
    // if the final word is incomplete, fill it
    if (bitsBuffered > 0)
      addBits(0, 32 - bitsBuffered);
    flushWords();
  }

  /** Returns the number of bits written to the stream. */
  size_t size() {
    return bitsWritten;
  }

  // return a mask for the given number of low bits
  static unsigned getMask(int bitCount) {
    if (bitCount <= 0) return 0;
    if (bitCount >= 32) return 0xffffffff;
    return (1 << bitCount) - 1;
  }

 private:

  void clearBuffer() {
    memset(buffer, 0, sizeof(unsigned) * capacity);
  }

  /** Assumes count + bitsBuffered <= 32. If not, bits will be lost. */
  void addBits(unsigned bits, int count) {
    bits &= getMask(count);
    buffer[wordsBuffered] |= bits << bitsBuffered;
    bitsBuffered += count;
    if (bitsBuffered >= 32) {
      wordsBuffered++;
      bitsBuffered = 0;
      if (wordsBuffered == capacity) flushWords();
    }
  }

  /** Write all the full words of data in the buffer to the output stream */
  void flushWords() {
    if (wordsBuffered == 0) return;

    fwrite(buffer, sizeof(unsigned), wordsBuffered, outf);

    // if the buffer wasn't full, the last word might have some
    // data in it, so copy it to the top
    unsigned leftoverBits = 0;
    if (wordsBuffered < capacity)
      leftoverBits = buffer[wordsBuffered];

    clearBuffer();
    buffer[0] = leftoverBits;
    wordsBuffered = 0;
  }
      
};


/**
   This is included as a simple unoptimized test case against which
   the fancier and buggier code can be tested.
   On my test machine, this runs at 8.4 MB/s, and BitStreamWriter
   runs at 128 MB/sec.
*/
class BitStreamWriterSimple {
  FILE *outf;
  unsigned buffer;
  int bitsBuffered;

 public:
  BitStreamWriterSimple(FILE *outf_, int unused=0)
    : outf(outf_), buffer(0), bitsBuffered(0) {
  }

  ~BitStreamWriterSimple() {
    flush();
  }
    

  /** Add 32 or less bits to the stream. The bits are contains in 'bits',
      and 'count' is the number of bits being written. Least significant
      bits go first. */
  void write(unsigned bits, int count) {
    for (int i=0; i < count; i++) {
      addBit(bits);
      bits >>= 1;
    }
  }      


  /** Add any number of bits to the stream. */
  void write(unsigned *bitArray, int count) {
    while (count > 32) {
      write(*bitArray++, 32);
      count -= 32;
    }
    write(*bitArray, count);
  }

  void flush() {
    if (bitsBuffered) {
      fwrite(&buffer, sizeof(unsigned), 1, outf);
      bitsBuffered = 0;
      buffer = 0;
    }
  }

  void addBit(unsigned bit) {
    bit &= 1;
    if (bitsBuffered == 0) {
      buffer = bit;
      bitsBuffered++;
    } else {
      buffer |= bit << bitsBuffered;
      bitsBuffered++;
      if (bitsBuffered == 32) flush();
    }
  }
};    
  
  

/**
   Ditto BitStreamWriterSimple. This is the simple version.
*/
class BitStreamReaderSimple {
  FILE *inf;
  unsigned buffer;
  int bitsBuffered;
  bool eof;

 public:
  BitStreamReaderSimple(FILE *inf_) : inf(inf_), buffer(0), bitsBuffered(0),
    eof(false) {

    refill();

  }

  // return the next 32 bits in the data
  // unsigned peek();
  // void peek(unsigned *bitArray, int count);

  // skip ahead this number of bits
  // void skip(int count);

  // return the number of bits consumed
  // size_t bitsConsumed();

  unsigned readBit() {
    if (bitsBuffered == 0) refill();
    unsigned result = buffer & 1;
    buffer >>= 1;
    bitsBuffered--;
    return result;
  } 

  unsigned read(int bitCount) {
    unsigned result = 0;
    if (bitCount > 32) bitCount = 32;
    for (int i=0; i < bitCount; i++) {
      result |= readBit() << i;
    }
    return result;
  }

  void read(unsigned *dest, int bitCount) {
    while (bitCount >= 32) {
      *dest++ = read(32);
      bitCount -= 32;
    }
    *dest++ = read(bitCount);
  }

  bool isEof() {
    return eof;
  }

 private:
  bool refill() {
    if (eof) return false;

    if (fread(&buffer, sizeof(unsigned), 1, inf) != 1) {
      eof = true;
      buffer = 0;
      bitsBuffered = 0;
      return false;
    } else {
      bitsBuffered = 32;
      return true;
    }
  }
};


class BitStreamReader {
  FILE *inf;
  unsigned *buffer;

  /** Size of the 'buffer' array */
  int capacity;
  
  /** Number of words in 'buffer' that contain data. For example,
      if the capacity is 1000 words but the input file contains only
      320 bytes, after filling the buffer this will be 80. */
  int wordsBuffered;

  /** Number of words and bits of the buffer that have been returned already
      by one of the 'read' methods. */
  int bitsUsed, wordsUsed;

  /** true iff EOF has been reached on 'inf' */
  bool eof;

 public:
  /** capacity_ is the size of the buffer in words (not bytes) */
  BitStreamReader(FILE *inf_, int capacity_=1024)
    : inf(inf_), capacity(capacity_), wordsBuffered(0),
      bitsUsed(0), wordsUsed(0), eof(false) {

    buffer = new unsigned[capacity];
    fillBuffer();
  }

  ~BitStreamReader() {
    delete[] buffer;
  }

  /** Read 32 or less bits. If there are not that many bits remaining
      in the data file, the missing bits will be 0, and eof will be set. */
  unsigned read(int bitsWanted) {
    if (bitsWanted <= 0 || isEmpty()) return 0;
    if (bitsWanted > 32) bitsWanted = 32;

    // the current word has enough bits
    if (bitsWanted <= 32-bitsUsed) {
      return getBits(bitsWanted);
    }

    // otherwise get some from this word, and some from the next
    else {
      int part = 32-bitsUsed;
      unsigned result = getBits(part);
      result |= getBits(bitsWanted - part) << part;
      return result;
    }
  }

  /** Read one bit. If empty, return 0. */
  unsigned readBit() {
    if (isEmpty()) return 0;

    unsigned bit = (buffer[wordsUsed] >> bitsUsed) & 1;

    bitsUsed++;
    if (bitsUsed == 32) {
      bitsUsed = 0;
      wordsUsed++;
      if (wordsUsed == wordsBuffered) fillBuffer();
    }

    return bit;
  }
  

  /** Read 'bitsWanted' bits into 'array' */
  void read(unsigned *array, int bitsWanted) {
    while (bitsWanted >= 32) {
      *array++ = read(32);
      bitsWanted -= 32;
    }
    *array = read(bitsWanted);
  }

  bool isEmpty() {
    return wordsUsed == wordsBuffered;
  }

  static unsigned getMask(int bits) {
    return BitStreamWriter::getMask(bits);
  }

 private:
  /** The caller must insure that count is less than the number of
      bits remaining in the current word: 32-bitsUsed. */
  unsigned getBits(int count) {
    unsigned result = (buffer[wordsUsed] >> bitsUsed) & getMask(count);
    bitsUsed += count;
    if (bitsUsed == 32) {
      bitsUsed = 0;
      wordsUsed++;
      if (wordsUsed == wordsBuffered) fillBuffer();
    }
    return result;
  }


  bool fillBuffer() {
    wordsBuffered = fread(buffer, sizeof(unsigned), capacity, inf);
    if (wordsBuffered < capacity) eof = true;
    bitsUsed = wordsUsed = 0;
    return wordsBuffered != 0;
  }
    
};


#endif // __BIT_STREAM_H__
