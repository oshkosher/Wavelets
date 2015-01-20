#ifndef __BIT_STREAM_H__
#define __BIT_STREAM_H__

#include <cstdio>
#include <cstring>
#include <cassert>
#include <vector>


class BitStreamFileSink {
  FILE *outf;
  std::vector<unsigned> buffer;
  static const int BUFFER_SIZE = 1024*16;

 public:
  BitStreamFileSink(FILE *o) : outf(o) {
    buffer.reserve(BUFFER_SIZE);
  }
  
  void add(unsigned word) {
    buffer.push_back(word);

    if (buffer.size() == BUFFER_SIZE) flush();
  }

  void flush() {
    fwrite(buffer.data(), sizeof(unsigned), buffer.size(), outf);
    buffer.clear();
  }
};


class BitStreamMemorySink {
  std::vector<unsigned> *buffer;

 public:
  BitStreamMemorySink(std::vector<unsigned> *b) : buffer(b) {}

  std::vector<unsigned> *getBuffer() {return buffer;}
  
  void add(unsigned word) {
    buffer->push_back(word);
  }

  void flush() {}
};
    

/**
  Little-endian bit stream.

  Each word of data is filled starting with the least significant bits.
  For example, if the string [1, 0, 1, 1] is written, there will be one word
  of data in the output; an integer with a value of 13.

  WordSink is a class that defines:
    void add(unsigned word);
    void flush()
*/
template<class WordSink>
class BitStreamWriter {
  WordSink *wordSink;

  // 1-word buffer
  unsigned buffer;

  /** Number of bits full in the current word. */
  int bitsBuffered;

  /** Total number of bits passed to either of the write() functions. */
  size_t bitsWritten;

 public:
  /** capacity_ is the size of the buffer in words (not bytes) */
  BitStreamWriter(WordSink *ws)
    : wordSink(ws), bitsBuffered(0), bitsWritten(0) {
    buffer = 0;
  }

  ~BitStreamWriter() {
    wordSink->flush();
  }

  /** Add 32 or less bits to the stream. The bits are contains in 'bits',
      and 'count' is the number of bits being written. Least significant
      bits go first. */
  void write(unsigned bits, int count) {
    if (count <= 0) return;

    bitsWritten += count;

    writeInternal(bits, count);
  }      

 private:
  void writeInternal(unsigned bits, int count) {
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

 public:
  /** Add any number of bits to the stream.
      XXX - this could be optimized a bit. */
  void write(unsigned *bitArray, int count) {
    if (count <= 0) return;

    bitsWritten += count;

    int i=0;
    while (count >= 32) {
      writeInternal(bitArray[i++], 32);
      count -= 32;
    }
    writeInternal(bitArray[i], count);
  }

  void flush() {
    // if the final word is incomplete, fill it
    if (bitsBuffered > 0) {
      wordSink->add(buffer);
      bitsBuffered = 0;
      buffer = 0;
    }
    wordSink->flush();
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

  /** Assumes count + bitsBuffered <= 32. If not, bits will be lost. */
  void addBits(unsigned bits, int count) {
    bits &= getMask(count);
    buffer |= bits << bitsBuffered;
    bitsBuffered += count;
    if (bitsBuffered >= 32) {
      assert(bitsBuffered == 32);
      wordSink->add(buffer);
      buffer = 0;
      bitsBuffered = 0;
    }
  }
      
};


class BitStreamFileSource {
  FILE *inf;
  size_t pos;
  std::vector<unsigned> buffer;

 public:
  BitStreamFileSource(FILE *i) : buffer(1024*16,0) {
    inf = i;
    pos = 0;
    fillBuffer();
  }
  
  bool get(unsigned &value) {
    if (pos == buffer.size()) {
      fillBuffer();
    }
    if (pos >= buffer.size()) return false;
    value = buffer[pos++];
    return true;
  }

 private:
  void fillBuffer() {
    size_t numRead = fread(buffer.data(), sizeof(unsigned), buffer.capacity(),
                           inf);
    if (numRead < buffer.capacity()) buffer.resize(numRead);
    pos = 0;
  }
};


class BitStreamMemorySource {
  const std::vector<unsigned> *buffer;
  size_t pos = 0;

 public:
  BitStreamMemorySource(const std::vector<unsigned> *b) : buffer(b) {}

  bool get(unsigned &word) {
    if (pos < buffer->size()) {
      word = (*buffer)[pos++];
      return true;
    } else {
      word = 0;
      return false;
    }
  }

};
  

template<class WordSource>
class BitStreamReader {
  WordSource *wordSource;

  // 1-word buffer
  unsigned buffer;

  /** Number of bits of the buffer that have been returned already
      by one of the 'read' methods. */
  int bitsUsed;

  /** true iff EOF has been reached on 'inf' */
  bool eof;

 public:
  /** capacity_ is the size of the buffer in words (not bytes) */
 BitStreamReader(WordSource *ws) :
  wordSource(ws), buffer(0), bitsUsed(0), eof(false) {
    fillBuffer();
  }

  /** Read one bit. If empty, return 0. */
  unsigned readBit() {

    unsigned bit = (buffer >> bitsUsed) & 1;

    bitsUsed++;
    if (bitsUsed == 32) fillBuffer();

    return bit;
  }

  /** Read 32 or less bits. If there are not that many bits remaining
      in the data file, the missing bits will be 0, and eof will be set. */
  unsigned read(int bitsWanted) {
    assert(bitsWanted >= 0 && bitsWanted <= 32);
    if (eof || bitsWanted==0) return 0;

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
  

  /** Read 'bitsWanted' bits into 'array' */
  bool read(unsigned *array, int bitsWanted) {
    if (eof) return false;

    while (bitsWanted >= 32) {
      *array++ = read(32);
      bitsWanted -= 32;
    }
    *array = read(bitsWanted);
    return true;
  }

  static unsigned getMask(int bits) {
    return BitStreamWriter<void>::getMask(bits);
  }

  bool isEmpty() {return eof;}

 private:
  /** The caller must insure that count is less than the number of
      bits remaining in the current word: 32-bitsUsed. */
  unsigned getBits(int count) {
    assert(count <= 32-bitsUsed);
    unsigned result = (buffer >> bitsUsed) & getMask(count);
    bitsUsed += count;
    if (bitsUsed == 32) fillBuffer();
    return result;
  }


  bool fillBuffer() {
    if (wordSource->get(buffer)) {
      bitsUsed = 0;
      return true;
    } else {
      buffer = 0;
      bitsUsed = 0;
      eof = true;
      return false;
    }
  }
    
};


#endif // __BIT_STREAM_H__
