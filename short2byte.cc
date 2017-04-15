#include <cstdio>

int main() {

  unsigned char buf[2];
  bool isLittleEndian = true;

  while (true) {

    if (2 == fread(buf, 1, 2, stdin)) {
      unsigned value = isLittleEndian
        ? (buf[0] + (buf[1] << 8))
        : (buf[1] + (buf[0] << 8));
      value >>= 8;
      putchar(value);
    } else {
      break;
    }

  }

  return 0;
}
