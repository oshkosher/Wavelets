#!/usr/bin/python

import sys

input_file = 'zero_runs_chicago'
REPLACEMENT_LEN = 12
COUNT_LEN = 10

"""
Take the output containing all the quantized values
head -n -4 foo | tail -n +10 | ./rle.py | sort -n | uniq -c > zero_runs_chicago


replacement len 6:
  chicago: 17..28  (142000 bits, 18kb saved, 4000 copies)
  collage: 17..28  (8500000 bits,  saved, 40000 copies)

replacement len 8:
  chicago: 27..35  (127000 bits, 16kb saved, 3000 copies)
  collage: 25..32  (7600000 bits saved, 33000 copies)

replacement len 9:
  chicago: 27..42 (121000 bits, 15kb saved, 2800 copies)  *this
  collage: 

replacement len 10:
  chicago: 30..42  (117000 bits, 14kb saved, 2600 copies)
  collage: 31..38  (6970000 bits saved, 27000 copies)

tag + 16 bit length
  8 bit quantization
  chicago: save 18kb  (replacement = 12 bits + count)
  collage: save 995kb  (replacement = 16 bits + count)

  7 bit quant
"""

run_lengths = []

for line in open(input_file):
  (count, length) = line.split()
  run_lengths.append((int(count), int(length)))

# replace a fixed number of copies with a reserved sequence
"""
for test_len in range(REPLACEMENT_LEN+1, 1000):
  saved = 0
  replacement_count = 0
  for (count, length) in run_lengths:
    if length > test_len:
      test_count = length // test_len
      savings_each = test_count * (test_len - REPLACEMENT_LEN)
      replacement_count += test_count
      saved += count * savings_each

  print('%d: %d %d' % (test_len, saved, replacement_count))
"""

# replace any number of copies with a reserved sequence and a 16-bit length
saved = 0
replacement_count = 0
for (count, length) in run_lengths:
  if length > REPLACEMENT_LEN + COUNT_LEN:
    savings_each = length - (REPLACEMENT_LEN + COUNT_LEN)
    saved += count * savings_each
    replacement_count += 1

print('%d bits (%d bytes) saved, %d instances' %
      (saved, saved/8, replacement_count))
