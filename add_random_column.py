#!/usr/bin/python

import random, sys

random.seed()

for line in sys.stdin:
  if line[-1] == '\n': line = line[:-1]
  line = line + '\t' + ('%.6f' % random.random())
  print line



