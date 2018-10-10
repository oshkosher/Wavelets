#!/usr/bin/python

import sys

zero = 127
track_zero = True
run_length = 0

for line in sys.stdin:
  value = int(line)
  
  if track_zero:
    
    if value == zero:
      run_length += 1
    else:
      if run_length > 0:
        print run_length
      run_length = 0
    
  else:
    
    if value == zero:
      if run_length > 0:
        print run_length
      run_length = 0
    else:
      run_length += 1

if run_length > 0:
  print run_length
