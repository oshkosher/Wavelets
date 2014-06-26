#!/usr/bin/python

import sys, random

def main(args):
  if len(args) != 2:
    print """
  random_matrix.py <width> <height>
"""
    return 1

  (width, height) = (int(args[0]), int(args[1]))

  print '%d cols %d rows' % (width, height)
  for y in xrange(height):
    for x in xrange(1, width):
      print '%.5f' % random.random(),
    print '%.5f' % random.random()



if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
