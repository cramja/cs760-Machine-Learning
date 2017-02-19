#!/usr/bin/env python

import sys
from naive_learner import learn as naive_learn
from tan_learner import learn as tan_learn

def main():
  if len(sys.argv) < 4:
    print "usage: <train> <test> <n|t>"
  trainfile = sys.argv[1]
  testfile = sys.argv[2]
  algo = sys.argv[3]
  if algo == 'n':
    naive_learn(trainfile, testfile)
  elif algo == 't':
    tan_learn(trainfile, testfile)
  else:
    print "unknown algorithm: {}".format(algo)

if __name__ == '__main__':
  main()