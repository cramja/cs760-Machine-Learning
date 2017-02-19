#!/usr/bin/env python

import sys
import arff
from naive_learner import learn_arff as naive_learn
from tan_learner import learn_arff as tan_learn

def main():
  trainfile = sys.argv[1]
  testfile = sys.argv[2]
  algo = sys.argv[3]

  ntrials = 4
  splitsize = [25, 50, 100]

  train_arff = arff.read_arff(trainfile)
  test_arff = arff.read_arff(testfile)
  fn = naive_learn if algo == 'n' else tan_learn
  print fn

  accuracies = []
  for size in splitsize:
    t_accuracy = []
    for trial in range(ntrials):
      t_accuracy.append(fn(train_arff.sample(size), test_arff))
    accuracies.append((float(sum(t_accuracy))/ntrials)/len(test_arff.data))
  for acc in accuracies:
    print "{}".format(acc)

if __name__ == '__main__':
  main()