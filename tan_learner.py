#!/usr/bin/env python

import sys
import arff

class CPT:
  # a 2 attribute cpt

  def __init__(self, arff, attr0, attr1):
    self.arff = arff
    self.attr0 = attr0
    self.attr1 = attr1
    self.nattr0 = len(arff.mapped_values[attr0])
    self.nattr1 = len(arff.mapped_values[attr1])
    self.tbl = []
    self.tobs = self.nattr0 * self.nattr1 * 2.0
    for i in range(2):
      self.tbl.append([1.0] * (self.nattr0 * self.nattr1))

  def __eq__(self, other):
    return self.attr0 == other.attr0 and self.attr1 == other.attr1

  def __hash__(self):
    return hash(self.attr0 * 1000 + self.attr1)

  def __str__(self):
    name_cls = self.arff.mapped_values[len(self.arff.mapped_values) - 1]
    name_attr0 = self.arff.mapped_values[self.attr0]
    name_attr1 = self.arff.mapped_values[self.attr1]
    tpairs = self.nattr0 * self.nattr1
    attr_str = "{:>5}|" * self.nattr0 * self.nattr1
    str_self = "{:<3}\\{:>5}|{}\n".format(
      "", 
      self.arff.header[self.attr0][:5], 
      attr_str.format(
          *[name_attr0[i/len(name_attr1)][:5] for i in range(tpairs)]))
    str_self += "{:<3}\\{:>5}|{}\n".format(
      "cls", 
      self.arff.header[self.attr1][:5], 
      attr_str.format(
          *[name_attr1[i%len(name_attr1)][:5] for i in range(tpairs)]))
    for i in range(2):
      str_self += "{:>9}|".format(i)
      for j in range(tpairs):
        str_self += "{:>5}|".format(self.tbl[i][j])
      str_self += "\n"
    return str_self

  def observe(self, row):
    c = row[-1]
    a = row[self.attr0] * self.nattr1 + row[self.attr1]
    self.tbl[c][a] += 1.0
    self.tobs += 1.0

  def prob_a01_given_c(self, c):
    # prob(a0, a1 | c)
    return self.tbl[c][a0 * self.nattr1 + a1] / sum(self.tbl[c])

  def prob_a0_given_c(self, a0, c):
    # prob(a0 | c)
    return sum(self.tbl[c][a0 * self.nattr1 : (a0 + 1) * self.nattr1]) / sum(self.tbl[c])

  def prob_a1_given_c(self, a1, c):
    # prob(a1 | c)
    return sum([i for i in self.tbl[c] if i % a1 == 0]) / sum(self.tbl[c])

  def prob_a01c(self, a0, a1, c):
    # prob(a0, a1, c)
    return self.tbl[c][a0 * self.nattr1 + a1] / self.tobs

  def prob(self, a0 = -1, a1 = -1, c = -1):
    def given(a):
      return a != -1

    def ignore(a):
      return a == -1

    a_idx = []
    if ignore(a0) and ignore(a1) and given(c):
      # prob c
      return sum(tbl[c]) / self.tobs
  
    elif ignore(a0) and given(a1) and given(c):
      # prob(a1|c)
      a_idx = [i for i in range(self.nattr0 * self.nattr1) if i % self.nattr1 == a1]

    elif given(a0) and ignore(a1) and given(c):
      # prob(a0|c)
      a_idx = [self.nattr1 * a0 + i for i in range(self.nattr1)]

    elif given(a0) and ignore(a1) and ignore(c):
      # prob(a0)
      s = 0
      for i in range(2):
        s += sum(self.tbl[i][a0 * self.nattr1 : (a0 + 1) * self.nattr1])
      return s / self.tobs

    elif ignore(a0) and given(a1) and ignore(c):
      # prob(a1)
      s = 0
      for i in range(2):
        s += sum([j for j in self.tbl[i] if j % a1 == 0])
      return s / self.tobs

    elif given(a0) and given(a1) and given(c):
      # prob(a0,a1,c)
      return self.tbl[c][a0 * self.nattr1 + a1]/self.tobs

    else:
      raise ValueError("improper arguments: {}, {}, {}".format(a0, a1, c))

    tc = 0.0
    tcc = 0.0
    for ic in range(2):
      for ia in a_idx:
        tc += self.tbl[ic][ia]
        if ic == c:
          tcc += self.tbl[ic][ia]

    return tcc / tc


def create_cpts(arff):
  # create every combination of attribute pairs into a cpt
  cpts = []
  tattr = len(arff.header) - 1
  for i in range(tattr):
    for j in xrange(i + 1, tattr):
      cpts.append(CPT(arff, i, j))
  for row in arff.data:
    for cpt in cpts:
      cpt.observe(row)
  return cpts


def main():
  if(len(sys.argv) == 1):
    print "useage: {} <train arff> <test arff>".format(sys.argv[0])
    exit(0)
  train_arff = arff.read_arff(sys.argv[1])
  test_arff = arff.read_arff(sys.argv[2])

  cpts = create_cpts(train_arff)
  cpt = cpts[0]
  print cpt
  print cpt.tobs

  print "1,1,1: {}".format(cpt.prob(1,1,1)) 
  print "1,1,-1: {}".format(cpt.prob(1,1,-1)) 


if __name__ == '__main__':
  main()