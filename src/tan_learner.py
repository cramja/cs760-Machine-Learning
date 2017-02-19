#!/usr/bin/env python

import sys
import arff
import math
from decimal import *

DEBUG=False

class CPT:
  # a 3 attribute cpt
  def __init__(self, arff, attrs):
    # attrs is a list of 3 attribute ids in the arff
    self.arff = arff
    self.attr = attrs
    self.nattr = [len(arff.mapped_values[attr]) for attr in self.attr]
    self.tbl = []
    self.tobs = 0.0
    for i in range(self.nattr[-1]):
      self.tbl.append([0.0] * (self.nattr[0] * self.nattr[1]))

  def __eq__(self, other):
    return self.attr == other.attr

  def __hash__(self):
    return hash(self.attr[0] * 1e6 + self.attr[1] * 1e3 + self.attr[2])

  def __str__(self):
    name_attrs = [self.arff.mapped_values[a] for a in self.attr]
    tpairs = self.nattr[0] * self.nattr[1]
    attr_str = "{:>5}|" * self.nattr[0] * self.nattr[1]
    str_self = "{:<3}\\{:>5}|{}\n".format(
      "", 
      self.attr[0],
      attr_str.format(*[i/self.nattr[1] for i in range(tpairs)]))
    str_self += "{:<3}\\{:>5}|{}\n".format(
      self.attr[2], 
      self.attr[1],
      attr_str.format(*[i%self.nattr[1] for i in range(tpairs)]))
    for i in range(self.nattr[2]):
      str_self += "{:>9}|".format(i)
      for j in range(tpairs):
        str_self += "{:>5}|".format(self.tbl[i][j])
      str_self += "{:>5}\n".format(sum(self.tbl[i]))
    return str_self

  def observe(self, row):
    c = row[self.attr[2]]
    a = (row[self.attr[0]] * self.nattr[1]) + row[self.attr[1]]
    self.tbl[c][a] += 1.0
    self.tobs += 1.0

  def prob_a0a1_given_c(self, a0, a1, c):
    # prob(a0, a1 | c)
    return (self.tbl[c][a0 * self.nattr[1] + a1] + 1.0) / (sum(self.tbl[c]) + (self.nattr[0] * self.nattr[1]))

  def prob_a0_given_c(self, a0, c):
    # prob(a0 | c)
    return (sum(self.tbl[c][a0 * self.nattr[1] : (a0 + 1) * self.nattr[1]]) + 1.0) / (sum(self.tbl[c]) + self.nattr[0])

  def prob_a1_given_c(self, a1, c):
    # prob(a1 | c)
    return (sum([v for idx, v in enumerate(self.tbl[c]) if (idx - a1) % self.nattr[1] == 0]) + 1.0) / (sum(self.tbl[c]) + self.nattr[1])

  def prob_c_given_a0a1(self, c, a0, a1):
    cnt_given = Decimal(0.0)
    cnt_c = Decimal(0.0)
    for i in range(self.nattr[2]):
      cnt_given += Decimal(self.tbl[i][a0 * self.nattr[1] + a1])
      cnt_c += Decimal(0.0) if i != c else Decimal(self.tbl[i][a0 * self.nattr[1] + a1])
    return (cnt_c + Decimal(1.0)) / (cnt_given + self.nattr[2])


  def prob_a0a1c(self, a0, a1, c):
    # prob(a0, a1, c)
    return (self.tbl[c][a0 * self.nattr[1] + a1] + 1.0) / (self.tobs + (self.nattr[2] * self.nattr[1] * self.nattr[0]))

def create_cpts(arff):
  # create every combination of attribute pairs into a cpt
  cpts = {}
  tattr = len(arff.header) - 1
  for i in range(tattr):
    for j in xrange(i + 1, tattr):
      cpts[(i,j)] = CPT(arff, [i,j,len(arff.header) - 1])
  for row in arff.data:
    for k, cpt in cpts.iteritems():
      cpt.observe(row)
  return cpts

def CMutualInfo(cpts, a0, a1):
  # conditional mutual information
  cpt = cpts[(a0, a1)]
  s = 0.0
  for atr0 in range(cpt.nattr[0]):
    for atr1 in range(cpt.nattr[1]):
      for c in range(cpt.nattr[2]):
        numer = cpt.prob_a0a1_given_c(atr0, atr1, c)
        denom = cpt.prob_a0_given_c(atr0, c) * cpt.prob_a1_given_c(atr1, c)
        s += (cpt.prob_a0a1c(atr0, atr1, c) * math.log(numer/denom, 2))
  return s

class Edge:
  def __init__(self, n0, n1, weight):
    self.n0 = n0
    self.n1 = n1
    self.directed = False
    self.weight = weight

  def __str__(self):
    return "{:>2}{:>3}{:>2}, w: {}".format(
      self.n0, 
      "->" if self.directed else "<->",
      self.n1,
      self.weight)

  def __repr__(self):
    return "({},{})".format(self.n0, self.n1)

class Node:
  def __init__(self, attr):
    self.attr = attr
    self.edges = []

class Graph:

  def __init__(self, edges, nodes = -1):
    self.edges = edges
    self.nodes = nodes if nodes != -1 else []
    if len(self.nodes) == 0:
      for e in edges:
        if e.n0 not in self.nodes:
          self.nodes.append(e.n0)
        if e.n1 not in self.nodes:
          self.nodes.append(e.n1)

  def __str__(self):
    return "\n".join([str(e) for e in self.edges])

  def get_node(self, n):
    for node in self.nodes:
      if node.attr == n:
        return node

  def prims(self):
    frontier = [e for e in self.nodes[0].edges]
    maxEdges = []
    covered = [self.nodes[0].attr]
    while len(covered) != len(self.nodes):
      me = frontier[0]
      for e in frontier[1:]:
        if e.weight > me.weight:
          me = e
        elif e.weight == me.weight:
          # Selection criteria given by the assignment:
          srcE = me.n0 
          dstE = me.n1 
          if me.n1 in covered:
            tmp = srcE; srcE = dstE; dstE = tmp;
          srce = e.n0
          dste = e.n1
          if e.n1 in covered:
            tmp = srce; srce = dste; dste = tmp;
          if srcE > srce:
            me = e
          elif srcE == srce:
            if dstE > dste:
              me = e

      if me.n0 not in covered:
        covered.append(me.n0)
        tmp = me.n0; me.n0 = me.n1; me.n1 = tmp;
      else:
        covered.append(me.n1)
      me.directed = True
      maxEdges.append(me)

      frontier = []
      for n in covered:
        for e in self.get_node(n).edges:
          if e.n0 not in covered or e.n1 not in covered:
            frontier.append(e)

    return Graph(maxEdges)

def graph_from_cpts(cpts, arff):
  nodes = []
  edges = []
  for attr in range(len(arff.header) - 1):
    nodes.append(Node(attr))
  for i, n0 in enumerate(nodes):
    for j in xrange(i + 1, len(arff.header) - 1):
      n1 = nodes[j]
      edges.append(Edge(i,j,CMutualInfo(cpts, n0.attr, n1.attr)))
      n0.edges.append(edges[-1])
      n1.edges.append(edges[-1])
  return Graph(edges, nodes)

class TAN:
  def __init__(self, mst, arff):
    self.arff = arff
    self.mst = mst
    idx_cls = len(arff.header) - 1
    self.cnt_cls = [0.0] * len(arff.mapped_values[idx_cls])
    self.root_attr = self.mst.edges[0].n0
    nroot_attr = len(arff.mapped_values[self.root_attr])
    self.cnt_root = [[0.0]*nroot_attr, [0.0]*nroot_attr]
    for e in self.mst.edges:
      e.cpt = CPT(arff, [idx_cls, e.n0, e.n1])
    for row in arff.data:
      self.cnt_cls[row[-1]] += 1.0
      self.cnt_root[row[-1]][row[self.root_attr]] += 1.0
      for e in self.mst.edges:
        e.cpt.observe(row)

  def predict(self, ex):
    rel_proba = []
    for cls in range(2):
      proba = Decimal(self.cnt_cls[cls] + 1.0) / Decimal(sum(self.cnt_cls) + 2.0)
      for edge in self.mst.edges:
        proba *= edge.cpt.prob_c_given_a0a1(ex[edge.n1], cls, ex[edge.n0])
      proba *= Decimal(self.cnt_root[cls][ex[self.root_attr]] + 1.0) / Decimal(sum(self.cnt_root[cls]) + len(self.cnt_root[cls]))
      rel_proba.append(proba)
    proba = [i/sum(rel_proba) for i in rel_proba]
    maxprob = max(proba)
    idxmax = proba.index(maxprob)
    return (idxmax, maxprob)

def learn_arff(train_arff, test_arff):
  cpts = create_cpts(train_arff)
  mst = graph_from_cpts(cpts, train_arff).prims()
  tan = TAN(mst, train_arff)
  correct = 0
  for row in test_arff.data:
    cls, proba = tan.predict(row)
    correct += 1 if cls == row[-1] else 0
  return correct

def learn(trainfile, testfile):
  train_arff = arff.read_arff(trainfile)
  test_arff = arff.read_arff(testfile)
  cpts = create_cpts(train_arff)

  if DEBUG:
    s = ""
    for i in xrange(0,len(train_arff.header) - 1):
      for j in xrange(0, len(train_arff.header) - 1):
        a0 = i if i < j else j
        a1 = i if i > j else j
        s += "{}, ".format("-1" if a0 == a1 else CMutualInfo(cpts, a0, a1))
      s += "\n"
    print s

  mst = graph_from_cpts(cpts, train_arff).prims()

  if DEBUG:
    print mst
    s = "["
    for e in mst.edges:
      s += "{},".format((e.n0, e.n1))
    print "{}]".format(s)

  tan = TAN(mst, train_arff)

  if DEBUG:
    # print cpts like david page does
    print mst.edges
    for edge in tan.mst.edges:
      print ""
      cpt = edge.cpt
      for c in range(cpt.nattr[2]):
        for a0 in range(cpt.nattr[0]):
          for a1 in range(cpt.nattr[1]):
            print "Pr({}={} | {}={},{}={}) = {}".format(
              cpt.attr[2], c,
              cpt.attr[0], a0,
              cpt.attr[1], a1, 
              cpt.prob_c_given_a0a1(c, a0, a1))


  # print tree structure
  print "{} class".format(train_arff.header[mst.edges[0].n0])
  for a in xrange(1, len(train_arff.header) - 1):
    for edge in tan.mst.edges:
      if edge.n1 == a:
        print "{} {} class".format(train_arff.header[edge.n1], train_arff.header[edge.n0])
        break
  print ""

  correct = 0
  for row in test_arff.data:
    cls, proba = tan.predict(row)
    print "{} {} {:<.16f}".format(
      train_arff.mapped_values[-1][cls],
      train_arff.mapped_values[-1][row[-1]],
      proba)
    correct += 1 if cls == row[-1] else 0

  print "\n{}".format(correct)

def main():
  if len(sys.argv) == 1:
    print "useage: {} <train arff> <test arff>".format(sys.argv[0])
    exit(0)
  learn(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
  main()