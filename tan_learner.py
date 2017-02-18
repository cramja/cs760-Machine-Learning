#!/usr/bin/env python

import sys
import arff
import math

DEBUG=False

class CPT:
  # a 2 attribute cpt

  def __init__(self, arff, attr0, attr1):
    self.arff = arff
    self.attr0 = attr0
    self.attr1 = attr1
    self.nattr0 = len(arff.mapped_values[attr0])
    self.nattr1 = len(arff.mapped_values[attr1])
    self.tbl = []
    self.tobs = 0.0 #self.nattr0 * self.nattr1 * 2.0
    for i in range(2):
      self.tbl.append([0.0] * (self.nattr0 * self.nattr1))

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
      # attr_str.format(
      #     *[name_attr0[i/len(name_attr1)][:5] for i in range(tpairs)]))
      attr_str.format(
          *[i/self.nattr1 for i in range(tpairs)]))
    str_self += "{:<3}\\{:>5}|{}\n".format(
      "cls", 
      self.arff.header[self.attr1][:5], 
      # attr_str.format(
      #     *[name_attr1[i%len(name_attr1)][:5] for i in range(tpairs)]))
      attr_str.format(
          *[i%self.nattr1 for i in range(tpairs)]))
    for i in range(2):
      str_self += "{:>9}|".format(i)
      for j in range(tpairs):
        str_self += "{:>5}|".format(self.tbl[i][j])
      str_self += "{:>5}\n".format(sum(self.tbl[i]))
    return str_self

  def observe(self, row):
    c = row[-1]
    a = (row[self.attr0] * self.nattr1) + row[self.attr1]
    self.tbl[c][a] += 1.0
    self.tobs += 1.0

  def prob_a0a1_given_c(self, a0, a1, c):
    # prob(a0, a1 | c)
    return (self.tbl[c][a0 * self.nattr1 + a1] + 1.0) / (sum(self.tbl[c]) + (self.nattr0 * self.nattr1))

  def prob_a0_given_c(self, a0, c):
    # prob(a0 | c)
    return (sum(self.tbl[c][a0 * self.nattr1 : (a0 + 1) * self.nattr1]) + 1.0) / (sum(self.tbl[c]) + self.nattr0)

  def prob_a1_given_c(self, a1, c):
    # prob(a1 | c)
    return (sum([v for idx, v in enumerate(self.tbl[c]) if (idx - a1) % self.nattr1 == 0]) + 1.0) / (sum(self.tbl[c]) + self.nattr1)

  def prob_a0a1c(self, a0, a1, c):
    # prob(a0, a1, c)
    return (self.tbl[c][a0 * self.nattr1 + a1] + 1.0) / (self.tobs + (2 * self.nattr1 * self.nattr0))

def create_cpts(arff):
  # create every combination of attribute pairs into a cpt
  cpts = {}
  tattr = len(arff.header) - 1
  for i in range(tattr):
    for j in xrange(i + 1, tattr):
      cpts[(i,j)] = CPT(arff, i, j)
  for row in arff.data:
    for k, cpt in cpts.iteritems():
      cpt.observe(row)
  return cpts

def CMutualInfo(cpts, a0, a1):
  # conditional mutual information
  cpt = cpts[(a0, a1)]
  s = 0.0
  for atr0 in range(cpt.nattr0):
    for atr1 in range(cpt.nattr1):
      for c in range(2):
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
            swap(srcE, dstE)
          srce = e.n0
          dste = e.n1
          if e.n1 in covered:
            swap(srce, dste)
          if srcE > srce:
            me = e
          elif srcE == srce:
            if dstE > dste:
              me = e

      if me.n0 not in covered:
        covered.append(me.n0)
      else:
        covered.append(me.n1)

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

def ghetto_proba_a0a1_given_c(arff, a0, a0v, a1, a1v, cv):
  total_c = 0.0
  total_a0a1 = 0.0
  for line in arff.data:
    if line[-1] == cv:
      total_c += 1.0
      total_a0a1 += 1.0 if line[a0] == a0v and line[a1] == a1v else 0.0
  # laplace
  total_c += len(arff.mapped_values[a0]) * len(arff.mapped_values[a1])
  total_a0a1 += 1
  return total_a0a1 / total_c

def ghetto_proba_a0a1c(arff, a0, a0v, a1, a1v, cv):
  total = 0.0
  total_a0a1c = 0.0
  for line in arff.data:
    total_a0a1c += (1.0 if line[a0] == a0v and line[a1] == a1v and line[-1] == cv else 0.0)
    total += 1
  # laplace
  total += len(arff.mapped_values[a0]) * len(arff.mapped_values[a1]) * 2
  total_a0a1c += 1
  return total_a0a1c / total

def ghetto_proba_a0_given_c(arff, a0, a0v, cv):
  total_c = 0.0
  total_a0 = 0.0
  for line in arff.data:
    if line[-1] == cv:
      total_c += 1.0
      total_a0 += 1.0 if line[a0] == a0v else 0.0
  # laplace
  total_c += len(arff.mapped_values[a0])
  total_a0 += 1
  return total_a0 / total_c

def ghetto_cmi(arff, a0, a1):
  na0 = len(arff.mapped_values[a0])
  na1 = len(arff.mapped_values[a1])
  s = 0.0
  for i in range(na0):
    for j in range(na1):
      for c in range(2):
        n = ghetto_proba_a0a1_given_c(arff,a0,i,a1,j,c)
        d = ghetto_proba_a0_given_c(arff, a0, i, c) * ghetto_proba_a0_given_c(arff, a1, j, c)
        s += (ghetto_proba_a0a1c(arff, a0, i, a1, j, c) * math.log(n/d,2.0))
  return s

def main():
  if(len(sys.argv) == 1):
    print "useage: {} <train arff> <test arff>".format(sys.argv[0])
    exit(0)
  train_arff = arff.read_arff(sys.argv[1])
  test_arff = arff.read_arff(sys.argv[2])

  print train_arff.mapped_values

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
    print gms
    s = "["
    for e in gm.edges:
      s += "{},".format((e.n0, e.n1))
    print "{}]".format(s)

  



if __name__ == '__main__':
  main()