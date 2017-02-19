#!/usr/bin/env python

import sys
import arff

DEBUG=False

class CPT:
  def __init__(self, idx_class, idx_data, n_class, n_data):
    self.idx_class = idx_class
    self.idx_data = idx_data
    self.n_class = n_class
    self.n_data = n_data
    self.t_obs = 0
    self.tab = []
    for i in range(n_class):
      self.tab.append([0.0] * n_data)

  def __str__(self):
    attr_str = "{:>5}|" * self.n_class
    str_self = "{:<4}\\{:<4}|{}\n".format(
      "atr", "cls", attr_str.format(*range(self.n_class)))
    for d in range(self.n_data):
      str_self += "{:>9}|".format(d)
      for c in range(self.n_class):
        str_self += "{:>5}|".format(self.tab[c][d])
      str_self += "\n"
    return str_self

  def laplace_est(self):
    self.t_obs += self.n_class * self.n_data
    for i in range(self.n_class):
      for j in range(self.n_data):
        self.tab[i][j] += 1.0

  def observe(self, obs):
    # obs is a row of data for which to update the cpt
    self.tab[obs[self.idx_class]][obs[self.idx_data]] += 1.0
    self.t_obs += 1.0

  def proba(self, obs, cls = -1):
    # returns probability of the obs given class
    #   cls is optional, and otherwise taken from the observation
    if cls == -1:
      cls = obs[self.idx_class]
    t_obs_cls = sum(self.tab[cls])
    t_cond_obs = self.tab[cls][obs[self.idx_data]]
    return t_cond_obs/t_obs_cls

def create_cpts(arff):
  # returns a tuple: (cpts list, class counts[+,-])
  #   assumes the class attribute is the final attribute
  cpts = []
  n_class = len(arff.mapped_values[-1])
  n_attr = len(arff.header) - 1
  for attr in range(n_attr):
    n_attr_vals = len(arff.mapped_values[attr])
    cpts.append(CPT(n_attr, attr, n_class, n_attr_vals))
    cpts[-1].laplace_est()

  # feed data into cpts
  cls_counts = [1.0] * len(arff.mapped_values[-1])
  for row in arff.data:
    cls_counts[row[-1]] += 1.0
    for cpt in cpts:
      cpt.observe(row)

  return (cpts,cls_counts)

def predict(cpt_tpl, arff):
  # do predictions
  #   cpt_tpl: (cpts list, class counts[+,-])
  #            assumes 2 classes
  # returns a list of lists of predictions [[cls1, cls2],... ]

  cpts = cpt_tpl[0]
  cls_counts = cpt_tpl[1]
  ncls = len(cls_counts)
  tcls = sum(cls_counts)
  pcls = [float(i)/tcls for i in cls_counts]

  preds = []
  for row in arff.data:
    prod_pcond = [1.0] * ncls
    for cpt in cpts:
      for cls in range(ncls):
        prod_pcond[cls] *= cpt.proba(row, cls)
    
    numers = [pcls[i] * prod_pcond[i] for i in range(ncls)]
    tnumers = sum(numers)
    preds.append([i / tnumers for i in numers])
  return preds

def learn_arff(train_arff, test_arff):
  cpt_tpl = create_cpts(train_arff)
  preds = predict(cpt_tpl, test_arff)
  correct = 0
  for idx, pred in enumerate(preds):
    proba = max(pred)
    cls_pred = pred.index(proba)
    correct += 1 if cls_pred == test_arff.data[idx][-1] else 0
  return correct

def learn(trainfile, testfile):
  train_arff = arff.read_arff(trainfile)
  test_arff = arff.read_arff(testfile)

  cpt_tpl = create_cpts(train_arff)

  if DEBUG:
    for tbl in cpt_tpl[0]:
      print "{}".format(train_arff.header[tbl.idx_data])
      print tbl
      print ""

  preds = predict(cpt_tpl, test_arff)

  for attr in train_arff.header[:-1]:
    print "{} class".format(attr)
  print ""

  correct = 0
  for idx, pred in enumerate(preds):
    proba = max(pred)
    cls_pred = pred.index(proba)
    correct += 1 if cls_pred == test_arff.data[idx][-1] else 0
    print "{} {} {:<.16}".format(
      train_arff.mapped_values[-1][cls_pred],
      train_arff.mapped_values[-1][test_arff.data[idx][-1]],
      proba)

  print "\n{}".format(correct)

def main():
  if(len(sys.argv) == 1):
    print "useage: {} <train arff> <test arff>".format(sys.argv[0])
    exit(0)

  learn(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
  main()