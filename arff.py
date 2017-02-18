#!/usr/bin/env python

import sys

def main():
  if(len(sys.argv) == 1):
    print "useage: {} [input file]".format(sys.argv[0])
    exit(0)
  print read_arff(sys.argv[1])

class Arff:
  def __init__(self):
    self.relation = ""
    self.header = []
    self.mapped_values = []
    self.data = []

  def __repr__(self):
    self_str = "{}\n\n".format(self.relation)
    for attr in self.header:
      self_str += (attr + "\n")
    return self_str

  def __str__(self):
    return self.__repr__()

  def name(self, attr = -1):
    if attr == -1:
      return self.relation
    return self.header[attr]


def read_data_line(lnum, line, arffobj):
  line_seg = line.split(",")
  if len(line_seg) != len(arffobj.header):
    raise ValueError("[{}] wrong num attrs ({} vs {}) in data: \"{}\"".format(
      lnum, len(line_seg), len(arffobj.header), line))
  row = []
  for i, attr in enumerate(line_seg):
    for j, val in enumerate(arffobj.mapped_values[i]):
      if val == attr:
        row.append(j)
        break
    if len(row) == i:
      raise ValueError("[{}] unknown attr found ({}) in data: \"{}\"".format(lnum, attr, line))
  arffobj.data.append(row)
  

def read_attribute_line(lnum, line, arffobj):
  # 'lymphatics' { normal, arched, deformed, displaced}
  line_seg = line.split()
  if len(line_seg) < 2:
    raise ValueError("[{}] unknown attr line fmt: \"{}\"".format(lnum, line))
  attr_n = line_seg[1][1:-1]
  mapped = []
  for attr in line_seg[2:]:
    if attr.startswith(",") or attr.startswith("{"):
      attr = attr[1:]
    elif attr.endswith(",") or attr.endswith("}"):
      attr = attr[:-1]
    if len(attr) == 0:
      continue
    mapped.append(attr)
  arffobj.header.append(attr_n)
  arffobj.mapped_values.append(mapped)


def read_arff(f):
  text = open(f, 'r').read()
  arffobj = Arff()
  data_mode = False
  lnum = -1
  for line in text.split("\n"):
    lnum += 1
    line = line.strip()
    if len(line) == 0 or line.startswith("%"):
      continue

    if data_mode:
      read_data_line(lnum, line, arffobj)
    else:
      if line.startswith('@relation'):
        arffobj.relation = line.split(" ")[-1]
      elif line.startswith('@attribute'):
        read_attribute_line(lnum, line, arffobj)
      elif line == "@data":
        data_mode = True
      else:
        raise ValueError("[{}] unknown line header: \"{}\"".format(lnum, line))
  return arffobj



if __name__ == '__main__':
  main()