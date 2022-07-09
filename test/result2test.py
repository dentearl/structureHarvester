#!/usr/bin/env python3
""" result2test
dent earl, dearl a soe ucsc edu
2014

script to take a structure result directory and kick to stdout
output that can be pasted into test_data.py and used in unittests
"""
from argparse import ArgumentParser, ArgumentTypeError
from glob import glob
import os
import sys


def DirType(d):
  """ given a string path to a directory, D, verify it can be used.
  """
  d = os.path.abspath(d)
  if not os.path.exists(d):
    raise ArgumentTypeError('DirType:%s does not exist' % d)
  if not os.path.isdir(d):
    raise ArgumentTypeError('DirType:%s is not a directory' % d)
  if os.access(d, os.R_OK):
    return d
  else:
    raise ArgumentTypeError('DirType:%s is not a readable dir' % d)


def initializeArguments(parser):
  parser.add_argument('--resultDir', type=DirType)
  parser.add_argument('--name', type=str)


def checkArguments(args, parser):
  pairs = tuple((item, getattr(args, item) )for item in ['resultDir'])
  for name, value in pairs:
    if value is None:
      parser.error('Specify --%s' % name)


def processDir(directory, name):
  if name is None:
    name = os.path.basename(directory)
  print('import test_harvesterCore as t_hc\n')
  print("e = t_hc.Example('%s'," % name)
  print("                 [")
  for path in glob(os.path.join(directory, '*_f')):
    with open(path, 'r') as f:
      print("                   ('%s'," % os.path.basename(path))
      print("                    '''")
      for line in f:
        line = line.strip()
        print(line)
      print("'''),")
  print("])")


def main():
  parser = ArgumentParser()
  initializeArguments(parser)
  args = parser.parse_args()
  checkArguments(args, parser)
  processDir(args.resultDir, args.name)


if __name__ == '__main__':
  main()
