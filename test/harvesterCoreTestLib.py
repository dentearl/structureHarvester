"""
Simple data structure for holding test data.
"""

class Example(object):
  def __init__(self, name, files):
    self.name = name  # results directory name
    self.files = files  # tuple, [0] filename [1] entirety of file
