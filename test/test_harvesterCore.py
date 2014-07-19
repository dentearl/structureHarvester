""" unit tests for the Structure Harvester.
dent earl, dearl a soe ucsc edu
July 2014


"""
from glob import glob
import imp
import os
import shutil
import string
import subprocess
import sys
import unittest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]))))
import harvesterCore as hc
import test_data as td


class Example(object):
  def __init__(self, name, files):
    self.name = name  # results directory name
    self.files = files  # tuple, [0] filename [1] entirety of file


def makeTempDirParent():
  """ make the parent temp dir directory
  """
  if not os.path.exists(os.path.join(os.curdir, '.tempTestDir')):
    os.mkdir(os.path.join(os.curdir, '.tempTestDir'))
  return os.path.join(os.curdir, '.tempTestDir')


def removeTempDirParent():
  """ remove the parent temp dir directory
  """
  if os.path.exists(os.path.join(os.curdir, '.tempTestDir')):
    shutil.rmtree(os.path.join(os.curdir, '.tempTestDir'))

def removeDirIfEmpty(d):
  """ remove the parent temp dir directory
  """
  if glob(os.path.join(d, '*')) == []:
    shutil.rmtree(d)


def makeTempDir(name=None):
  """ make the typical directory where all temporary test files will be stored.
  """
  makeTempDirParent()
  charSet = string.ascii_lowercase + '123456789'
  if name is None:
    while True:
      name = '%s_%s' % (''.join(random.choice(charSet) for x in xrange(4)),
                        ''.join(random.choice(charSet) for x in xrange(4)))
      if not os.path.exists(os.path.join(os.curdir, '.tempTestDir', name)):
        break
  if not os.path.exists(os.path.join(os.curdir, '.tempTestDir', name)):
    os.mkdir(os.path.join(os.curdir, '.tempTestDir', name))
  return os.path.join(os.curdir, '.tempTestDir', name)


def removeDir(dirpath):
  """ destroy a directory
  """
  if os.path.exists(dirpath):
    shutil.rmtree(dirpath)
  if glob(os.path.join(os.path.dirname(dirpath), '*')) == []:
    # if this is the last tempDir to be destroyed, destroy the parent.
    removeTempDirParent()


def which(program):
  """which() acts like the unix utility which, but is portable between os.
  If the program does not exist in the PATH then 'None' is returned.
  """
  def is_exe(fpath):
    return os.path.exists(fpath) and os.access(fpath, os.X_OK)

  fpath, fname = os.path.split(program)
  if fpath != '':
    if is_exe(program):
      return os.path.abspath(program)
  else:
    for path in [os.path.dirname(os.path.abspath(os.curdir)), os.curdir]:
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return os.path.abspath(exe_file)
    for path in os.environ['PATH'].split(os.pathsep):
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return os.path.abspath(exe_file)
  return None


def runCommands(cmds, localTempDir, inPipes=None, outPipes=None, errPipes=None):
  """ Run commands from CMDS list.
  """
  if inPipes is None:
    inPipes = [None] * len(cmds)
  if outPipes is None:
    outPipes = [None] * len(cmds)
  if errPipes is None:
    errPipes = [None] * len(cmds)
  for i, c in enumerate(cmds, 0):
    assert c is not None
    if inPipes[i] is None:
      sin = None
    else:
      sin = subprocess.PIPE
    if outPipes[i] is None:
      sout = None
    else:
      sout = subprocess.PIPE
    if errPipes[i] is None:
      serr = None
    else:
      serr = subprocess.PIPE
    p = subprocess.Popen(c, cwd=localTempDir, stdin=sin,
                         stdout=sout, stderr=serr)
    if inPipes[i] is None:
      sin = None
    else:
      if not os.path.exists(inPipes[i]):
        raise IOError('Unable to locate inPipe file: %s for command %s'
                      % (inPipes[i], ' '.join(c)))
      sin = open(inPipes[i], 'r').read()
    if outPipes[i] is None:
      pout, perr = p.communicate(sin)
      handleReturnCode(p.returncode, cmds[i])
    else:
      with open(outPipes[i], 'w') as f:
        f.write(p.communicate(sin)[0])
      handleReturnCode(p.returncode, cmds[i])


def handleReturnCode(retcode, cmd):
  """ handle the return codes from runCommands
  """
  if not isinstance(retcode, int):
    raise TypeError('handleReturnCode takes an integer for '
                    'retcode, not a %s.' % retcode.__class__)
  if not isinstance(cmd, list):
    raise TypeError('handleReturnCode takes a list for '
                    'cmd, not a %s.' % cmd.__class__)
  if retcode:
    if retcode < 0:
      raise RuntimeError('Experienced an error while trying to execute: '
                         '%s SIGNAL:%d' %(' '.join(cmd), -retcode))
    else:
      raise RuntimeError('Experienced an error while trying to execute: '
                         '%s retcode:%d' %(' '.join(cmd), retcode))


def populateExample(example, dirpath):
  """ Given an example, expand it into dirpath
  """
  for filename, data in example.files:
    if not os.path.exists(os.path.join(dirpath, 'results')):
      os.mkdir(os.path.join(dirpath, 'results'))
    with open(os.path.join(dirpath, 'results', filename), 'w') as f:
      f.write(data)


class ClumppPriorPopInfo(unittest.TestCase):
  def test_failure(self):
    """ Harvester should raise exception if clumpp run on prior pop info files.
    """
    tmpDirParent = makeTempDirParent()
    tmpDirs = []
    for example in td.clumppPriorPopInfo:
      tmpDirs.append(
        os.path.abspath(makeTempDir('clumppPriorPopInfo_%s' % example.name)))
      results = os.path.join(tmpDirs[-1], 'results')
      out = os.path.join(tmpDirs[-1], 'out')
      if not os.path.exists(out):
        os.mkdir(out)
      populateExample(example, tmpDirs[-1])
      files = glob(os.path.join(results, '*_f'))
      for f in files:
        data = hc.Data()
        data.records = {} # key is K, value is an array
        run, errorString = hc.readFile(f, data)
        self.assertTrue(run is not None)
        data.records.setdefault(run.k, []).append(run)
      data.sortedKs = data.records.keys()
      data.sortedKs.sort()
      self.assertRaises(hc.ClumppPriorPopInfo, hc.clumppGeneration,
                        results, out, data)
    [self.addCleanup(removeDir, d) for d in tmpDirs]
    self.addCleanup(removeDirIfEmpty, tmpDirParent)


class ClumppRegExFailure(unittest.TestCase):
  def test_failure(self):
    """ Harvester should raise exception if clumpp regular expressions fail
    """
    tmpDirParent = makeTempDirParent()
    tmpDirs = []
    for example in td.clumppRegEx:
      tmpDirs.append(
        os.path.abspath(makeTempDir('clumppRegEx_%s' % example.name)))
      results = os.path.join(tmpDirs[-1], 'results')
      out = os.path.join(tmpDirs[-1], 'out')
      if not os.path.exists(out):
        os.mkdir(out)
      populateExample(example, tmpDirs[-1])
      files = glob(os.path.join(results, '*_f'))
      for f in files:
        data = hc.Data()
        data.records = {} # key is K, value is an array
        run, errorString = hc.readFile(f, data)
        self.assertTrue(run is not None)
        data.records.setdefault(run.k, []).append(run)
      data.sortedKs = data.records.keys()
      data.sortedKs.sort()
      self.assertRaises(hc.ClumppRegEx, hc.clumppGeneration,
                        results, out, data)
    [self.addCleanup(removeDir, d) for d in tmpDirs]
    self.addCleanup(removeDirIfEmpty, tmpDirParent)


class ClumppLineStructureFailure(unittest.TestCase):
  def test_failure(self):
    """ Harvester should raise exception if clumpp generation sees weird line structure.
    """
    tmpDirParent = makeTempDirParent()
    tmpDirs = []
    for example in td.clumppLineStructure:
      tmpDirs.append(
        os.path.abspath(makeTempDir('clumppLineStructure_%s' % example.name)))
      results = os.path.join(tmpDirs[-1], 'results')
      out = os.path.join(tmpDirs[-1], 'out')
      if not os.path.exists(out):
        os.mkdir(out)
      populateExample(example, tmpDirs[-1])
      files = glob(os.path.join(results, '*_f'))
      for f in files:
        data = hc.Data()
        data.records = {} # key is K, value is an array
        run, errorString = hc.readFile(f, data)
        self.assertTrue(run is not None)
        data.records.setdefault(run.k, []).append(run)
      data.sortedKs = data.records.keys()
      data.sortedKs.sort()
      self.assertRaises(hc.ClumppLineStructure, hc.clumppPopFile,
                        results, out, data)
    [self.addCleanup(removeDir, d) for d in tmpDirs]
    self.addCleanup(removeDirIfEmpty, tmpDirParent)


class UnexpectedValue(unittest.TestCase):
  def test_failure(self):
    """ Harvester should raise exception if the input contains unexpected values.
    """
    tmpDirParent = makeTempDirParent()
    tmpDirs = []
    for example in td.readFileUnexpectedValue:
      tmpDirs.append(
        os.path.abspath(makeTempDir('readFileUnexpected_%s' % example.name)))
      results = os.path.join(tmpDirs[-1], 'results')
      out = os.path.join(tmpDirs[-1], 'out')
      if not os.path.exists(out):
        os.mkdir(out)
      populateExample(example, tmpDirs[-1])
      files = glob(os.path.join(results, '*_f'))
      for f in files:
        data = hc.Data()
        data.records = {} # key is K, value is an array
        self.assertRaises(hc.UnexpectedValue, hc.readFile, f, data)
    [self.addCleanup(removeDir, d) for d in tmpDirs]
    self.addCleanup(removeDirIfEmpty, tmpDirParent)


class KnownGood(unittest.TestCase):
  def test_fullRun(self):
    """ HarvesterCore should read and process these examples without error.
    """
    tmpDirParent = makeTempDirParent()
    tmpDirs = []
    for example in td.knownGood:
      tmpDirs.append(os.path.abspath(makeTempDir('readFile_%s' % example.name)))
      populateExample(example, tmpDirs[-1])
      cmd = [which('structureHarvester.py'),
             '--dir', os.path.join(tmpDirs[-1], 'results'),
             '--out', os.path.join(tmpDirs[-1], 'out'),
             '--evanno', '--clumpp']
      runCommands([cmd], tmpDirs[-1])
    [self.addCleanup(removeDir, d) for d in tmpDirs]
    self.assertTrue(True)
    self.addCleanup(removeDirIfEmpty, tmpDirParent)



if __name__ == '__main__':
  unittest.main()
