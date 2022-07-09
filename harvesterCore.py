""" harvesterCore
dent earl, dearl (a) soe ucsc edu
March 2012

Functions shared between the web and stand-alone versions
of Structure Harvester.

http://www.structureharvester.com/
https://github.com/dentearl/structureHarvester/
http://taylor0.biology.ucla.edu/structureHarvester/
http://users.soe.ucsc.edu/~dearl/software/structureHarvester/


##############################
CITATION
Earl, Dent A. and vonHoldt, Bridgett M. (2012)
STRUCTURE HARVESTER: a website and program for visualizing
STRUCTURE output and implementing the Evanno method.
Conservation Genetics Resources 4(2) 359-361. DOI: 10.1007/s12686-011-9548-7

##############################
REFERENCES

Evanno et al., 2005.  Detecting the number of clusters of individuals using
the software STRUCTURE: a simulation study. Molecular Ecology 14, 2611-2620.

Jakobsson M., Rosenberg N. 2007. CLUMPP: a cluster matching and permutation
program for dealing with label switching and multimodality in analysis
of population structure. Bioinformatics 23(14): 1801-1806.

Pritchard J., Stephens M., Donnelly. P. 2000. Genetics 155:945-959.

##############################
LICENSE

Copyright (C) 2007-2022 by
Dent Earl (dearl (a) soe ucsc edu, dentearl (a) gmail com)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""
import math
import os
import re
import time

__version__ = 'vB.1 July 2022' # alpha.number convention for core
EPSILON = 0.0000001 # for determining if a stdev ~ 0


class HarvesterError(Exception): pass
class HarvesterErrorDetail(HarvesterError):
  def __init__(self, filename, valuename, value, data):
    self.filename = filename
    self.valuename = valuename
    self.value = value
    self.data = data
class PriorPopulation(HarvesterError): pass
class UnexpectedValue(HarvesterErrorDetail):
  def __str__(self):
    return repr('Error: One of your files, %s, contains an unexpected value:'
                ' [%s = %s] '
                'Generally these problems can be resolved by discarding the '
                'file and re-running STRUCTURE for this value of K.' %
                (self.filename, self.valuename, self.value))
class ClumppRegEx(HarvesterError):
  def __init__(self, filename, regexs, lineno, line, data):
    self.filename = filename
    self.regexs = regexs
    self.lineno = lineno
    self.line = line
    self.data = data
  def __str__(self):
    expressions = ',\n'.join(self.regexs)
    return repr('Error, clumpp generation for file %s failed to '
                'match any of the regular expressions %s '
                'line number:%d %s '
                'Please verify your input file\'s characters and '
                'try again after making a correction.' %
                (self.filename, expressions, self.lineno, self.line))
class ClumppLineStructure(HarvesterError):
  def __init__(self, filename, regex, lineno, line, data):
    self.filename = filename
    self.regex = regex
    self.lineno = lineno
    self.line = line
    self.data = data
  def __str__(self):
    return repr('Error, clumpp popfile generation failed. Your file %s '
                'has an unexepected line structure! Job failed to match '
                'the regular expression: %s on line number %d: %s' %
                (self.filename, self.regex, self.lineno, self.line))
class ClumppPriorPopInfo(HarvesterError):
  def __init__(self, filename ,data):
    self.filename = filename
    self.data = data
  def __str__(self):
    return repr('Error, clumpp generation for file %s failed: '
                'STRUCTURE was run using prior population information '
                'and it is not clear how to convert his into a stand '
                'alone Q matrix.' % self.filename)


class Data:
  """ Data will be used to store variables in a fashion
  that allows for easy access between functions.
  """
  def __init__(self):
    self.evanno          = False
    self.uniqueName      = ''
    self.records         = None
    self.sortedKs        = None
    self.localURL        = None
    self.homeURL         = None
    self.estLnProbMeans  = None
    self.estLnProbStdevs = None
    self.LnPK            = None
    self.LnPPK           = None
    self.deltaK          = None


class RunRecord:
  """ Stores the output of a single structure run.
  """
  def __init__(self):
    self.name      = ''
    self.k         = -1
    self.indivs    = -1
    self.loci      = -1
    self.burnin    = -1
    self.reps      = -1
    self.estLnProb = -1
    self.meanLlh   = -1
    self.varLlh    = -1
    self.runNumber = -1


def addAttribute(p, value, run, data):
  """ add a particular attribute to the run object. p is a string with the type
  value is a string that contains the value to add, run is the run object, data
  is the data object and badValue is a function that takes the run name, the
  bad value's string name, the value, and the data object.
  """
  if p == 'indivs':
    try:
      run.indivs = int(value)
    except ValueError:
      raise UnexpectedValue(run.name, 'individuals', value, data)
  elif p == 'loci':
    try:
      run.loci = int(value)
    except ValueError:
      raise UnexpectedValue(run.name, 'loci', value, data)
  elif p == 'k':
    try:
      run.k = int(value)
    except ValueError:
      raise UnexpectedValue(run.name, 'populations assumed', value, data)
  elif p == 'burnin':
    try:
      run.burnin = int(value)
    except ValueError:
      raise UnexpectedValue(run.name, 'Burn-in period', value, data)
  elif p == 'reps':
    try:
      run.reps = int(value)
    except ValueError:
      raise UnexpectedValue(run.name, 'Reps', value, data)
  elif p == 'lnprob':
    if value == 'nan':
      raise UnexpectedValue(run.name, 'Estimated Ln Prob of Data', value, data)
    try:
      run.estLnProb = float(value)
    except ValueError:
      raise UnexpectedValue(run.name, 'Estimated Ln Prob of Data', value, data)
  elif p == 'meanln':
    if value == 'meanln':
      raise UnexpectedValue(run.name, 'Estimated Ln Prob of Data', value, data)
    try:
      run.meanLlh = float(value)
    except ValueError:
      raise UnexpectedValue(run.name, 'Mean value of ln likelihood', value, data)
  elif p == 'varln':
    if value == 'nan':
      raise UnexpectedValue(run.name, 'Estimated Ln Prob of Data', value, data)
    try:
      run.varLlh = float(value)
    except ValueError:
      raise UnexpectedValue(run.name, 'Variance of ln likelihood', value, data)
  else:
    sys.stderr.write('Error, %s unknown pattern type %s\n'
                     % (data.uniqueName, p))


def validateRecord(run):
  """ validateRecord checks the run object to make sure it contains all
  necessary information. If the run is valid it returns the run,
  else it returns None.
  """
  errorString = ''
  if run.name == '':
    return None
  dataNames = ['k', 'indivs', 'loci', 'burnin', 'reps',
               'estLnProb', 'meanLlh', 'varLlh']
  for i, r in enumerate([run.k, run.indivs, run.loci, run.burnin, run.reps,
               run.estLnProb, run.meanLlh, run.varLlh]):
    if r == -1:
      errorString = ('StructureHarvester: In file %s, unable to read value '
                     'for %s\n'
                     'This can be caused by forgetting to check the box '
                     '"Compute the '
                     'probability of the data (for estimating K)" in '
                     'STRUCTURE.\n'
                     % (run.name, dataNames[i]))
      return None, errorString
  return run, errorString


def writeRawOutputToFile(filename, data):
  file = open(filename, 'w')
  file.write('# This document produced by structureHarvester core version %s\n'
             % __version__)
  file.write('# http://www.structureharvester.com/\n')
  file.write('# https://github.com/dentearl/structureHarvester/\n')
  file.write('# http://taylor0.biology.ucla.edu/structureHarvester\n')
  file.write('# http://users.soe.ucsc.edu/~dearl/software/structureHarvester\n')
  file.write('# Written by Dent Earl, dearl (a) soe ucsc edu.\n')
  file.write('# CITATION:\n# Earl, Dent A. and vonHoldt, Bridgett M. (2012)\n'
             '# STRUCTURE HARVESTER: a website and program for visualizing\n'
             '# STRUCTURE output and implementing the Evanno method.\n'
             '# Conservation Genetics Resources 4(2) 359-361. '
             'DOI: 10.1007/s12686-011-9548-7\n'
             '# Core version: %s\n'
             % __version__)
  file.write('# File generated at %s\n'
             % (time.strftime('%Y-%b-%d %H:%M:%S %Z', time.localtime())))
  file.write('#\n')
  file.write('\n##########\n')
  file.write('# K\tReps\t'
        'mean est. LnP(Data)\t'
        'stdev est. LnP(Data)\n')
  for i in range(0, len(data.sortedKs)):
    k = data.sortedKs[i]
    file.write('%d\t%d\t%f\t%f\n' % (k, len(data.records[k]),
                                     data.estLnProbMeans[k],
                                     data.estLnProbStdevs[k]))
  file.write('\n##########\n')
  file.write('# File name\tRun #\t'
         'K\tEst. Ln prob. of data\t'
         'Mean value of Ln likelihood\t'
         'Variance of Ln likelihood\n')
  for k in data.records:
    for r in data.records[k]:
      if r.runNumber != -1:
        numStr = r.runNumber
      else:
        numStr = ''
      file.write('%s\t%s\t%d\t%.1f\t%.1f\t%.1f\n'
                 % (r.name, numStr, r.k, r.estLnProb, r.meanLlh, r.varLlh))
  file.write('\nSoftware written by Dent Earl while at EEB Dept, UCLA, '
             'BME Dept UCSC\n'
             'dearl (a) soe ucsc edu\n')
  file.write('CITATION:\nEarl, Dent A. and vonHoldt, Bridgett M. (2012)\n'
             'STRUCTURE HARVESTER: a website and program for visualizing\n'
             'STRUCTURE output and implementing the Evanno method.\n'
             '# Conservation Genetics Resources 4(2) 359-361. '
             'DOI: 10.1007/s12686-011-9548-7\n'
             'Core version: %s' % __version__)
  file.close()


def readFile(filename, data):
  run = RunRecord()
  run.name = os.path.basename(filename)
  m = re.search(r'.*_(\d+)_f$', run.name)
  if m != None:
    run.runNumber = m.group(1)
  regExs = {'indivs' : r'^(\d+) individuals.*$',
              'loci'   : r'^(\d+) loci.*$',
              'k'      : r'^(\d+) populations assumed.*$',
              'burnin' : r'^(\d+) Burn\-in period.*$',
              'reps'   : r'^(\d+) Reps.*$',
              # nan, inf
              'lnprob' : r'^Estimated Ln Prob of Data\s+=\s+([\d.$nainf-]+).*$',
              'meanln' : r'^Mean value of ln likelihood\s+=\s+([\d.$nainf-]+).*$',
              'varln'  : r'^Variance of ln likelihood\s+=\s+([\d.$nainf-]+).*$',
            }
  pats = {}
  for r in regExs:
    pats[r] = re.compile(regExs[r])
  isDone = False
  with open(filename, 'r') as infile:
    for line in infile:
      if isDone:
        break
      line = line.strip()
      for p in pats:
        m = re.match(pats[p], line)
        if m != None:
          addAttribute(p, m.group(1), run, data)
          if p == 'varln':
            isDone = True
  return validateRecord(run)


def clumppGeneration(inputDir, outputDir, data):
  regex1 = r'''^(\d+)                  # individual number
         \s+(\S+)          # Label
         \s+(\S+)          # %Miss
         \s+([\-\d]+)\s+:  # Pop, can be -9
         \s+(.*?)$         # Inferred clusters'''
  regex2 = r'''^(\d+)              # individual number
         \s+(\S+)      # Label
         \s+(\S+)\s+:  # %Miss. if only 1 inferred cluster, then 1 is not printed
         \s+(.*?)$     # Inferred clusters'''
  # regex2 also covers the new style of confidence intervals
  # following the standard values
  regex3 = r'''^(\d+)                # individual number
         \s+(\S+)\s+:    # %Miss
         \s+([\d. ]*?)$  # Inferred clusters'''
  regexIntervals = r'\(\S+,\S+\)' # used to trim off the confidence intervals
  pats = []
  for r in [regex1, regex2, regex3]:
      pats.append(re.compile(r, re.X))
  patInterval = re.compile(regexIntervals)
  for k in data.sortedKs:
    with open(os.path.join(outputDir, 'K%d.indfile' % k), 'w') as f:
      for r in data.records[k]:
        with open(os.path.join(inputDir, r.name), 'r') as inFile:
          printing = False
          for lineno, line in enumerate(inFile, 1):
            line = line.strip()
            if line == '':
              printing = False
              continue
            if line.startswith('Probability of being from assumed population | '):
              raise ClumppPriorPopInfo(r.name, data)
            d = line.split()
            if (d[0] == 'Label' and d[1] == '(%Miss)') or d[0] == '(%Miss)':
              printing = True
              continue
            # find all confidence intervals and remove them
            m = re.findall(patInterval, line)
            for i in m:
              line = line.replace(i, '')
            line = line.strip()
            if printing:
              m1 = re.match(pats[0], line)
              m2 = re.match(pats[1], line)
              m3 = re.match(pats[2], line)
              if m1 is None and m2 is None and m3 is None:
                # We failed to match any of the Clumpp Regular Exceptions
                raise ClumppRegEx(r.name, [regex1, regex2, regex3],
                                  lineno, line, data)
              else:
                if m1 != None:
                  f.write('%3d %3d %s %2d : %s\n'
                          % (int(m1.group(1)), int(m1.group(1)), m1.group(3),
                             int(m1.group(4)), m1.group(5)))
                elif m2 != None:
                  f.write('%3d %3d %s 1 : %s\n'
                          % (int(m2.group(1)), int(m2.group(1)),
                             m2.group(3), m2.group(4)))
                elif m3 != None:
                  f.write('%3d %3d %s 1 : %s\n'
                          % (int(m3.group(1)), int(m3.group(1)),
                             m3.group(2), m3.group(3)))
          f.write('\n')


def clumppPopFile(inputDir, outputDir, data):
  regex1 = r'''^(\d+): # Given Pop
         \s+([\d. ]*?) # All Inferred Cluster columns
         \s+(\d+)$ # # Number of Individuals'''
  pats = []
  f = None
  for r in [regex1]:
    pats.append(re.compile(r, re.X))
  for k in data.sortedKs:
    for r in data.records[k]:
      with open(os.path.join(inputDir, r.name), 'r') as inFile:
        printing = False
        skipLine = False
        workComplete = False
        for lineno, line in enumerate(inFile, 1):
          if workComplete:
            continue
          line = line.strip()
          d = line.split()
          if printing and line.startswith('--------------------'):
            workComplete = True
            printing = False
            continue
          if len(d) < 3:
            continue
          if d[0] == 'Given' and d[1] == 'Inferred' and d[2] == 'Clusters':
            printing = True
            skipLine = True
            maybeCloseF(f)
            f = open(os.path.join(outputDir, 'K%d.popfile' % k), 'a')
            continue
          if skipLine:
            skipLine = False
            continue
          if printing:
            m1 = re.match(pats[0], line)
            if m1 is None:
              maybeCloseF(f)
              raise ClumppLineStructure(r.name, regex1, lineno, line, data)
            else:
              f.write('%3d:' % (int(m1.group(1))))
              for i in range(2, len(m1.groups()) + 1):
                f.write('\t%s' % m1.group(i))
              f.write('\n')
        maybeCloseF(f)


def maybeCloseF(f):
  """ Helper function to try to clean up calls in clumppPopFile().
  """
  if f != None:
    if not f.closed:
      f.write('\n')
      f.close()


def evannoTests(data, isWeb=False):
  """ sucess returns None, failure returns
  a string which can be directly printed in to the output.
  """
  fail = False
  numRuns = 0
  for k in data.records:
    numRuns += len(data.records[k])
  numReps = numRuns / float(len(data.records))
  out  = '<p>Stats: number of runs = %d, ' % numRuns
  out += 'number of replicates = %.2f, ' % numReps
  out += 'minimum K tested = %d, ' % data.sortedKs[0]
  out += 'maximum K tested = %d.</p>\n' % data.sortedKs[-1]
  # Must test at least three values of K
  if len(data.sortedKs) >= 3:
    out += '<p style="color:gray;">Test: '
    out += 'You must test at least 3 values of K. PASS</p>\n'
  else:
    fail = True
    out += '<p style="color:red;">Test: '
    out += 'You must test at least 3 values of K. FAIL</p>\n'

  # K values must be sequential
  if ((data.sortedKs[-1] + 1) - data.sortedKs[0]) == len(data.sortedKs):
    out += '<p style="color:gray;">Test: '
    out += 'K values must be sequential. PASS</p>\n'
  else:
    fail = True
    out += '<p style="color:red;">Test: '
    out += 'K values must be sequential. FAIL</p>\n'
    out += '<p style="margin-left:5em;">'
    prevK = data.sortedKs[0] - 1
    i = 0
    for k in data.sortedKs:
      i += 1
      if i == len(data.sortedKs):
        spacer = ''
      else:
        spacer = ', '
      if k != prevK + 1:
        out += '<span style="color:red;">missing: %d</span>%s' % (k, spacer)
      else:
        out += '%d%s' % (k, spacer)
      prevK = k
    out += '</p>\n'
  # Number of replicates must be greater than 2 (stdev)
  if numReps > 1:
    out += '<p style="color:gray;">Test: '
    out += 'The number of replicates per K > 1. PASS</p>\n'
  else:
    fail = True
    out += '<p style="color:red;">Test: '
    out += 'The number of replicates per K > 1. FAIL</p>\n'
  # Standard Devation for a K (but not the first or last K) is zero
  for i in range(1, len(data.sortedKs) - 1):
    k = data.sortedKs[i]
    if data.estLnProbStdevs[k] < EPSILON: # our epsilon
      fail = True
      out += ('<p style="color:red;">Test: Standard devation of est. '
              'Ln Pr(Data) '
              'is less than 0.0000001. FAIL for K = %d. The Evanno test '
              'requires division '
              'by the standard deviation of the est. Ln Pr(Data) values for '
              'all K between '
              'the first and last K value, and thus when the standard '
              'deviation is within '
              'Epsilon (%f) of zero we cannot perform the test. (Try '
              'running more '
              'replicates to hopefully increase the standard deviation for '
              'that value of '
              ' K.)</p>\n' % (k, EPSILON))
  if fail:
    if isWeb:
      return out
    else:
      re.sub('&[^&]+?;', '', out)
      return re.sub('<[^<]+?>', '', out)
  return None


def calculateMeansAndSds(data):
  data.estLnProbMeans = {}
  data.estLnProbStdevs = {}
  for k in data.records:
    data.estLnProbMeans[k] = 0
    data.estLnProbStdevs[k] = 0
    for r in data.records[k]:
      data.estLnProbMeans[k] += r.estLnProb
    data.estLnProbMeans[k] /= len(data.records[k])
    if len(data.records[k]) > 1:
      for r in data.records[k]:
        data.estLnProbStdevs[k] += (r.estLnProb - data.estLnProbMeans[k]) ** 2
      data.estLnProbStdevs[k] /= (len(data.records[k]) - 1)
      data.estLnProbStdevs[k] = math.sqrt(data.estLnProbStdevs[k])


def calculatePrimesDoublePrimesDeltaK(data):
  """ This function takes in the data object and uses the
  estimated log probability means dictionary (data.estLnProbMeans) and the
  estimated log probability standard deviations dictionary
  (data.estLnProbStdevs) to calculate dictionaries keyed on K values (ints)
  for the three Evanno quantities of interest:
  L'(K) : data.LnPK
  L''(K) : data.LnPP(K)
  delta K : data.deltaK
  Note that to calculate the deltaK for 'thisK' you need estimated log prob mean
  values for both the previous K, 'prevK' and the next K, 'nextK'. So if you run
  Structure for K = 1..20, you'll only get delta K for K = 2..19.
  """
  data.LnPK = {}
  data.LnPPK = {}
  data.deltaK = {}
  for i in range(1, len(data.sortedKs)):
    thisK = data.sortedKs[i]
    prevK = data.sortedKs[i - 1]
    data.LnPK[thisK] = data.estLnProbMeans[thisK] - data.estLnProbMeans[prevK]
  for i in range(1, len(data.sortedKs) - 1):
    prevK = data.sortedKs[i - 1]
    thisK = data.sortedKs[i]
    nextK = data.sortedKs[i + 1]
    data.LnPPK[thisK] = abs(data.LnPK[nextK] - data.LnPK[thisK])
    # data.deltaK[thisK] = data.LnPPK[thisK] / float(data.estLnProbStdevs[thisK])
    data.deltaK[thisK] = (abs(data.estLnProbMeans[nextK] -
                              2.0 * data.estLnProbMeans[thisK] +
                              data.estLnProbMeans[prevK]) /
                          float(data.estLnProbStdevs[thisK]))
