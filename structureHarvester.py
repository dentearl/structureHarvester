#!/usr/bin/env python2.6
""" structureHarvester
2007-2012
dent earl, dearl (a) soe ucsc edu

http://users.soe.ucsc.edu/~dearl/software/structureHarvester/
http://taylor0.biology.ucla.edu/structureHarvester/

##############################
CITATION
Earl, Dent A. and vonHoldt, Bridgett M. (2011)                
STRUCTURE HARVESTER: a website and program for visualizing    
STRUCTURE output and implementing the Evanno method.          
Conservation Genetics Resources DOI: 10.1007/s12686-011-9548-7

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

Copyright (C) 2007-2012 by 
Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)

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
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

"""
import glob
import harvesterCore as hc
from optparse import OptionParser
import os
import re
import sys
import time

__version__ = 'v0.6.92 March 2012'
EPSILON = 0.0000001 # for determining if a stdev ~ 0

try:
    eval("enumerate([1, 2, 3], 1)")
except TypeError:
    raise ImportError('This script requires Python version 2.5 <= v < 3.0. '
                      'This version %d.%d' % (sys.version_info[0], sys.version_info[1]))

def initOptions(parser):
   parser.add_option('--dir',dest='resultsDir',
                     type ='string', help='The structure Results/ directory.')
   parser.add_option('--out',dest='outDir',
                     type ='string', help='The out directory. If it does not exist, it will be created. Output written to summary.txt')
   parser.add_option('--evanno',dest='evanno', action='store_true', default=False,
                     help='If possible, performs the Evanno 2005 method. Written to evanno.txt. default=%default')
   parser.add_option('--clumpp',dest='clumpp', action='store_true', default=False,
                     help='Generates one K*.indfile for each value of K run, for use with CLUMPP. default=%default')

def checkOptions(parser, options):
   if not options.resultsDir:
       parser.error('Error, specify --dir \n')
   if not os.path.exists(options.resultsDir):
       parser.error('Error, --dir %s does not exist!\n' % options.resultsDir)
   if not os.path.isdir(options.resultsDir):
       parser.error('Error, --dir %s is not a directory!\n' % options.resultsDir)
   options.resultsDir = os.path.abspath(options.resultsDir)
   if not options.outDir:
       parser.error('Error, specify --out \n')
   if os.path.exists(options.outDir) and not os.path.isdir(options.outDir):
       parser.error('Error, --out %s already exists but is not a directory!\n' % options.outDir)
   if not os.path.exists(options.outDir):
       os.mkdir(options.outDir)

def badValue(filename, valuename, value, data):
    sys.stderr.write('Error: %s contains an unexpected value:\n'
                     '    %s = %s\n'
                     'Generally these problems can be resolved by discarding the file and '
                     're-running STRUCTURE for this value of K.'
                     % (filename, valuename, value))
    sys.exit(1)
   
def harvestFiles(data, options):
    files = glob.glob(os.path.join(options.resultsDir, '*_f'))
    if len(files) < 1:
        sys.stderr.write('Error, unable to locate any _f files in '
                         'the results directory %s' % options.resultsDir)
        sys.exit(1)
    data.records = {} # key is K, value is an array
    for f in files:
        run, errorString = hc.readFile(f, data, badValue)
        if run is not None:
            if run.k not in data.records:
                data.records[run.k] = []
            data.records[run.k].append(run)
        else:
            sys.stderr.write('Error, unable to extract results from file %s.\n' % f)
            sys.stderr.write('%s\n' % errorString)
            sys.exit(1)
    data.sortedKs = data.records.keys()
    data.sortedKs.sort()

def evannoMethod(data, options):
    if not options.evanno:
        return
    value = hc.evannoTests(data)
    if value is not None:
        sys.stderr.write('Unable to perform Evanno method for the follow reason(s):\n')
        sys.stderr.write(value)
        sys.exit(1)
    hc.calculatePrimesDoublePrimesDeltaK(data)
    writeEvannoTableToFile(data, options)

def writeEvannoTableToFile(data, options):
    file = open(os.path.join(options.outDir, 'evanno.txt'), 'w')
    file.write('# This document produced by structureHarvester.py %s core %s\n' % (__version__, hc.__version__))
    file.write('# http://users.soe.ucsc.edu/~dearl/struct_harvest\n')
    file.write('# http://taylor0.biology.ucla.edu/struct_harvest\n')
    file.write('# Written by Dent Earl, dearl (a) soe ucsc edu.\n')
    file.write('# CITATION:\n# Earl, Dent A. and vonHoldt, Bridgett M. (2011)\n'
               '# STRUCTURE HARVESTER: a website and program for visualizing\n'
               '# STRUCTURE output and implementing the Evanno method.\n'
               '# Conservation Genetics Resources DOI: 10.1007/s12686-011-9548-7\n'
               '# Stand-alone version: %s\n'
               '# Core version: %s\n'
               % (__version__, hc.__version__))
    file.write('# File generated at %s\n' % (time.strftime('%Y-%b-%d %H:%M:%S %Z', time.localtime())))
    file.write('#\n')
    file.write('\n##########\n')
    file.write('# K\tReps\t'
                'Mean LnP(K)\tStdev LnP(K)\t'
                'Ln\'(K)\t|Ln\'\'(K)|\tDelta K\n')
    for i in xrange(0, len(data.sortedKs)):
        k = data.sortedKs[i]
        if k in data.LnPK:
            LnPKstr = '%f' % data.LnPK[k]
        else:
            LnPKstr = 'NA'
        if k in data.LnPPK:
            LnPPKstr = '%f' % data.LnPPK[k]
        else:
            LnPPKstr = 'NA'
        if k in data.deltaK:
            deltaKstr = '%f' % data.deltaK[k]
        else:
            deltaKstr = 'NA'
        file.write('%d\t'
                    '%d\t%.4f\t%.4f\t'
                    '%s\t%s\t%s\n' % (k, len(data.records[k]),
                                      data.estLnProbMeans[k],
                                      data.estLnProbStdevs[k],
                                      LnPKstr, LnPPKstr, deltaKstr))
    file.close()

def failHandler(message):
    sys.stderr.write(message)
    sys.exit(1)

def main():
    usage = ('usage: %prog --dir=path/to/dir/ --out=path/to/dir/ [options]\n\n'
              '%prog takes a STRUCTURE results directory (--dir) and an\n'
              'output directory (--out will be created if it does not exist) and then\n'
              'depending on the other options selected harvests data from the results\n'
              'directory and performs the selected analyses')
    data = hc.Data()
    parser = OptionParser(usage=usage, version='%prog '+__version__+' core '+hc.__version__)
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(parser, options)

    harvestFiles(data, options)
    hc.calculateMeansAndSds(data)
    if options.clumpp:
        hc.clumppGeneration(options.resultsDir, options.outDir, data, failHandler)
        hc.clumppPopFile(options.resultsDir, options.outDir, data, failHandler)
    evannoMethod(data, options)
    hc.writeRawOutputToFile(os.path.join(options.outDir, 'summary.txt'), data)

if __name__ == '__main__':
    main()
