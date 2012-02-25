# Structure Harvester
(c) 2007-2012 The Author, see LICENSE for details.

## Author
[Dent Earl](https://github.com/dentearl/)

## Dependencies
* Python version 2.5 < x < 3.0

## Installation
1. Download or <code>clone</code> the package.
2. Place <code>structureHarvester.py</code> and <code>harvesterCore.py</code> in the same directory (or leave them where they are).
3. There is no step 3.

## Usage
    Usage: structureHarvester.py --dir=path/to/dir/ --out=path/to/dir/ [options]

    structureHarvester.py takes a STRUCTURE results directory (--dir) and an
    output directory (--out will be created if it does not exist) and then
    depending on the other options selected harvests data from the results
    directory and performs the selected analyses

    Options:
      --version         show program's version number and exit
      -h, --help        show this help message and exit
      --dir=RESULTSDIR  The structure Results/ directory.
      --out=OUTDIR      The out directory. If it does not exist, it will be created. Output written to summary.txt
      --evanno          If possible, performs the Evanno 2005 method. Written to evanno.txt. default=False
      --clumpp          Generates one K*.indfile for each value of K run, for use with CLUMPP. default=False
