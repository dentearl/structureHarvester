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

## References
* J. Pritchard, M. Stephens, P. Donnelly. 2000. Genetics 155:945-959. http://pritch.bsd.uchicago.edu/structure.html http://www.genetics.org/cgi/content/full/155/2/945
* Evanno et al., 2005. Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study. Molecular Ecology 14 , 2611 - 2620 http://www3.interscience.wiley.com/journal/118706173/abstract
* M. Jakobsson, N. Rosenberg 2007. CLUMPP: a cluster matching and permutation program for dealing with label switching and multimodality in analysis of population structure. Bioinformatics 23(14): 1801-1806. http://rosenberglab.bioinformatics.med.umich.edu/clumpp.html http://bioinformatics.oxfordjournals.org/cgi/content/full/23/14/1801

## Citation
* Earl, Dent A. and vonHoldt, Bridgett M. (2012) STRUCTURE HARVESTER: a website and program for visualizing STRUCTURE output and implementing the Evanno method. Conservation Genetics Resources vol. 4 (2) pp. 359-361. doi: 10.1007/s12686-011-9548-7 http://www.springerlink.com/content/jnn011511h415358/
