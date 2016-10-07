### Python Dependencies

Picopili is built from a combination of Python, Perl, R, and *nix shell scripts.

Most scripts depend only on packages from the Python Standard Library
([[https://docs.python.org/2/library/]]). In addition, `admix_rel.py` depends
on numpy. We strongly support using Anaconda ([[https://www.continuum.io/downloads]])
to manage Python package dependencies, but a barebones installation of Python 2.X + numpy
should be sufficient for picopili in most cases.

Scripts are primarily tested under Python 2.7 and Anaconda 2.1.0, but should be broadly
compatible with most Python 2.X versions. If you encounter compaitibility issues,
please contact rwalters(at)broadinstitute.org and we would be happy to assist.

##### Full list of package dependencies

* argparse 
* cPickle 
* copy 
* distutils 
* glob 
* gzip 
* math 
* numpy
* os 
* random 
* re 
* string 
* subprocess 
* sys 
* textwrap 
* time 
* warnings 

