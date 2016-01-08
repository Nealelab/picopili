#! /usr/bin/env python

####################################
# impute_rel.py
# written by Raymond Walters, January 2016
"""
Impute GWAS data with related individuals
"""
# 
# Is a wrapper for:
# - shape_rel.py
# - imp2_rel.py
# - bg_imp.py
# - agg_imp.py
# 
####################################



import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import os
import subprocess
import warnings
from args_impute import *
from py_helpers import unbuffer_stdout, read_conf, file_tail, link, warn_format
unbuffer_stdout()
warnings.formatwarning = warn_format

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############
parser = argparse.ArgumentParser(prog='impute_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parserphase, parserimpute, parserchunk, parserref, parserbg, parsercluster])
                    
args = parser.parse_args()


# TODO: print args




#############
print '\n...Checking dependencies...'
#############




# TODO: here



#############
print '\n...Submitting first task...'
#############




# TODO: here




# finish
print '\n############'
print '\n'
print 'SUCCESS!'
print 'All jobs submitted.\n'
exit(0)


# eof