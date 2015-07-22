#! /usr/bin/env python

####################################
# admix_rel.py
# written by Raymond Walters, July 2015
"""
Estimates relatedness for admixed sample using ADMIXTURE, REAP
"""
# Overview:
# 1) Input source and target plink files
#    - assume matching QCed LD pruned SNPs, source is unrelated subset of target IDs
# 2) Run ADMIXTURE on unrelated set
# 3) Select population exemplars based on admixture proportions
# 4) Run ADMIXTURE on full data in supervised mode using pop exemplars
# 5) Estimate relatedness with REAP
# 6) Generate diagnostic plots
#    - exemplars on PCA (if PCA available)
#    - final admixture on PCA (if PCA available)
#    - IBD0/IBD1 for REAP
#
####################################


####################################
# Setup
# a) load python dependencies
# b) get variables/arguments
# c) read config file
# d) check dependencies
####################################

import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import os
import subprocess
from distutils import spawn
import argparse
from py_helpers import unbuffer_stdout
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='admix_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))


# import

# argparse

# check dependencies

# run admix

# select exemplars in admix output (error out if a pop has < 10 exemplars)

# run supervised admix

# run reap using admix results

# Generate diagnostic plots
# - exemplars on PCA (if PCA available)
# - final admixture on PCA (if PCA available)
# - IBD0/IBD1 for REAP

# final cleanup

# finish
