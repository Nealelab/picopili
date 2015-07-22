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
from py_helpers import read_conf, unbuffer_stdout
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='admix_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))

arg_base = parser.add_argument_group('Basic Arguments')
arg_admix = parser.add_argument_group('Admixture Settings')
arg_reap = parser.add_argument_group('Relatedness Settings')
arg_exloc = parser.add_argument_group('Software Executable Locations')

arg_base.add_argument('--unrel-bfile', 
                    type=str,
                    metavar='FILESTEM',
                    help='File stem for plink bed/bim/fam files ' + \
                         'with unrelated individuals to estimate admixture.',
                    required=True)
arg_base.add_argument('--target-bfile', 
                    type=str,
                    metavar='FILESTEM',
                    help='file stem for plink bed/bim/fam files. ' + \
                         'Relatedness will be estimated for these samples. ' + \
                         'All individuals from --unrel-bfile should also be present in this data.',
                    required=True)
arg_base.add_argument('--outname',
                    type=str,
                    metavar='OUTNAME',
                    help='base name for output files; recommend 4 character stem to match ricopili',
                    required=True)
arg_base.add_argument('--outdir',
                    type=str,
                    metavar='DIRNAME',
                    help='Directory for output files. Will create if needed. ' + \
                    'Uses ./OUTNAME_admix_rel if unspecified',
                    required=False)
arg_base.add_argument('--no-cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')

arg_admix.add_argument('--npops',
                    type=int,
                    metavar='INT',
                    help='Number of ancestral populations for admixture',
                    required=False,
                    default=4)
arg_admix.add_argument('--prop-th',
                    type=float,
                    metavar='FLOAT',
                    help='Minimum admixture proportion to select an individual ' + \
                         'from the unrelated set as an exemplar for the ancestral population',
                    required=False,
                    default=.95)
arg_admix.add_argument('--min-exemplar',
                    type=int,
                    metavar='INT',
                    help='Minimum number of individuals from unrelated set ' + \
                         'assigned as exemplars for each ancestral population',
                    required=False,
                    default=20)
arg_admix.add_argument('--multithread-cores',
                       type=int,
                       metvar='INT',
                       help='Number of cores to use for multi-threading in admixture analysis',
                       required=False,
                       default=1)

arg_reap.add_argument('--min-rel',
                      type=float,
                      metavar='FLOAT',
                      help='Minimum pi-hat relatedness level to include in output. ' + \
                           'Default is halfway between 3rd and 4th degree relatives.',
                      required=False,
                      default=.09375)

arg_exloc.add_argument('--rscript-ex',
                    type=str,
                    metavar='PATH',
                    help='path to Rscript executable, tries reading from PATH if unspecified',
                    required=False,
                    default=None)
arg_exloc.add_argument('--admixture-ex',
                    type=str,
                    metavar='PATH',
                    help='path to ADMIXTURE executable',
                    required=False,
                    default="/humgen/atgu1/fs03/shared_resources/shared_software/bin/admixture")
arg_exloc.add_argument('--reap-ex',
                    type=str,
                    metavar='PATH',
                    help='path to REAP executable',
                    required=False,
                    default="/humgen/atgu1/fs03/shared_resources/shared_software/bin/REAP")


# import

# argparse

# check dependencies

# run admix

# select exemplars in admix output (error out if a pop has < 10 exemplars)

# run supervised admix

# prep files for reap
# - tped
# ./plink --file mydata --recode12 --output-missing-genotype 0 --transpose --out newfile

# run reap using admix results
# use -m to reduce output
# see sec. 6 of reap documentation


# Generate diagnostic plots
# - exemplars on PCA (if PCA available)
# - final admixture on PCA (if PCA available)
# - IBD0/IBD1 for REAP

# final cleanup

# finish
