#! /usr/bin/env python

####################################
# ped_confirm.py
# written by Raymond Walters, November 2015
"""
Confirm reported pedigrees consistent with relatedness-based reconstruction 
"""
# Overview:
# 1) Run PRIMUS to estimate pedigrees from relatedness
# 2) Compare estimated to reported pedigrees
#       - split fam file by FID
#       - add dummy parents if needed
#       - run find_expected_pedigree.pl from PRIMUS
# 3) Parse logs to find possible mismatches
# 4) Ouput results
#
####################################


####################################
# Setup
# a) load python dependencies
# b) get variables/arguments
####################################

import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
# import os
# import subprocess
import argparse
# from string import ascii_uppercase
# from glob import glob
# from numpy import digitize
# import random
# import warnings
from args_ped import *
from py_helpers import unbuffer_stdout
# file_len, test_exec, read_conf, find_from_path, link, gz_confirm
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='ped_confirm.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parseribd,parserexloc])

args = parser.parse_args()





# PRIMUS
# --input FILE=(reap, zip not allowed) IBD0=6 IBD1=7 IBD2=8 RELATEDNESS=9
# --rel_threshold .09375
# --output_dir
# --no_IMUS
# --sexes FILE=.fam SEX=5 
# --max_gens 5
# --affections FILE=fam AFFECTION=6 

# awk '{print $1,$2,"fid"$1"father","fid"$1"mother",$5,$6 > "pitt_clean_pruned_famex_fid"$1".fam"}' pitt_clean_pruned_famex.fam

# for fi in `wc -l pitt_clean_pruned_famex_fid*.fam | tr -s ' ' | sed 's/^ //' | awk '{if($1>1) print $2}'`; do echo "#### $fi ####" >> pitt_confirm_peds.log; ~/PRIMUS_v1.8.0/bin/find_expected_pedigree.pl ../pitt_clean_pruned_famex_PRIMUS/Summary_pitt_clean_pruned_famex_cleaned.genome.txt ${fi} | tail -n 2 >> pitt_confirm_peds.log; echo "" >> pitt_confirm_peds.log; echo "" >> pitt_confirm_peds.log; done

### sample output

#dataset_dir = ./pitt3.REAP_pairs_relatedness.txt_PRIMUS/
#Ref_fam file = pitt_clean_pruned_famex_fid1358.fam
#Family network of interest: 2
#
#
#intersection: 135803 135804
#fam file: ./pitt3.REAP_pairs_relatedness.txt_PRIMUS//pitt3.REAP_pairs_relatedness.txt_network2/pitt3.REAP_pairs_relatedness.txt_network2_1.fam
#MATCH!!! network2 fam 1: ./pitt3.REAP_pairs_relatedness.txt_PRIMUS//pitt3.REAP_pairs_relatedness.txt_network2/pitt3.REAP_pairs_relatedness.txt_network2_1.fam

###### or 

#dataset_dir = ./pitt3.REAP_pairs_relatedness.txt_PRIMUS/
#Ref_fam file = pitt_clean_pruned_famex_fid1358.fam
#Family network of interest: 2
#
#
#intersection: 135803 135804
#fam file: ./pitt3.REAP_pairs_relatedness.txt_PRIMUS//pitt3.REAP_pairs_relatedness.txt_network2/pitt3.REAP_pairs_relatedness.txt_network2_1.fam
#NO MATCHING PEDIGREE FOUND


# grep -B2 "NO MATCH" pitt_confirm_peds.log
