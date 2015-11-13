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
import subprocess
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


# print settings
print 'Using settings:'
print '--input-ibd '+str(args.input_ibd)
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
print '--format '+str(args.format)
print '--min-rel '+str(args.min_rel)

# verify input files exist
assert os.path.isfile(args.input_ibd), "IBD/relatedness file does not exist (%r)" % args.input_ibd
assert os.path.isfile(str(args.bfile)+'.fam'), "Plink fam file does not exist (%s)" % str(args.bfile)+'.fam'

# test executables
test_exec(args.primus_ex, 'PRIMUS')
test_exec(args.findped_ex, 'PRIMUS pedigree matching script')

# unzip relatedness file if needed
if args.input_ibd.endswith('.gz'):
    print 'Unzipping IBD relatedness file'
    ibd_txtfile = str(args.input_ibd) + '.txt'
    subprocess.check_call(['gunzip','-c',str(args.input_ibd),'>',str(ibd_txtfile)])
else:
    ibd_txtfile = str(args.input_ibd)

assert os.path.isfile(ibd_txtfile), "Failed to extract IBD/relatedness file (%r)" % args.input_ibd



print '\n'
print '############'
print 'Begin!'
print '############'


#############
print '\n...Constructing pedigrees from IBD/relatedness...'
# - Run using primus
# - args depend on ibd file format (default here: reap)
#############

primus_peds_log = open(str(args.out)+'.primus_peds.log', 'w')

if(args.format=='reap'):
    pr_input_text = '--input '+str(ibd_textfile)+' IBD0=6 IBD1=7 IBD2=8 RELATEDNESS=9'
else:
    raise ValueError('Unsupported IBD file format (%s)' % str(args.format))

subprocess.check_call([str(args.primus_ex),
                       str(pr_input_text),
                       '--rel_threshold',str(args.min_rel),
                       '--no_IMUS',
                       '--max_gens',str(5),
                       '--sexes',str(args.bfile)+'.fam','SEX=5',
                       '--affections',str(args.bfile)+'.fam','AFFECTION=6',
                       '--output_dir',str(args.out)+'_primus_peds'],
                       stderr=subprocess.STDOUT,
                       stdout=primus_peds_log)

primus_pr_log.close()



#############
print '\n...Preparing fam files...'
# - Format fam files for PRIMUS pedigree matching
# - Requires separate fam file per FID
# - Sibs must have parents
# TODO: find general solution for adding dummy parents
# - currently: add father/mother if >1 IID in FID and no IIDs have listed parents
#############



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
