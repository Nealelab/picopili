#! /usr/bin/env python

####################################
# gwas_dfam.py
# written by Raymond Walters, December 2015
"""
Runs GWAS of plink data using plink DFAM
"""
# Overview:
# 1) argparse / file checks
# 2) gwas
# 
# TODO: effective sample size?
#
####################################


# ish benchmarks: 
#   50 snp, n=473, 1 covar -> 2s
#   50 snp, n=473, no covar -> 2s
#   5k snp, n=473, no covar -> 2m, 30s
# 110k snp, n=473, no covar -> ~1 hr (cf. near-instant w/ dfam)


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

import os
import subprocess
import argparse
# from glob import glob
from args_gwas import *
from py_helpers import unbuffer_stdout, test_exec, find_exec
# , read_conf, link
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='gwas_dfam.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parsergwas,parsersoft])

args = parser.parse_args()

# sanity check covariate specification, keep/exclude and extract/exclude
if (args.covar is not None) or (args.covar is not None):
    raise ValueError('Covariates not allowed for Plink DFAM analysis')

if (args.keep is not None) and (args.remove is not None):
    raise ValueError('Specifying both \'--keep\' and \'--remove\' is redundant. Please verify.')

if (args.extract is not None) and (args.exclude is not None):
    raise ValueError('Specifying both \'--extract\' and \'--exclude\' is redundant. Please verify.')


### print settings in use
print 'Basic Settings:'
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
print '--addout '+str(args.addout)
if args.no_cleanup:
    print '--no-cleanup '+str(args.no_cleanup)
    
if args.pheno is not None and str(args.pheno) != "None":
    print '--pheno '+str(args.pheno)

print '\nAnalysis Subset:'
if args.keep is not None:
    print '--keep '+str(args.keep)
else:
    print '--remove '+str(args.remove)
if args.extract is not None:
    print '--extract '+str(args.extract)
else:
    print '--exclude '+str(args.exclude)

print '\nSoftware Settings:'
print '--rplink-ex '+str(args.rplink_ex)



#############
print '\n...Checking dependencies...'
# check exists, executable
#############

# R-compatible plink
if args.rplink_ex is None or args.rplink_ex == "None":
    args.rplink_ex = find_exec('plink',key='rplloc')

test_exec(args.rplink_ex, 'Plink')

# verify bfiles are files, not paths
assert '/' not in args.bfile, "--bfile must specify only a file stem, not a path"

# verify input files exist
if args.keep is not None:
    assert os.path.isfile(args.keep), "ID inclusion file does not exist (%r)" % args.keep
if args.remove is not None:
    assert os.path.isfile(args.remove), "ID exclusion file does not exist (%r)" % args.remove
if args.extract is not None:
    assert os.path.isfile(args.extract), "SNP inclusion file does not exist (%r)" % args.extract
if args.exclude is not None:
    assert os.path.isfile(args.exclude), "SNP exclusion file does not exist (%r)" % args.exclude
if args.pheno is not None:
    assert os.path.isfile(args.pheno), "Phenotype file does not exist (%r)" % args.pheno


print '\n'
print '############'
print 'Begin!'
print '############'

#############
print '\n...Running GWAS...'
#############

# setup text for keep/remove, extract/exclude
if args.pheno is not None and str(args.pheno) != "None":
    pheno_txt = ['--pheno', str(args.pheno)]
else:
    pheno_txt = ['']

if args.keep is not None:
    keep_txt = ['--keep', str(args.keep)]
elif args.remove is not None:
    keep_txt = ['--remove', str(args.remove)]
else:
    keep_txt = ['']

if args.extract is not None:
    extract_txt = ['--extract', str(args.extract)]
elif args.exclude is not None:
    extract_txt = ['--exclude', str(args.exclude)]
else:
    extract_txt = ['']

# setup output name
if args.addout is not None:
    outstem = 'dfam.' + str(args.out) + '.' + str(args.addout)
else:
    outstem = 'dfam.' + str(args.out)

# assemble plink call
dfam_call = [str(args.rplink_ex)] + \
                ['--bfile', str(args.bfile)] + \
                keep_txt + \
                extract_txt + \
                pheno_txt + \
                ['--dfam', '--prune'] + \
                ['--silent', '--memory', str(2000)] + \
                ['--out', outstem]
dfam_call = filter(None,dfam_call)

# run
print ' '.join(dfam_call)
print ' '
subprocess.check_call(dfam_call)

# check proper output 
assert os.path.isfile(outstem +'.dfam'), 'Failed to generated GWAS results file %s' % outstem +'.dfam'


print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)