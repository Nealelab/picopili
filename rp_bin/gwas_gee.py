#! /usr/bin/env python

####################################
# gwas_gee.py
# written by Raymond Walters, December 2015
"""
Runs GWAS of plink data using GEE and sandwich standard error
"""
# Overview:
# 1) argparse / file checks
# 2) trigger rserve
# 3) gwas
# 
# TODO: effective sample size? (see notes in gee_logit_covar.R)
# TODO: support continuous phenotypes
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
from warnings import warn
# from glob import glob
from args_gwas import *
from py_helpers import unbuffer_stdout, test_exec, find_from_path, file_len
# , read_conf, link
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='gwas_gee.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parsergwas,parsersoft])

args = parser.parse_args()

# formatting on covar-number
if args.covar_number is not None:
    cov_num_txt = ['--covar-number'] + args.covar_number
else:
    cov_num_txt = ['']

# sanity check covariate specification, keep/exclude and extract/exclude
if (args.covar_number is not None) and (args.covar is None):
    raise ValueError('Covariate number specified without a covariate file')

if (args.keep is not None) and (args.remove is not None):
    raise ValueError('Specifying both \'--keep\' and \'--remove\' is redundant. Please verify.')

if (args.extract is not None) and (args.exclude is not None):
    raise ValueError('Specifying both \'--extract\' and \'--exclude\' is redundant. Please verify.')

# get R if not provided
if args.r_ex == None or args.rscript_ex == "None":
    args.r_ex = find_from_path('R', 'R')

### print settings in use
print 'Basic Settings:'
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
print '--addout '+str(args.addout)
if args.no_cleanup:
    print '--no-cleanup '+str(args.no_cleanup)

print '\nCovariates:'
print '--covar '+str(args.covar)
if args.covar is not None:
    print ' '.join(cov_num_txt)

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
if args.rserve_active:
    print '--rserve-active '+str(args.rserve_active)
else:
    print '--r-ex '+str(args.r_ex)
print '--rplink-ex '+str(args.rplink_ex)


##############
#print '\n...Reading ricopili config file...'
##############
#
#### read plink loc from config
#conf_file = os.environ['HOME']+"/ricopili.conf"
#configs = read_conf(conf_file)


#############
print '\n...Checking dependencies...'
# check exists, executable
#############

# verify executables
test_exec(args.rplink_ex, 'Plink')
if not args.rserve_active:
    test_exec(args.r_ex, 'R')
# TODO: find a way to test Rserve available?

# check required R scripts present
rp_bin = os.path.dirname(os.path.realpath(__file__)) # use location of current script to get rp_bin
if args.covar is None:
    R_gee = rp_bin + '/gee_logit_nocov.R'
else:
    R_gee = rp_bin + '/gee_logit_covar.R'
assert os.path.isfile(R_gee), 'Failed to find R GEE script %s' % str(R_gee)

# verify bfiles are files, not paths
assert '/' not in args.bfile, "--bfile must specify only a file stem, not a path"

# verify input files exist
if args.covar is not None:
    assert os.path.isfile(args.covar), "Covariate file does not exist (%r)" % args.covar
if args.keep is not None:
    assert os.path.isfile(args.keep), "ID inclusion file does not exist (%r)" % args.keep
if args.remove is not None:
    assert os.path.isfile(args.remove), "ID exclusion file does not exist (%r)" % args.remove
if args.extract is not None:
    assert os.path.isfile(args.extract), "SNP inclusion file does not exist (%r)" % args.extract
if args.exclude is not None:
    assert os.path.isfile(args.exclude), "SNP exclusion file does not exist (%r)" % args.exclude


# warn if data is large
if args.extract is not None:
    nsnp = file_len(str(args.extract))
elif args.exclude is not None:
    nsnp = file_len(str(args.bfile)+'.bim') - file_len(str(args.exclude))
else:
    nsnp = file_len(str(args.bfile)+'.bim')

if nsnp > 50000:
    warn('Large number of SNPs present for analysis (%d). Consider splitting for efficiency.' % int(nsnp))


print '\n'
print '############'
print 'Begin!'
print '############'

#############
# Start Rserve to allow Plink R-plugins
if not args.rserve_active:
    print '\n...Starting Rserve...'
    #############

    subprocess.check_call([str(args.r_ex), 'CMD', 'Rserve', '--no-save'])
    

#############
print '\n...Running GWAS...'
#############

# setup text for covar, keep/remove, extract/exclude
if args.covar is not None:
    covar_txt = ['--covar', str(args.covar)] + cov_num_txt
else:
    covar_txt = ['']

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
    outstem = 'gee.' + str(args.out) + '.' + str(args.addout)
else:
    outstem = 'gee.' + str(args.out)

# assemble plink call
gee_call = [str(args.rplink_ex)] + \
                ['--bfile', str(args.bfile)] + \
                keep_txt + \
                extract_txt + \
                covar_txt + \
                ['--family'] + \
                ['--R', str(R_gee)] + \
                ['--out', outstem]
gee_call = filter(None,gee_call)

# run
print ' '.join(gee_call)
print ' '
subprocess.check_call(gee_call)

# check proper output 
assert os.path.isfile(outstem +'.auto.R'), 'Failed to generated GWAS results file %s' % outstem +'.auto.R'


print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)