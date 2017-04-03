#! /usr/bin/env python

####################################
# gwas_logis.py
# written by Raymond Walters, March 2017
"""
Runs GWAS of plink data using logistic regression
"""
# Overview:
# 1) argparse / file checks
# 2) gwas
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

import os
import subprocess
import argparse
from warnings import warn
from args_gwas import parserbase, parsergwas
from py_helpers import unbuffer_stdout, file_len, find_exec
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='gwas_logis.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parsergwas])

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


#############
print '\n...Checking dependencies...'
# check exists, executable
#############

plinkx = find_exec('plink',key='p2loc')

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
if args.pheno is not None:
    assert os.path.isfile(args.pheno), "Phenotype file does not exist (%r)" % args.pheno

# warn if data is large
if args.extract is not None:
    nsnp = file_len(str(args.extract))
elif args.exclude is not None:
    nsnp = file_len(str(args.bfile)+'.bim') - file_len(str(args.exclude))
else:
    nsnp = file_len(str(args.bfile)+'.bim')

if nsnp > 1000000:
    warn('Large number of SNPs present for analysis (%d). Consider splitting for efficiency.' % int(nsnp))


print '\n'
print '############'
print 'Begin!'
print '############'


#############
print '\n...Running GWAS...'
#############

# setup text for covar, keep/remove, extract/exclude
if args.covar is not None:
    covar_txt = ['--covar', str(args.covar)] + cov_num_txt
    hide_cov_txt = ['hide-covar']
else:
    covar_txt = ['']
    hide_cov_txt = ['']

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
    outstem = 'logis.' + str(args.out) + '.' + str(args.addout)
else:
    outstem = 'logis.' + str(args.out)

# assemble plink call
gwas_call = [str(plinkx)] + \
                ['--bfile', str(args.bfile)] + \
		['--memory', str(2000)] + \
                keep_txt + \
                extract_txt + \
                covar_txt + \
                pheno_txt + \
                ['--logistic', 'beta'] + hide_cov_txt + \
                ['--ci', '0.95'] + \
                ['--silent'] + \
                ['--out', outstem]
gwas_call = filter(None,gwas_call)

# run
print ' '.join(gwas_call)
print ' '
subprocess.check_call(gwas_call)


#############
print '\n...Verifying output exists...'
#############

# check proper output 
assert file_len(outstem +'.assoc.logistic') >= nsnp, 'Failed to generate full GWAS results file %s' % outstem +'.assoc.logistic'


print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)
