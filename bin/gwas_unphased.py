#! /usr/bin/env python

####################################
# gwas_unphased.py
# written by Raymond Walters, Nov 2021
"""
Runs GWAS of plink data using UNPHASED
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
import pandas as pd
from warnings import warn
from args_gwas import parserbase, parsergwas, parsersoft
from py_helpers import unbuffer_stdout, file_len, find_exec
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='gwas_unphased.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parsergwas,parsersoft])
parser.add_argument('--phenofile',
			type=str,
			help='formatted fam/ped file, expected to be created in gwas_rel.py',
			required=True)
parser.add_argument('--covfile',
                        type=str,
			help='formatted covariate file, expected to be created in gwas_rel.py',
			required=False)

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
print '--phenofile ' + str(args.phenofile)

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

unphasedx = find_exec('unphased.lnx',key='unphasedloc')

# verify bfiles are files, not paths
assert '/' not in args.bfile, "--bfile must specify only a file stem, not a path"

# verify input files exist
if args.covar is not None:
    assert os.path.isfile(args.covfile), "Covariate file does not exist (%r)" % args.covfile

assert os.path.isfile(args.phenofile), "Phenotype file does not exist (%r)" % args.phenofile


# warn if data is large
if args.extract is not None:
    nsnp = file_len(str(args.extract))
elif args.exclude is not None:
    nsnp = file_len(str(args.bfile)+'.bim') - file_len(str(args.exclude))
else:
    nsnp = file_len(str(args.bfile)+'.bim')

if nsnp > 10000:
    warn('Large number of SNPs present for analysis (%d). Consider splitting for efficiency.' % int(nsnp))


print '\n'
print '############'
print 'Begin!'
print '############'

#############
print '\n...Getting covariates...'
#############

cov_in = pd.read_csv(str(args.covfile), header=0, delim_whitespace=True, dtype=str, nrows=2)
ncovs = len(cov_in.columns)
covnames = [str(cov_in.columns.values.tolist()[x]) for x in range(5,ncovs)]

#############
print '\n...Running GWAS...'
#############

# setup text for covar
if args.covar is not None:
    if args.covar_number is not None:
        args.covar_number = ','.join(args.covar_number)
        # print(args.covar_number)
        covnums = []
	cov_sets = str(args.covar_number).split(",")
	for cset in cov_sets:
	    crange = cset.split("-")
	    if len(crange) == 1:
	        covnums = covnums + crange
	    elif len(crange) == 2:
	        covnums = covnums + [x for x in range(int(crange[0]), int(crange[1])+1)]
	    else:
	        raise argparse.ArgumentError(args.covar_number,"Unexpected part of --covar-number: "+str(cset))

	filter(None,covnums)
        covar_txt = ['-phenofile'] + [args.covfile] + ['-confounder'] + [str(covnames[int(x)-1]) for x in covnums]

    else:
        covar_txt = ['-phenofile'] + [args.covfile] + ['-confounder'] + covnames


else:
    covar_txt = ['']

# setup output name
if args.addout is not None:
    outstem = 'unphased.' + str(args.out) + '.' + str(args.addout)
else:
    outstem = 'unphased.' + str(args.out)

# setup unphased modelling args
unphased_args = ['-individual', '-time']
if args.unphased_model is not None:
	unphased_args = unphased_args + ['-model', str(args.unphased_model)]
else:
        unphased_args = unphased_args + ['-model', 'commonmain']

if not args.unphased_no_missing:
        unphased_args = unphased_args + ['-missing']

if not args.unphased_linkage:
        unphased_args = unphased_args + ['-nolinkage']

if not args.unphased_no_parentrisk:
        unphased_args = unphased_args + ['-parentrisk']

# assemble UNPHASED call
gwas_call = [str(unphasedx)] + \
                ['-bedfile', str(args.bfile)+'.bed'] + \
		['-mapfile', str(args.bfile)+'.bim'] + \
		['-pedfile', str(args.phenofile)] + \
		covar_txt + \
		unphased_args + \
		['-zero', '0'] + \
		['-output', str(outstem)+'.log'] + \
		['-tabularfile', str(outstem)] # + \
		# ['-marker', 'rs1229984:100239319:T:C'] # 'rs10435479']

# tabularfile = stem of <filename>.Disease.{family,unrelated}.out

gwas_call = filter(None,gwas_call)

# run
print ' '.join(gwas_call)
print ' '
# print gwas_call
subprocess.check_call(gwas_call)


#############
print '\n...Verifying output exists...'
#############

# check proper output 
assert file_len(outstem +'.Disease.family.out') >= (3*nsnp + 1), 'Failed to generate full GWAS results file %s' % outstem +'.Disease.family.out'


print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)
