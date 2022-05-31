####################
#
# args_gwas.py
# By Raymond Walters, December 2015
#
# Shared argparse arguments for GWAS. Extracted here
# to centralize arguments for the tasks and the top-level 
# driver script.
#
# Also define groups within parsers for nicer help print format.
#
# See also: args_pca.py, args_ped.py, args_qc.py
#
# TODO: deduplicate?
#
####################

# imports
import argparse
import os


############
#
# Basic Arguments
# standard I/O args shared by all modules
#
############
parserbase = argparse.ArgumentParser(add_help=False)
arg_base = parserbase.add_argument_group('Basic Arguments')

arg_base.add_argument('--bfile', 
                    type=str,
                    metavar='FILESTEM',
                    help='file stem for input plink bed/bim/fam',
                    required=True)
arg_base.add_argument('--out',
                    type=str,
                    metavar='OUTNAME',
                    help='base name for output; recommend 4 character stem to match ricopili',
                    required=True)
arg_base.add_argument('--addout',
                    type=str,
                    metavar='STR',
                    help='additional output string; intended for labelling subsets, secondary analyses, etc',
                    required=False)
arg_base.add_argument('--no-cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')


############
#
# GWAS Analysis Arguments
# - covariate settings
# - analysis subset
#
############

parsergwas = argparse.ArgumentParser(add_help=False)
arg_test = parsergwas.add_argument_group('Association Analysis')
arg_subset = parsergwas.add_argument_group('Analysis Subset')


arg_test.add_argument('--model', 
                    type=str.lower,
                    choices=['dfam','gee','gmmat','gmmat-fam','unphased','logistic','linear'],
                    help='which GWAS testing method to use for family data. Current options are plink \'--dfam\' (generalized TDT-alike), GEE (generalized estimating equations), GMMAT (logistic mixed model), UNPHASED (TDT extension by Dudbridge), or classic linear or logistic regression (i.e. without controlling for family structure)',
                    required=False,
                    default='gee')
arg_test.add_argument('--pheno', 
                    type=str,
                    metavar='FILE',
                    help='file containing alternate phenotype. Passed directly to plink. If unspecified, phenotype from fam file is used.',
                    required=False)
arg_test.add_argument('--covar', 
                    type=str,
                    metavar='FILE',
                    help='file containing analysis covariates. Passed directly to plink where applicable, or else parsed for GMMAT/UNPHASED.',
                    required=False)
arg_test.add_argument('--covar-number',
                    nargs='+',
                    metavar='COL',
                    help='which columns to use from covariate file (numbered from third column). Passed directly to plink where applicable, or else parsed for GMMAT/UNPHASED.',
                    required=False)
arg_test.add_argument('--strict-bfile',
                      type=str,
		      metavar='FILESTEM',
		      help='(GMMAT model only) file stem for input plink bed/bim/fam used to compute GRM. ' + \
		           'Defaults to using --bfile with MAF > .01, missing < .01 if unspecified.',
		      required=False)
arg_subset.add_argument('--keep',
                    type=str,
                    metavar='FILE',
                    help='file of individuals to keep for analysis. Passed directly to plink.',
                    required=False)
arg_subset.add_argument('--remove',
                    type=str,
                    metavar='FILE',
                    help='file of individuals to remove from analysis. Passed directly to plink.',
                    required=False)
arg_subset.add_argument('--extract',
                    type=str,
                    metavar='FILE',
                    help='file of SNPs to keep for analysis. Passed directly to plink.',
                    required=False)
arg_subset.add_argument('--exclude',
                    type=str,
                    metavar='FILE',
                    help='file of SNPs to remove from analysis. Passed directly to plink.',
                    required=False)


############
#
# UNPHASED
#
############

arg_unphased = parsergwas.add_argument_group('UNPHASED Model Settings')

arg_unphased.add_argument('--unphased-model',
                    type=str.lower,
                    choices=['commonmain','full','haplomain','allelemain','gxg','null'],
                    help='which model to use within UNPHASED. See UNPHASED documentation for details. Options other than \'commonmain\' are currently untested, and may break output formatting currently assumed by picopili here.',
                    required=False,
                    default='commonmain')
arg_unphased.add_argument('--unphased-no-missing',
                    action='store_true',
		    help='run UNPHASED without \'-missing\', i.e. will exclude individuals with missing genotypes instead of averaging over possible genotypes in family. Dropping \'-missing\' should be faster, at the cost of reduced robustness to missingness.')
arg_unphased.add_argument('--unphased-linkage',
                    action='store_true',
		    help='run UNPHASED without \'-nolinkage\', i.e. do not assume siblings in a family have independent transmissions from one another. Assuming no linkage will be faster and more powerful, but may be incorrect for regions with linkage (e.g. in presence of association if siblings were selected based on case status)')
arg_unphased.add_argument('--unphased-no-parentrisk',
                    action='store_true',
		    help='run UNPHASED without \'-parentrisk\', i.e. do not model the genetic risk of parents. Dropping \'-parentrisk\' should run faster and increase power, but at the cost of less valid inference in the presence of missing parents or multiple siblings in a family.')


############
#
# Parallelization settings
#
############

parserchunk = argparse.ArgumentParser(add_help=False)
arg_snpchunk = parserchunk.add_argument_group('Parallel Jobs')

arg_snpchunk.add_argument('--snp-chunk', 
                    type=int,
                    metavar='INT',
                    help='Number of SNPs to analyze in each parallel chunk',
                    required=False,
                    default=9950)



############
#
# Aggregation settings
#
############

parseragg = argparse.ArgumentParser(add_help=False)
arg_agg = parseragg.add_argument_group('GWAS Results Filtering')

arg_agg.add_argument('--maf-a-th', 
                    type=float,
                    metavar='FLOAT',
                    help='Threshold for minor allele frequency in cases to include in GWAS results',
                    required=False,
                    default=.005)
arg_agg.add_argument('--maf-u-th', 
                    type=float,
                    metavar='FLOAT',
                    help='Threshold for minor allele frequency in controls to include in GWAS results',
                    required=False,
                    default=.005)
arg_agg.add_argument('--info-file', 
                    type=str,
                    metavar='FILE',
                    help='File containing info scores for imputed SNPs',
                    required=False)
arg_agg.add_argument('--info-th', 
                    type=float,
                    metavar='FLOAT',
                    help='Threshold for imputation info score to include in GWAS results',
                    required=False,
                    default=.6)
arg_agg.add_argument('--p-th2', 
                    type=float,
                    metavar='FLOAT',
                    help='p-value threshold for inclusion in the p-sorted top results output',
                    required=False,
                    default=1e-3)
arg_agg.add_argument('--max-se', 
                    type=float,
                    metavar='FLOAT',
                    help='Maximum SE allowed for GWAS results. Only applies to GEE model. Useful for filtering out numerically unstable results.',
                    required=False,
                    default=100.0)

############
#
# Software settings
#
############

parsersoft = argparse.ArgumentParser(add_help=False)
arg_soft = parsersoft.add_argument_group('Software')
arg_clust = parsersoft.add_argument_group('Cluster Settings')
arg_exloc = parsersoft.add_argument_group('Executable Locations')

#arg_soft.add_argument('--rserve-active',
#                    action='store_true',
#                    help='skip launching Rserve. Without this argument, will try \'R CMD Rserve\' to enable Plink-R plugin interface.')
arg_clust.add_argument('--sleep', 
                    type=int,
                    metavar='SEC',
                    help='Number of seconds to delay on start of cluster jobs',
                    required=False,
                    default=30)
arg_exloc.add_argument('--r-ex',
                    type=str,
                    metavar='PATH',
                    help='path to R executable, tries reading from PATH if unspecified',
                    required=False,
                    default=None)
arg_exloc.add_argument('--rscript-ex',
                    type=str,
                    metavar='PATH',
                    help='path to Rscript executable, tries reading from PATH if unspecified',
                    required=False,
                    default=None)
arg_exloc.add_argument('--rserve-ex',
                    type=str,
		    metavar='PATH',
		    help='path to Rserve executable. Defaults to picopili.conf value. Only required for \'--model gee\'.',
		    required=False,
		    default=None)
arg_exloc.add_argument('--rplink-ex',
                    type=str,
                    metavar='PATH',
                    help='path to plink executable with R plugin interface and \'--dfam\' (e.g. plink 1.7 or current versions of 1.9, but not older 1.9 betas). Default behavior relies on picopili.conf, but this arg is kept for legacy support.',
                    required=False,
                    default=None)
arg_soft.add_argument('--plink-mem',
                    type=int,
		    metavar='MB',
		    help='Memory to use for top level plink calls (freq and chunking). Default ' + \
		    'should be fine for most datasets, may need to increase for large datasets.',
		    default=2000,
		    required=False)
arg_clust.add_argument('--port',
                    type=int,
                    metavar='TCP_PORT',
                    help='TCP port to use for Rserve (for \'--model gee\' only)',
                    required=False,
                    default=6311)

# eof
