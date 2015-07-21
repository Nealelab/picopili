####################
#
# args_pca.py
# By Raymond Walters, Jul 2015
#
# Shared argparse arguments for PCA tasks. Extracted here
# to centrslize arguments for the tasks and the top-level 
# driver script.
#
# Also define groups within parsers for nicer help print format.
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
arg_base.add_argument('--no-cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')



############
#
# QC Arguments
# Arguments for the strict QC module
# 
# Includes:
# - standard SNP QC criteria
# - LD pruning crieria
# - long-range LD exclusion regions
# - Args for indels, strand-ambiguous, and non-autosomal SNPs
#
############

parserqc = argparse.ArgumentParser(add_help=False)
arg_qcth = parserqc.add_argument_group('QC Thresholds')
arg_ldparam = parserqc.add_argument_group('LD Pruning Parameters')
arg_snpcrit = parserqc.add_argument_group('Additional SNP Inclusion Criteria')

arg_qcth.add_argument('--mind-th',
                    type=float,
                    metavar='FLOAT',
                    help='individual missingness threshold',
                    required=False,
                    default=0.95)
arg_qcth.add_argument('--maf-th',
                    type=float,
                    metavar='FLOAT',
                    help='minor allele frequency threshold',
                    required=False,
                    default=0.05)
arg_qcth.add_argument('--hwe-th',
                    type=float,
                    metavar='FLOAT',
                    help='Hardy-Weinberg p-value threshold',
                    required=False,
                    default=1e-4)
arg_qcth.add_argument('--miss-th',
                    type=float,
                    metavar='FLOAT',
                    help='SNP missingness threshold',
                    required=False,
                    default=0.02)
arg_ldparam.add_argument('--ld-th',
                    type=float,
                    metavar='FLOAT',
                    help='LD pruning threshold',
                    required=False,
                    default=0.2)                    
arg_ldparam.add_argument('--ld-wind',
                    type=int,
                    metavar='INT',
                    help='LD pruning window size',
                    required=False,
                    default=200)
# arg_ldparam.add_argument('--ld-wind-move',
#                     type=int,
#                     metavar='INT',
#                     help='LD pruning window movement rate',
#                     required=False,
#                     default=100)              
arg_ldparam.add_argument('--keep-mhc',
                    action='store_true',
                    help='do not remove MHC region (chr 6: 25-35 Mb)')   
arg_ldparam.add_argument('--keep-chr8inv',
                    action='store_true',
                    help='do not remove chr. 8 inversion region (chr 8: 7-13 Mb)')
arg_ldparam.add_argument('--extra-ld-regions',
                    nargs='?',
                    metavar='FILE',
                    const='price_2008_ld_regions.txt',
                    default=None,
                    help='exclude LD regions other than than the MHC and the \
                    chr 8 inversion. If FILE specified, \
                    each line should have 3 columns indicating chromosome, \
                    starting base pair, and end base pair of region to exclude. \
                    If no file, uses regions from Price et al. (2008, AJHG).')            
arg_snpcrit.add_argument('--keep-indels',
                    action='store_true',
                    help='do not remove indels, i.e. variants with alleles I, D, -, or multiple bases')
arg_snpcrit.add_argument('--keep-strand-ambiguous',
                    action='store_true',
                    help='do not remove strand ambiguous SNPs (i.e. A/T or G/C)')
arg_snpcrit.add_argument('--all-chr',
                    action='store_true',
                    help='keep all chromosomes, instead of autosomes only')
                    
  

############
#
# PCA Arguments
# Arguments for the IMUS/PCA module
# 
# Includes:
# - relatedness threshhold for defining IMUS set
# - Number of PCs to compute
# - PCA output controls (directory, number of PCs to plot)
# - File paths for external software not provided by ~/ricopili.conf
#
############                  
                    
parserpca = argparse.ArgumentParser(add_help=False)
arg_imus = parserpca.add_argument_group('Unrelated Set')
arg_pca = parserpca.add_argument_group('Principal Components')
arg_exloc = parserpca.add_argument_group('Software Executable Locations')

arg_imus.add_argument('--rel-deg',
                    type=int,
                    metavar='INT',
                    help='relatedness degree threshold for defining \"unrelated\" set',
                    required=False,
                    default=3)
arg_pca.add_argument('--npcs',
                    type=int,
                    metavar='INT',
                    help='number of principal components to compute',
                    required=False,
                    default=10)
arg_pca.add_argument('--plot-all',
                    action='store_true',
                    help='plot all pairs of PCs, instead of top 6')
arg_pca.add_argument('--pcadir',
                    type=str,
                    metavar='DIRNAME',
                    help='name for PCA output directory, defaults to OUTNAME_imus_pca',
                    required=False)
# arg_exloc.add_argument('--plink-ex',
#                    type=str,
#                    metavar='PATH',
#                    help='path to plink executable, read from ~/ricopili.conf if unspecified',
#                    required=False)
arg_exloc.add_argument('--rscript-ex',
                    type=str,
                    metavar='PATH',
                    help='path to Rscript executable, tries reading from PATH if unspecified',
                    required=False,
                    default=None)
arg_exloc.add_argument('--primus-ex',
                    type=str,
                    metavar='PATH',
                    help='path to PRIMUS executable',
                    required=False,
                    default=os.environ['HOME']+"/PRIMUS_v1.8.0/bin/run_PRIMUS.pl")
# arg_exloc.add_argument('--smartpca-ex',
#                    type=str,
#                    metavar='PATH',
#                    help='path to smartpca executable',
#                    required=False,
#                    default="/humgen/atgu1/fs03/shared_resources/shared_software/EIG6.0beta_noreq/bin/smartpca")




############
#
# Grid Arguments
# Arguments for the cluster job submission
# 
# Includes:
# - pre-execution sleep time
# - test submission, to not actually run jobs
#
############

parsergrid = argparse.ArgumentParser(add_help=False)
arg_grid = parserpca.add_argument_group('Cluster Job Submission')

arg_grid.add_argument('--sleep',
                    type=int,
                    metavar='INT',
                    help='Pre-execution sleep time for each job',
                    required=False,
                    default=10)
arg_grid.add_argument('--test-sub',
                       action='store_true',
                       help='Test run without submitting jobs',
                       required=False)
