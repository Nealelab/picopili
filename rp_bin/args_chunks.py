####################
#
# args_chunks.py
# By Raymond Walters, December 2015
#
# Shared argparse arguments for chunking data for parallelization. Extracted
# here to centralize arguments for the tasks and the top-level 
# driver scripts.
#
# Also define groups within parsers for nicer help print format.
#
# See also: args_pca.py, args_ped.py, args_qc.py, etc
#
# TODO: deduplicate?
#
####################

# imports
import argparse
# import os



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
#arg_base.add_argument('--no-cleanup',
#                    action='store_true',
#                    help='skip cleanup of interim files')



############
#
# SNP Chunking Parameters
# settings to specify SNP chunking priorities
#
############

parsersnpchunk = argparse.ArgumentParser(add_help=False)
arg_snpchunk = parsersnpchunk.add_argument_group('SNP Chunking Parameters')

arg_snpchunk.add_argument('--Mb-size', 
                    type=float,
                    metavar='FLOAT',
                    help='Minimum size of chunk, in Mb (megabases)',
                    required=True,
                    default=3.0)
arg_snpchunk.add_argument('--snp-size', 
                    type=int,
                    metavar='INT',
                    help='Minimum size of chunk, in number of SNPs',
                    required=True,
                    default=50)
arg_snpchunk.add_argument('--max-chunks', 
                    type=int,
                    metavar='INT',
                    help='Maximum number of chunks to create. Will throw ' + \
                         'a warning if maximum is reached before all SNPs are chunked.',
                    required=True,
                    default=5000)
arg_snpchunk.add_argument('--ignore-centromeres', 
                    action='store_true',
                    help='Allow chunks to span centromere boundaries')
arg_snpchunk.add_argument('--chr-info-file', 
                    type=str,
                    metavar='FILE',
                    help='file with chromosome length and centromere locations. ' + \
                         'Default file is retrieved from installation location.',
                    required=False,
                    default='hg19_ucsc_chrinfo.txt')
arg_snpchunk.add_argument('--allow-small-chunks', 
                    action='store_true',
                    help='Allow chunks with fewer than \'--snp-size\' markers without a warning. ' + \
                         'Such chunks may occur due to sparse data (e.g. few SNPs ' + \
                         'on the short arm of chr21) or could indicate bad chromosome build information.')

# eof