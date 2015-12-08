####################
#
# args_impute.py
# By Raymond Walters, December 2015
#
# Shared argparse arguments for imputation. Extracted
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
                    metavar='STRING',
                    help='additional output string; intended for labelling subsets, secondary analyses, etc',
                    required=False)
#arg_base.add_argument('--no-cleanup',
#                    action='store_true',
#                    help='skip cleanup of interim files')




############
#
# Prephasing Parameters
# settings for prephasing with shapeit
#
############

parserphase = argparse.ArgumentParser(add_help=False)
arg_align = parserphase.add_argument_group('Reference Alignment Settings')
arg_shape = parserphase.add_argument_group('SHAPEIT Arguments')
arg_refloc = parserphase.add_argument_group('Reference File Locations')
arg_submit = parserphase.add_argument_group('Cluster Submission Settings')

arg_align.add_argument('--popname', 
                       type=str.lower,
#                       choices=['eur','afr','asn','eas','sas','amr','asw','mix'],
                       metavar='STR',
                        help='Ancestral population (e.g. \'afr\'; used for getting reference allele freqs)',
                        required=False,
                        default='eur')
arg_align.add_argument('--sfh', 
                       type=float,
                       metavar='FLOAT',
                       help='secure frequency margin for flipping strand ambiguous SNPs',
                       required=False,
                       default=0.2)
arg_align.add_argument('--fth', 
                       type=float,
                       metavar='FLOAT',
                       help='allowed frequency difference compared to the reference',
                       required=False,
                       default=0.15)
arg_align.add_argument('--refdir', 
                       type=str,
                       metavar='PATH',
                        help='Ricopili reference directory. Used to get allele frequencies during alignment',
                        required=False,
                        default='/psych/genetics_data/ripke/1000Genomes_reference/1KG_Oct14/1000GP_Phase3_sr/subchr/')
arg_shape.add_argument('--window',
                        type=int,
                        metavar='INT',
                        help='window size for shapeit, in megabases (Mb)',
                        required=False,
                        default=5)
# TODO: fix this
arg_shape.add_argument('--refstem',
                        type=int,
                        metavar='INT',
                        help='reference to use with shapeit. CURRENTLY UNUSED, hardcoded to shared 1KG Phase 3',
                        required=False)
arg_shape.add_argument('--seed',
                        type=int,
                        metavar='INT',
                        help='random seed for shapeit',
                        required=False,
                        default=12345)
arg_refloc.add_argument('--map-dir', 
                        type=str,
                        metavar='PATH',
                        help='Directory with genomic maps, per chromosome. ' + \
                             'Expected filenames are ./genetic_map_chr${i}_combined_b37.txt',
                        required=False,
                        default='/humgen/atgu1/fs03/shared_resources/1kG/shapeit/genetic_map')
arg_submit.add_argument('--sleep',
                        type=int,
                        metavar='INT',
                        help='wait time before executing UGER tasks, in seconds',
                        required=False,
                        default=30)
arg_submit.add_argument('--mem-req',
                        type=int,
                        metavar='INT',
                        help='memory to request for shapeit, per chromosome, in GB',
                        required=False,
                        default=4)
arg_submit.add_argument('--threads',
                        type=int,
                        metavar='INT',
                        help='parallel threads to use for shapeit, per chromosome',
                        required=False,
                        default=4)



############
#
# Imputation Parameters
#
############

parserimpute = argparse.ArgumentParser(add_help=False)
arg_imp = parserimpute.add_argument_group('IMPUTE2 Arguments')


# eof



