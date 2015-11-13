####################
#
# args_ped.py
# By Raymond Walters, Nov 2015
#
# Shared argparse arguments for pedigree checking tasks. Extracted here
# to centralize arguments for the tasks and the top-level 
# driver script.
#
# Also define groups within parsers for nicer help print format.
#
# See also: args_pca.py
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

arg_base.add_argument('--input-ibd', 
                    type=str,
                    metavar='FILE',
                    help='file containing IBD relatedness estimates',
                    required=True)
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
#arg_base.add_argument('--no-cleanup',
#                    action='store_true',
#                    help='skip cleanup of interim files')



############
#
# Genotyping Rate
# Optional file with genotyping rate (for filtering)
#
############
parsergeno = argparse.ArgumentParser(add_help=False)
arg_geno = parsergeno.add_argument_group('Genotyping Rate (Optional)')

arg_geno.add_argument('--geno', 
                      type=str,
                      metavar='FILE',
                      help='file with genotype missingness rate per individual ' + \
                           '(i.e. the .imiss file from plink --missing)',
                      required=False,
                      default='NONE')  
                    

############
#
# Relatedness file handling
# Settings governing use of relatedness files
#
############
parseribd = argparse.ArgumentParser(add_help=False)
arg_ibd = parseribd.add_argument_group('Relatedness File Settings')

arg_ibd.add_argument('--format',
                      type=str.lower,
                      choices=['reap'],
                      help='format of the input IBD file',
                      required=False,
                      default='reap') 
arg_ibd.add_argument('--min-rel',
                    type=float,
                    metavar='FLOAT',
                    help='Minimum pi-hat relatedness level to treat as related. ' + \
                         'Default is halfway between 3rd and 4th degree relatives.',
                    required=False,
                    default=.09375)


############
#
# Preference Weigths
# Weights for prioritizing keep/remove decisions when filtering cryptic relatedness
#
############
parserweights = argparse.ArgumentParser(add_help=False)
arg_prefWt = parserweights.add_argument_group('Filtering Preferences (Weights)')

arg_prefWt.add_argument('--case-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of case status when selecting individuals to ' + \
                             'keep among cryptically related pairs',
                        required=False,
                        default=5.0)
arg_prefWt.add_argument('--con-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of control status when selecting individuals to ' + \
                             'keep among cryptically related pairs',
                        required=False,
                        default=2.0)
arg_prefWt.add_argument('--miss-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of missing phenotype when selecting individuals to ' + \
                             'keep among cryptically related pairs',
                        required=False,
                        default=1.0)
arg_prefWt.add_argument('--fam-case-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of a related case in pedigree when selecting ' + \
                             'individuals to keep among cryptically related pairs',
                        required=False,
                        default=5.0)
arg_prefWt.add_argument('--fam-con-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of a related control in pedigree when selecting ' + \
                             'individuals to keep among cryptically related pairs',
                        required=False,
                        default=2.0)
arg_prefWt.add_argument('--fam-miss-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of a related individual in pedigree with a ' + \
                             'missing phenotype when selecting individuals to keep ' + \
                             'among cryptically related pairs',
                        required=False,
                        default=1.0)
arg_prefWt.add_argument('--cross-fid-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of relationship to individuals in a different ' + \
                             'pedigree when selecting individuals to keep among ' + \
                             'cryptically related pairs',
                        required=False,
                        default=-10.0)
arg_prefWt.add_argument('--geno-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of genotyping rate when selecting individuals ' + \
                             'to keep among cryptically related pairs.',
                        required=False,
                        default=0.1)
arg_prefWt.add_argument('--rand-weight',
                        type=float,
                        metavar='FLOAT',
                        help='When selecting individuals to keep among cryptically ' + \
                             'related pairs, the range of the random value used to ' + \
                             'break ties. NOTE: a large value here can override the ' + \
                             'preference from the case/control/pedigree weights.',
                        required=False,
                        default=1e-5)
arg_prefWt.add_argument('--seed',
                        type=int,
                        metavar='INT',
                        help='Random seed. Applies to random values for breaking ties ' + \
                             'among cryptically related individuals',
                        required=False,
                        default=123)
                        




############
#
# Software Executables
# Locations for software dependencies not in ricopili config file
#
############
parserexloc = argparse.ArgumentParser(add_help=False)
arg_exloc = parserexloc.add_argument_group('Software Executable Locations')

arg_exloc.add_argument('--primus-ex',
                    type=str,
                    metavar='PATH',
                    help='path to main PRIMUS executable',
                    required=False,
                    default=os.environ['HOME']+"/PRIMUS_v1.8.0/bin/run_PRIMUS.pl")
arg_exloc.add_argument('--findped-ex',
                    type=str,
                    metavar='PATH',
                    help='path to PRIMUS pedigree checker executable',
                    required=False,
                    default=os.environ['HOME']+"/PRIMUS_v1.8.0/bin/find_expected_pedigree.pl")                    

