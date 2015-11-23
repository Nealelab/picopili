####################
#
# args_qc.py
# By Raymond Walters, Nov 2015
#
# Shared argparse arguments for QC. Extracted here
# to centralize arguments for the tasks and the top-level 
# driver script.
#
# Also define groups within parsers for nicer help print format.
#
# See also: args_pca.py, args_ped.py
#
# TODO: deduplicate args_pca to use args from here
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
arg_base.add_argument('--no-cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')



############
#
# Pre-Proccessing Arguments
# Args for data labelling inherited from ricopili
# 
# Includes:
# - disease name
# - population name
# - platform checking
#
############

parsertag = argparse.ArgumentParser(add_help=False)
arg_tag = parsertag.add_argument_group('ID Tagging Information')

arg_tag.add_argument('--skip-fid-tags',
                     action='store_true',
                     help='stops ricopili-style tags from being added to individual FIDs')
arg_tag.add_argument('--disname',
                     type=str,
                     metavar='DIS',
                     help='disease name to use in FID, filename tags. For ' + \
                          'consistency with ricopili, a 3-character code is recommended.',
                     required=False,
                     default='dis')
arg_tag.add_argument('--popname',
                     type=str,
                     metavar='POP',
                     help='population name to use in FID, filename tags. ' + \
                          'Recommended options: mix, eur, asn, afr, aam.',
                     required=False,
                     default='mix')
arg_tag.add_argument('--skip-platform',
                     action='store_true',
                     help='skip ricopili platform guessing. If skipped, will denote ' + \
                          'platform as NONE in tags. Automatically skipped if ' + \
                          'FID tags are skipped.')



############
#
# QC Arguments
# Arguments for the general QC module
# 
# Includes:
# - ID QC criteria
# - SNP QC criteria
#
############

parserqc = argparse.ArgumentParser(add_help=False)
arg_indqc = parserqc.add_argument_group('Individual QC Criteria')
arg_snpqc = parserqc.add_argument_group('SNP QC Criteria')

arg_indqc.add_argument('--mind-th',
                       type=float,
                       metavar='FLOAT',
                       help='individual missingness threshold',
                       required=False,
                       default=0.95)
arg_indqc.add_argument('--het-th',
                       type=float,
                       metavar='FLOAT',
                       help='Absolute F_het threshold for heterozygosity',
                       required=False,
                       default=0.2)
arg_indqc.add_argument('--skip-sex-check',
                       action='store_true',
                       help='skip sex check with chrX heterozygosity'), 
arg_snpqc.add_argument('--pre-miss',
                       type=float,
                       metavar='FLOAT',
                       help='Threshold for pre-filter on SNP missingness (before QC of IDs)',
                       required=False,
                       default=0.05)
arg_snpqc.add_argument('--miss-th',
                       type=float,
                       metavar='FLOAT',
                       help='SNP missingness threshold',
                       required=False,
                       default=0.02)
arg_snpqc.add_argument('--diff-miss',
                       type=float,
                       metavar='FLOAT',
                       help='Absolute SNP differential missingness (case vs. control) threshold',
                       required=False,
                       default=0.02)
arg_snpqc.add_argument('--hwe-th-cas',
                       type=float,
                       metavar='FLOAT',
                       help='Hardy-Weinberg p-value threshold for HWE among cases',
                       required=False,
                       default=1e-10)
arg_snpqc.add_argument('--hwe-th-con',
                       type=float,
                       metavar='FLOAT',
                       help='Hardy-Weinberg p-value threshold for HWE among controls',
                       required=False,
                       default=1e-6)
arg_snpqc.add_argument('--hwe-th-all',
                       type=float,
                       metavar='FLOAT',
                       help='Hardy-Weinberg p-value threshold for HWE in all founders, including missing phenotypes',
                       required=False,
                       default=1e-6)
arg_snpqc.add_argument('--maf-th',
                       type=float,
                       metavar='FLOAT',
                       help='minor allele frequency threshold. Default is no filter',
                       required=False,
                       default=-1.0)



############
#
# Mendal Error Arguments
# Arguments for checking mendelian errors in family data
# 
# Includes:
# - which mendel error test to use
# - QC thresholds for IDs, SNPs
# - whether to zero out remaining mendal errors
#
############

parsermendel = argparse.ArgumentParser(add_help=False)
arg_mendel = parsermendel.add_argument_group('Mendelian Error Checks')

arg_mendel.add_argument('--mendel',
                        type=str.lower,
                        choices=['none', 'trios', 'duos', 'multigen'],
                        help='Mendel error testing method. All methods come from plink2, ' + \
                        'and differ according to behavior when parental genotypes are missing. ' + \
                        'See plink2 documentation on "--mendel" (trios), "--mendel-duos", and "--mendel-multigen".',  
                        default='none',
                        required=False)
arg_mendel.add_argument('--id-mendel-th',
                        type=float,
                        metavar='FLOAT',
                        help='Threshold for proportion of SNPs with mendel errors, otherwise individual excluded',
                        default=.02,
                        required=False)
arg_mendel.add_argument('--snp-mendel-th',
                        type=float,
                        metavar='FLOAT',
                        help='Threshold for proportion of individuals with mendel errors, otherwise SNP excluded',
                        default=.01,
                        required=False)
arg_mendel.add_argument('--keep-mendel',
                        action='store_true',
                        help='Prevents setting mendelian errors to missing')


# eof