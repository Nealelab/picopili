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
arg_submit = parserphase.add_argument_group('SHAPEIT Resource Requirements')

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
arg_align.add_argument('--ref-info', 
                       type=str,
                       metavar='PATH',
                        help='gzipped file of reference information, with columns ' + \
                                '"id","position","a0","a1", and $popname, where $popname' + \
                                'contains the allele frequency for the "a1" allele. Can ' + \
                                'include "###" in place of chromosome number (as in default).',
                        required=False,
                        default='/humgen/atgu1/fs03/shared_resources/1kG/shapeit/1000GP_Phase3_chr###.legend.gz')
arg_shape.add_argument('--window',
                        type=float,
                        metavar='FLOAT',
                        help='window size for shapeit, in megabases (Mb)',
                        required=False,
                        default=5.0)
arg_shape.add_argument('--no-duohmm',
                        action='store_true',
                        help='omit --duohmm flag for family-aware pre-phasing in shapeit',
                        required=False)
arg_shape.add_argument('--shape-seed',
                        type=int,
                        metavar='INT',
                        help='random seed for shapeit',
                        required=False,
                        default=12345)
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
# SNP Chunking Parameters
# 
# - Subset of of options from args_chunks.py
# - Remaining arguments (e.g. centromeres, short chunks) fixed for definining imputation chunks
#
############

parserchunk = argparse.ArgumentParser(add_help=False)
arg_snpchunk = parserchunk.add_argument_group('SNP Chunking')

arg_snpchunk.add_argument('--Mb-size', 
                    type=float,
                    metavar='FLOAT',
                    help='Minimum size of chunk, in Mb (megabases)',
                    required=False,
                    default=3.0)
arg_snpchunk.add_argument('--snp-size', 
                    type=int,
                    metavar='INT',
                    help='Minimum size of chunk, in number of SNPs',
                    required=False,
                    default=30)
arg_snpchunk.add_argument('--chr-info-file', 
                    type=str,
                    metavar='FILE',
                    help='file with chromosome length and centromere locations. ' + \
                         'Default file is retrieved from installation location.',
                    required=False,
                    default='hg19_ucsc_chrinfo.txt')

############
#
# Imputation Parameters
#
############

parserimpute = argparse.ArgumentParser(add_help=False)
arg_imp = parserimpute.add_argument_group('IMPUTE2 Arguments')

arg_imp.add_argument('--Ne',
                     type=int,
                     metavar='INT',
                     help='effective population size for imputation',
                     required=False,
                     default=20000)
arg_imp.add_argument('--buffer',
                     type=int,
                     metavar='KB',
                     help='size of buffer region, in kb, to use around target region when imputing genomic chunks',
                     required=False,
                     default=1000)
arg_imp.add_argument('--imp-seed',
                     type=int,
                     metavar='INT',
                     help='random seed for impute2',
                     required=False,
                     default=54321)


############
#
# Imputation Reference
#
############

parserref = argparse.ArgumentParser(add_help=False)
arg_ref = parserref.add_argument_group('Imputation Reference')

arg_ref.add_argument('--ref-maps', 
                     type=str,
                     metavar='FILENAME',
                     help='Genomic maps. To specify files split by chromosome, use "###" to indicate chromosome number (see default).',
                     required=False,
                     default='/humgen/atgu1/fs03/shared_resources/1kG/shapeit/genetic_map/genetic_map_chr###_combined_b37.txt')
arg_ref.add_argument('--ref-haps',
                     type=str,
                     metavar='FILENAME',
                     help='Imputation reference .hap.gz file for shapeit and impute2. Can use "###" to indicate chromosome number (see default).',
                     required=False,
                     default='/humgen/atgu1/fs03/shared_resources/1kG/shapeit/1000GP_Phase3_chr###.hap.gz')
arg_ref.add_argument('--ref-legs',
                     type=str,
                     metavar='FILENAME',
                     help='Imputation reference .legend.gz file for shapeit and impute2. Can use "###" to indicate chromosome number (see default).',
                     required=False,
                     default='/humgen/atgu1/fs03/shared_resources/1kG/shapeit/1000GP_Phase3_chr###.legend.gz')
arg_ref.add_argument('--ref-samps',
                     type=str,
                     metavar='FILENAME',
                     help='Imputation reference .sample file for shapeit and impute2. Can use "###" to indicate chromosome number.',
                     required=False,
                     default='/humgen/atgu1/fs03/shared_resources/1kG/shapeit/1000GP_Phase3.sample')

############
#
# Best guess filtering
#
############

parserbg = argparse.ArgumentParser(add_help=False)
arg_bg = parserbg.add_argument_group('Best-Guess Genotypes')

arg_bg.add_argument('--bg-th',
                        type=float,
                        metavar='FLOAT',
                        help="Minimum posterior probability for making best-guess calls",
                        required=False,
                        default=0.8)
arg_bg.add_argument('--info-th',
                        type=float,
                        metavar='FLOAT',
                        help="Info score threshold for filtering best-guess genotypes",
                        required=False,
                        default=0.6)
arg_bg.add_argument('--keep-mendel',
                        action='store_true',
                        help='Prevents setting mendelian errors to missing')
arg_bg.add_argument('--maf-th',
                        type=float,
                        metavar='FLOAT',
                        help="Minor allele frequency threshold for filtering best-guess genotypes",
                        required=False,
                        default=0.005)
arg_bg.add_argument('--miss-th',
                        type=float,
                        metavar='FLOAT',
                        help="SNP missingness threshold for filtering best-guess genotypes",
                        required=False,
                        default=0.02)
arg_bg.add_argument('--mac-th',
                        type=int,
                        metavar='INT',
                        help="Minor allele count threshold for filtering best-guess genotypes",
                        required=False,
                        default=None)
arg_bg.add_argument('--hard-call-th',
                        type=float,
                        metavar='FLOAT',
                        help="Maximum uncertainty for making best-guess calls. Passed directly to plink --hard-call-threshold.",
                        required=False,
                        default=None)
arg_bg.add_argument('--max-info-th',
                        type=float,
                        metavar='FLOAT',
                        help="Maximum info score threshold for filtering best-guess genotypes",
                        required=False,
                        default=2.0)
arg_bg.add_argument('--mendel',
                        type=str.lower,
                        choices=['none', 'trios', 'duos', 'multigen'],
                        help='Mendel error testing method. All methods come from plink2, ' + \
                        'and differ according to behavior when parental genotypes are missing. ' + \
                        'See plink2 documentation on "--mendel" (trios), "--mendel-duos", and "--mendel-multigen".',  
                        default='multigen',
                        required=False)


############
#
# Cluster Settings
#
############

parsercluster = argparse.ArgumentParser(add_help=False)
arg_clust = parsercluster.add_argument_group('Cluster Settings')

arg_clust.add_argument('--sleep', 
                    type=int,
                    metavar='SEC',
                    help='Number of seconds to delay on start of cluster jobs',
                    required=False,
                    default=30)
arg_clust.add_argument('--full-pipe', 
                    action='store_true',
                    help='Proceed through full imputation pipeline',
                    required=False)

# eof
