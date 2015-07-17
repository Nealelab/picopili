#! /usr/bin/env python

####################################
# strict_qc.py
# written by Raymond Walters, July 2015
"""
Runs strict QC for GWAS data
"""
# Overview:
# 1) Input QCed plink bed/bim/fam
# 2) Get QC metrics with plink
# 3) Get SNPs to exclude
#     - indels
#     - strand ambiguous
#     - long-range LD regions
#     - low MAF
#     - HWE failures
#     - low call rate
# 4) Remove failed SNPs
# 5) LD prune
# 6) Clean up files
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
from glob import glob
from py_helpers import file_len, read_conf



#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='strict_qc.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))

parser.add_argument('--bfile', 
                    type=str,
                    metavar='FILESTEM',
                    help='file stem for input plink bed/bim/fam',
                    required=True)
parser.add_argument('--out',
                    type=str,
                    metavar='OUTNAME',
                    help='base name for output; recommend 4 character stem to match ricopili',
                    required=True)
parser.add_argument('--mind-th',
                    type=float,
                    metavar='FLOAT',
                    help='individual missingness threshold',
                    required=False,
                    default=0.95)
parser.add_argument('--maf-th',
                    type=float,
                    metavar='FLOAT',
                    help='minor allele frequency threshold',
                    required=False,
                    default=0.05)
parser.add_argument('--hwe-th',
                    type=float,
                    metavar='FLOAT',
                    help='Hardy-Weinberg p-value threshold',
                    required=False,
                    default=1e-4)
parser.add_argument('--miss-th',
                    type=float,
                    metavar='FLOAT',
                    help='SNP missingness threshold',
                    required=False,
                    default=0.02)
parser.add_argument('--ld-th',
                    type=float,
                    metavar='FLOAT',
                    help='LD pruning threshold',
                    required=False,
                    default=0.2)                    
parser.add_argument('--ld-wind',
                    type=int,
                    metavar='INT',
                    help='LD pruning window size',
                    required=False,
                    default=200)
# parser.add_argument('--ld-wind-move',
#                     type=int,
#                     metavar='INT',
#                     help='LD pruning window movement rate',
#                     required=False,
#                     default=100)
parser.add_argument('--keep-indels',
                    action='store_true',
                    help='do not remove indels, i.e. variants with alleles I, D, -, or multiple bases')
parser.add_argument('--keep-strand-ambiguous',
                    action='store_true',
                    help='do not remove strand ambiguous SNPs (i.e. A/T or G/C)')
parser.add_argument('--keep-mhc',
                    action='store_true',
                    help='do not remove MHC region (chr 6: 25-35 Mb)')   
parser.add_argument('--keep-chr8inv',
                    action='store_true',
                    help='do not remove chr. 8 inversion region (chr 8: 7-13 Mb)')
parser.add_argument('--extra-ld-regions',
                    nargs='?',
                    metavar='FILE',
                    const='price_2008_ld_regions.txt',
                    default=None,
                    help='exclude LD regions other than than the MHC and the \
                    chr 8 inversion. If FILE specified, \
                    each line should have 3 columns indicating chromosome, \
                    starting base pair, and end base pair of region to exclude. \
                    If no file, uses regions from Price et al. (2008, AJHG).')
parser.add_argument('--all-chr',
                    action='store_true',
                    help='keep all chromosomes, instead of autosomes only')
parser.add_argument('--no-cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')

args = parser.parse_args()

# derived arguments
ld_move = int(args.ld_wind / 2)

# print settings
print 'Using settings:'
print '--bfile '+args.bfile
print '--out '+args.out
print '--mind-th '+str(args.mind_th)
print '--keep-indels '+str(args.keep_indels)
print '--keep-strand-ambiguous '+str(args.keep_strand_ambiguous)
print '--keep-mhc '+str(args.keep_mhc)
print '--keep-chr8inv '+str(args.keep_chr8inv)
print '--keep-indels '+str(args.keep_indels)
print '--extra-ld-regions '+str(args.extra_ld_regions)
print '--maf-th '+str(args.maf_th)
print '--hwe-th '+str(args.hwe_th)
print '--miss-th '+str(args.miss_th)
print '--ld-th '+str(args.ld_th)
print '--ld_wind '+str(args.ld_wind)
print '--all_chr '+str(args.all_chr)


 
#############
print '\n...Reading ricopili config file...'
#############

### read plink loc from config

conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

plinkx = configs['p2loc']+"plink"

# get directory containing current script
# (hack to help find ld region text file)
rp_bin = os.path.dirname(os.path.realpath(__file__))


#############
print '\n...Checking dependencies...'
# check exists, executable
#############

# plink
assert os.path.isfile(plinkx), "Plink not found at %r" % plinkx
assert os.access(plinkx, os.X_OK), "Plink not executable (%r)" % plinkx
print "Plink found: %s" % plinkx

# ld region file, if needed
# try in rp_bin/lib/ in addition to cwd
if args.extra_ld_regions != None:
    if os.path.isfile(args.extra_ld_regions):
        print "LD region file found: %s" %  args.extra_ld_regions
    elif os.path.isfile(str(rp_bin + '/lib/' + args.extra_ld_regions)):
        args.extra_ld_regions = str(rp_bin + '/lib/' + args.extra_ld_regions)
        print "LD region file found: %s" %  args.extra_ld_regions
    else:
        raise IOError("LD region file %s not found in current directory or %s." % (args.extra_ld_regions, str(rp_bin + '/lib/')))



print '\n'
print '############'
print 'Begin!'
print '############'



####################################
# Get QC metrics with plink
# - exclude mostly missing IDs (>95% missing)
# - allele freqs
# - missingness rate
# - hardy-weinberg
####################################

#############
print '\n...Getting QC metrics...'
#############

sumstat_out = args.out+".qcsumstat"

subprocess.check_call([str(plinkx), 
               "--bfile", args.bfile,
               "--mind", str(args.mind_th),
               "--freq",
               "--missing",
               "--hardy",
               "--silent",
               "--allow-no-sex",
               "--out", sumstat_out])



####################################
# Identify exclusions in bim file
# - indels (based on alleles of I, D, -, or multiple bases)
# - strand ambiguous
# - long LD regions
# write all to file with reason (snpout_nam)
####################################

#############
print '\n...Finding indels, strand ambiguous SNPs, and long LD regions...'
#############

### define indel alleles
indels = ['i','d','-']

### prep extra long LD regions
if args.extra_ld_regions != None:
    ld_reg = []
    with open(args.extra_ld_regions, 'r') as f:
        for line in f:
            ld_reg.append((line.split()))


### prep bim input, snp exclusion output
bim_in_nam = args.bfile + '.bim'
snpout_nam = args.out + '.exclude_snps.txt'

snp_in = open(bim_in_nam, 'r')
snp_out = open(snpout_nam, 'w')

for line in snp_in:
    (chrom, snp, cm, bp, a1, a2) = line.split()
    chrom = int(chrom)
    bp = int(bp)
    
    a1l = a1.lower()
    a2l = a2.lower()
    
    if not args.keep_indels:
        if len(a1l) > 1 or len(a2l) > 1:
            snp_out.write(snp + ' indel_allele\n')    
        elif any(a == a1l for a in indels) or any(a == a2l for a in indels):
            snp_out.write(snp + ' indel_allele\n')
   
   if no args.keep_strand_ambiguous:
        if (a1l=='a') and (a2l=='t'):
            snp_out.write(snp + ' strand_ambiguous\n')
        elif (a1l=='t') and (a2l=='a'):
            snp_out.write(snp + ' strand_ambiguous\n')
        elif (a1l=='g') and (a2l=='c'):
            snp_out.write(snp + ' strand_ambiguous\n')
        elif (a1l=='c') and (a2l=='g'):
            snp_out.write(snp + ' strand_ambiguous\n')

    if not args.keep_mhc:    
        if (chrom==6) and (bp > 25000000) and (bp < 35000000):
            snp_out.write(snp + ' mhc_region\n')
            
    if not args.keep_chr8inv:
        elif (chrom==8) and (bp > 7000000) and (bp < 13000000):
            snp_out.write(snp + ' chr8inv_region\n')
        
    elif args.extra_ld_regions != None:
        if any([chrom==int(ld_chr) and int(ld_bp_start) < bp < int(ld_bp_end) for (ld_chr, ld_bp_start, ld_bp_end) in ld_reg]):
            snp_out.write(snp + ' longLD_region\n')
            

snp_in.close()



####################################
# Identify low MAF SNPs
# - write to snp_out_nam
####################################

#############
print '\n...Finding low MAF SNPs...'
#############

frq_nam = sumstat_out + '.frq'
frqs = open(frq_nam, 'r')
dumphead = frqs.readline()

for line in frqs:
    (chrom, snp, a1, a2, maf, nobs) = line.split()

    if float(maf) < args.maf_th or (1.0-float(maf)) < args.maf_th:
        snp_out.write(snp + ' low_MAF\n')

frqs.close()


####################################
# Identify HWE failures
# - write to snp_out_nam
####################################

#############
print '\n...Finding HWE failures...'
#############

hwe_nam = sumstat_out + '.hwe'
hwes = open(hwe_nam, 'r')
dumphead = hwes.readline()

for line in hwes:
    (chrom, snp, test, a1, a2, geno, Ohet, Ehet, p) = line.split()
    
    if (test == "ALL") and ( float(p) < args.hwe_th ):
        snp_out.write(snp + ' HWE_fail\n')

hwes.close()


####################################
# Identify call rate failures
# - write to snp_out_nam
####################################

#############
print '\n...Finding call rate failures...'
#############

lmiss_nam = sumstat_out + '.lmiss'
lmiss = open(lmiss_nam, 'r')
dumphead = lmiss.readline()

for line in lmiss:
    (chrom, snp, nmiss, ngeno, fmiss) = line.split()
    
    if float(fmiss) > args.miss_th:
        snp_out.write(snp + ' high_missing\n')

lmiss.close()
snp_out.close()



####################################
# Remove QC failures using plink
# - autosomes only, unless otherwise specified
# - also removes to high missingness IDs
####################################

#############
print '\n...Removing filtered SNPs...'
#############

filtered_out = args.out+".strictqc"

if args.all_chr:
    subprocess.check_call([str(plinkx), 
                   "--bfile", args.bfile,
                   "--mind", str(args.mind_th),
                   "--exclude", snpout_nam,
                   "--make-bed",
                   "--silent",
                   "--allow-no-sex",
                   "--out", filtered_out])
else:
   subprocess.check_call([str(plinkx), 
               "--bfile", args.bfile,
               "--mind", str(args.mind_th),
               "--exclude", snpout_nam,
               "--autosome",
               "--make-bed",
               "--silent",
               "--allow-no-sex",
               "--out", filtered_out]) 



####################################
# LD prune the QC+ SNPs
# - using in-sample LD
# - for robustness, loop pruning until no further exclusions
# - once stable, extract LD pruned set
####################################

#############
print '\n...Beginning LD pruning...'
#############

# init
i = 1

subprocess.check_call([str(plinkx), 
               "--bfile", filtered_out,
               "--indep-pairwise", str(args.ld_wind), str(ld_move), str(args.ld_th),
               "--silent",
               "--allow-no-sex",
               "--out", args.out + '.prune' + str(i) + '.tmp' ])

# tracking number of SNPs before, after altest round of pruning
nprune_old = file_len(filtered_out + '.bim')
nprune_new = file_len(args.out + '.prune' + str(i) + '.tmp.prune.in')

# loop til no additional exclusions
while nprune_old > nprune_new:

    i += 1
    #############
    print 'Pruning pass ' + str(i)
    #############
    subprocess.check_call([str(plinkx), 
               "--bfile", filtered_out,
               "--extract", args.out + '.prune' + str(i-1) + '.tmp.prune.in',
               "--indep-pairwise", str(args.ld_wind), str(ld_move), str(args.ld_th),
               "--silent",
               "--allow-no-sex",
               "--out", args.out + '.prune' + str(i) + '.tmp' ])

    nprune_old = nprune_new
    nprune_new = file_len(args.out + '.prune' + str(i) + '.tmp.prune.in')  


#############
print '\n...Extracting LD pruned set...'
#############

subprocess.check_call([str(plinkx), 
               "--bfile", filtered_out,
               "--extract", args.out + '.prune' + str(i) + '.tmp.prune.in',
               "--make-bed",
               "--silent",
               "--allow-no-sex",
               "--out", args.out + '.strictqc.pruned' ])



####################################
# Clean up files
# - tar.gz the plink QC metrics, logs, and final LD pruning
# - zip SNP exclusion list with reasons
# - remove interim bed/bim/fam (QC+ pre-pruning), unused metrics, tar-ed files
# - remove interim pruning results
# - remove extraneuous .hh and .nosex files if present 
####################################

if not args.no_cleanup:
    
    #############
    print '\n...Cleaning up files...'
    #############
    
    #############
    print 'Zipping to ' + args.out + '.qc_files.tar.gz:'
    #############
    subprocess.check_call(["tar", "-zcvf",
                           args.out + '.qc_files.tar.gz',
                           sumstat_out + '.log',
                           frq_nam,
                           hwe_nam, 
                           lmiss_nam,
                           filtered_out + '.log',
                           args.out + '.prune' + str(i) + '.tmp.prune.in',
                           args.out + '.prune' + str(i) + '.tmp.log',
                           ])
    #############
    print 'Compress:'
    #############
    subprocess.check_call(["gzip", "-fv", snpout_nam])

    #############
    print 'Remove interim:'
    #############
    subprocess.check_call(["rm", "-v",
                           filtered_out + '.bed',
                           filtered_out + '.bim',
                           filtered_out + '.fam',
                           sumstat_out + '.imiss',
                           frq_nam,
                           hwe_nam,
                           lmiss_nam,
                           sumstat_out + '.log',
                           filtered_out + '.log',
                           ])
    
    subprocess.check_call(["rm", "-v"] + glob(args.out+".prune*.tmp.*"))
    
    #############
    print 'Remove if exist:'
    #############
    # allowing failure, since files may or may not exists
    subprocess.call(["rm", "-v",
                     sumstat_out + '.hh',
                     sumstat_out + '.nosex',
                     filtered_out + '.hh',
                     filtered_out + '.nosex',
                     args.out + '.strictqc.pruned.hh',
                     args.out + '.strictqc.pruned.nosex'])



print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)
