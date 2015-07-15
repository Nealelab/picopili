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
# 3) Get SNPs to exlcude
#     - not ACGT (e.g. indels)
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
parser.add_argument('--keep-indels',
                    action='store_true',
                    help='do not remove indels, i.e. variants with alleles I, D, -, or multiple bases')                    
parser.add_argument('--extra-ld-regions',
                    action='store_true',
                    help='exclude additional LD regions from Price et al. (2008, AJHG)')
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
parser.add_argument('--all-chr',
                    action='store_true',
                    help='keep all chromosomes, instead of autosomes only')
parser.add_argument('--no-cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')

args, pass_through_args = parser.parse_known_args()

# derived arguments
ld_move = int(args.ld_wind / 2)

# print settings
print 'Using settings:'
print '--bfile '+args.bfile
print '--out '+args.out
print '--mind-th '+str(args.mind_th)
print '--mind-th '+str(args.mind_th)
print '--npcs '+str(args.npcs)


 
#############
print '\n...reading ricopili config file...'
#############



print '...reading ricopili config file...'
### read plink loc from config

conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

plinkx = configs['p2loc']+"plink"



print '...getting descriptive with plink...'
### get descriptives, exclude high mind
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



print '...finding indels, strand ambiguous SNPs, and long LD regions...'
### get strand ambi list, mhc/etc liost
bim_in_nam = args.bfile + '.bim'
# ambiex_nam = args.out + '_ambiexclude.txt'
# ldex_nam = args.out + '_ldexclude.txt'
snpout_nam = args.out + '.exclude_snps.txt'

snp_in = open(bim_in_nam, 'r')
# ambiex_out = open(ldex_nam, 'w')
# ldex_out = open(ldex_nam, 'w')
snp_out = open(snpout_nam, 'w')

indels = ['i','d','-']

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
   
    if (a1l=='a') and (a2l=='t'):
#        ambiex_out.write(snp + '\n')
        snp_out.write(snp + ' strand_ambiguous\n')
    elif (a1l=='t') and (a2l=='a'):
#        ambiex_out.write(snp + '\n')
        snp_out.write(snp + ' strand_ambiguous\n')
    elif (a1l=='g') and (a2l=='c'):
#        ambiex_out.write(snp + '\n')
        snp_out.write(snp + ' strand_ambiguous\n')
    elif (a1l=='c') and (a2l=='g'):
#        ambiex_out.write(snp + '\n')
        snp_out.write(snp + ' strand_ambiguous\n')
        
    if (chrom==6) and (bp > 25000000) and (bp < 35000000):
#        ldex_out.write(snp + '\n')
        snp_out.write(snp + ' mhc_region\n')
    elif (chrom==8) and (bp > 7000000) and (bp < 13000000):
#        ldex_out.write(snp + '\n')
        snp_out.write(snp + ' chr8inv_region\n')
    elif args.extra_ld_regions:
        if (chrom==1) and (bp > 48000000) and (bp < 52000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==2) and (bp > 86000000) and (bp < 101000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==2) and (bp > 134000000) and (bp < 138000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==2) and (bp > 183000000) and (bp < 190000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==3) and (bp > 47000000) and (bp < 50000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==3) and (bp > 83000000) and (bp < 87000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==3) and (bp > 89000000) and (bp < 98000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==5) and (bp > 44000000) and (bp < 51000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==5) and (bp > 98000000) and (bp < 101000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==5) and (bp > 129000000) and (bp < 132000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==5) and (bp > 135000000) and (bp < 139000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==6) and (bp > 57000000) and (bp < 64000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==6) and (bp > 140000000) and (bp < 143000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==7) and (bp > 55000000) and (bp < 66000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==8) and (bp > 43000000) and (bp < 50000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==8) and (bp > 112000000) and (bp < 115000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==10) and (bp > 37000000) and (bp < 43000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==11) and (bp > 46000000) and (bp < 57000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==11) and (bp > 87000000) and (bp < 91000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==12) and (bp > 33000000) and (bp < 40000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==12) and (bp > 109000000) and (bp < 112000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')
        elif (chrom==20) and (bp > 32000000) and (bp < 35000000):
#            ldex_out.write(snp + '\n')
            snp_out.write(snp + ' longLD_region\n')

snp_in.close()



print '...finding low MAF SNPs...'
### get low maf
frq_nam = sumstat_out + '.frq'

frqs = open(frq_nam, 'r')
dumphead = frqs.readline()

for line in frqs:
    (chrom, snp, a1, a2, maf, nobs) = line.split()

    if float(maf) < args.maf_th or (1.0-float(maf)) < args.maf_th:
        snp_out.write(snp + ' low_MAF\n')

frqs.close()



print '...finding HWE failures...'
### get hwe failure
hwe_nam = sumstat_out + '.hwe'

hwes = open(hwe_nam, 'r')
dumphead = hwes.readline()

for line in hwes:
    (chrom, snp, test, a1, a2, geno, Ohet, Ehet, p) = line.split()
    
    if (test == "ALL") and ( float(p) < args.hwe_th ):
        snp_out.write(snp + ' HWE_fail\n')

hwes.close()



print '...Finding call rate failures...'
### get lmissing 
lmiss_nam = sumstat_out + '.lmiss'

lmiss = open(lmiss_nam, 'r')
dumphead = lmiss.readline()

for line in lmiss:
    (chrom, snp, nmiss, ngeno, fmiss) = line.split()
    
    if float(fmiss) > args.miss_th:
        snp_out.write(snp + ' high_missing\n')

lmiss.close()
snp_out.close()



print '...Removing filtered SNPs...'
### run plink to exclude failures
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




print '...beginning LD pruning...'
### ld prune (loop, apply)

# init
i = 1

subprocess.check_call([str(plinkx), 
               "--bfile", filtered_out,
               "--indep-pairwise", str(args.ld_wind), str(ld_move), str(args.ld_th),
               "--silent",
               "--allow-no-sex",
               "--out", args.out + '.prune' + str(i) + '.tmp' ])

nprune_old = file_len(filtered_out + '.bim')
nprune_new = file_len(args.out + '.prune' + str(i) + '.tmp.prune.in')

# loop til no additional exclusions
while nprune_old > nprune_new:
    i += 1
    print '...LD pruning pass ' + str(i) + '...'
    subprocess.check_call([str(plinkx), 
               "--bfile", filtered_out,
               "--extract", args.out + '.prune' + str(i-1) + '.tmp.prune.in',
               "--indep-pairwise", str(args.ld_wind), str(ld_move), str(args.ld_th),
               "--silent",
               "--allow-no-sex",
               "--out", args.out + '.prune' + str(i) + '.tmp' ])

    nprune_old = nprune_new
    nprune_new = file_len(args.out + '.prune' + str(i) + '.tmp.prune.in')  

print '...extracting LD pruned set...'
# apply
subprocess.check_call([str(plinkx), 
               "--bfile", filtered_out,
               "--extract", args.out + '.prune' + str(i) + '.tmp.prune.in',
               "--make-bed",
               "--silent",
               "--allow-no-sex",
               "--out", args.out + '.strictqc.pruned' ])



# cleanup
if not args.no_cleanup:
    print '...cleaning up files...'
    print 'zipping to ' + args.out + '.qc_files.tar.gz:'
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
    
    subprocess.check_call(["gzip", "-f", snpout_nam])

    print 'remove interim:'
    subprocess.check_call(["rm",
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
    
    subprocess.check_call(["rm"] + glob(args.out+".prune*.tmp.*"))
    
    # allowing failure, since files may or may not exists
    print 'remove if exist:' 
    subprocess.call(["rm", 
                     sumstat_out + '.hh',
                     sumstat_out + '.nosex',
                     filtered_out + '.hh',
                     filtered_out + '.nosex',
                     args.out + '.strictqc.pruned.hh',
                     args.out + '.strictqc.pruned.nosex'])

print '############'
print '\n'
print 'SUCCESS!\n'
exit(0)
