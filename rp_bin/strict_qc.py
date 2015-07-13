#! /usr/bin/env python

# load requirements
import os
import subprocess
import argparse


# parse arguments
parser = argparse.ArgumentParser(prog='strict_qc.py',
                                 formatter_class=lambda prog:
                                 argparse.HelpFormatter(prog, max_help_position=40))

parser.add_argument('--input', '-i', 
                    type=str,
                    metavar='FILESTEM',
                    help='file stem for input plink bed/bim/fam',
                    required=True)
parser.add_argument('--output', '-o', 
                    type=str,
                    metavar='NAME',
                    help='base name for output; recommend 4 character stem to match ricopili',
                    required=True)
parser.add_argument('--mind-th',
                    type=float,
                    metavar='FLOAT',
                    help='individual missingness threshold',
                    required=False,
                    default=0.95)
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
                    help='keep all chrosomes (default: only autosomes)')
parser.add_argument('--no-cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')

args, pass_through_args = parser.parse_known_args()

ld_move = int(args.ld_wind / 2)

print ld_move



### read plink loc from config

conf_file = os.environ['HOME']+"/ricopili.conf"

configs = {}
with open(conf_file, 'r') as f:
    for line in f:
        (key, val) = line.split()
        configs[str(key)] = val

plinkx = configs['p2loc']+"plink"

print plinkx



### get descriptives, exclude high mind
sumstat_out = args.output+".qcsumstat"

subprocess.check_call([str(plinkx), 
               "--bfile", args.input,
               "--mind", str(args.mind_th),
               "--freq",
               "--missing",
               "--hardy",
               "--out", sumstat_out])


### get strand ambi list, mhc/etc liost
bim_in_nam = args.input + '.bim'
# ambiex_nam = args.output + '_ambiexclude.txt'
# ldex_nam = args.output + '_ldexclude.txt'
snpout_nam = args.output + '_exclude_snps.txt'

snp_in = open(bim_in_nam, 'r')
# ambiex_out = open(ldex_nam, 'w')
# ldex_out = open(ldex_nam, 'w')
snp_out = open(snpout_nam, 'w')

for line in snp_in:
    (chrom, snp, cm, bp, a1, a2) = line.split()
    
    if (a1.lower()=='a') and (a2.lower()=='t'):
#        ambiex_out.write(snp + '\n')
        snp_out.write(snp + ' strand_ambiguous\n')
    elif (a1.lower()=='t') and (a2.lower()=='a'):
#        ambiex_out.write(snp + '\n')
        snp_out.write(snp + ' strand_ambiguous\n')
    elif (a1.lower()=='g') and (a2.lower()=='c'):
#        ambiex_out.write(snp + '\n')
        snp_out.write(snp + ' strand_ambiguous\n')
    elif (a1.lower()=='c') and (a2.lower()=='g'):
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



### get low maf
frq_nam = sumstat_out + '.frq'

frqs = open(frq_nam, 'r')
dumphead = frqs.readline()

for line in frqs:
    (chrom, snp, a1, a2, maf, nobs) = line.split()

    if maf < args.maf_th or maf > args.maf_th:
        snp_out.write(snp + ' lowMAF\n')

frqs.close()



# get hwe failure
hwe_nam = sumstat_out + '.hwe'

hwes = open(hwe_nam, 'r')
dumphead = hwes.readline()

for line in frqs:
    (chrom, snp, test, a1, a2, geno, Ohet, Ehet, p) = line.split()
    
    if (test == "ALL") and ( p < args.hwe_th ):
        snp_out.write(snp + ' HWE_fail\n')

hwes.close()



# get lmissing 


# combine exclude lists


# run plink to exclude


# ld prune (loop, apply)


# cleanup



# exit
