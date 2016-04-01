#! /usr/bin/env python

####################################
# chunk_snps.py
# written by Raymond Walters, December 2015
"""
Divide SNPs into genomic chunks
"""
#
# Overview:
# 1) Read reference info on chr lengths/centromeres
# 2) Read bim file sof SNPs to be chunked
# 3) Define chunks (see logic below)
#
# logic for creating chunks:
# for each chr, for each arm (see file of hg19 chromsome lengths and centromere )
# - init region
# -- get all SNPs in region as to-do list (in curr_snps dict)
# - while SNPs remain outside existing chunks
# -- init next chunk (by larger of Mb or #SNPs)
# -- if < minSNPs remain, set end to final bp of region
# -- remove SNPs in chunk from list
# -- record start/end bp for chunk
#
# TODO:
# for now: autosomes only, ignore chr 0, assume no duplicate snps
#
####################################


import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

import os
# import subprocess
import argparse
import copy
# from glob import glob
from args_chunks import *
from py_helpers import unbuffer_stdout, file_len, warn_format
unbuffer_stdout()
import warnings
warnings.formatwarning = warn_format

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='chunk_snps.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parsersnpchunk])

args = parser.parse_args()

# derived info
bases = int(args.Mb_size * float(1e6))

if args.addout is not None:
    outname = str(args.out)+'.'+str(args.addout)+'.chunks.txt'
else:
    outname = str(args.out)+'.chunks.txt'

# retrieve default chromosome info if needed
if args.chr_info_file == 'hg19_ucsc_chrinfo.txt':
    rp_bin = os.path.dirname(os.path.realpath(__file__))
    rp_dir = os.path.dirname(rp_bin)
    chrinfo = str(rp_dir)+'/lib/hg19_ucsc_chrinfo.txt'
else:
    chrinfo = args.chr_info_file




# TODO: print argument
# TODO: verify dependencies

#############
print '\n...Checking dependencies...'
# check exists, executable
#############








print '\n'
print '############'
print 'Begin!'
print '############'

#############
print '\n...Reading chromosome length/centromere information...'
#############

chr_ref = open(chrinfo, 'r')
chrmid = {}
chrend = {}
dumphead = chr_ref.readline()
for line in chr_ref:
    (chrom, centromere, length) = line.split()
    chrmid[str(chrom)] = int(centromere)
    chrend[str(chrom)] = int(length)
chr_ref.close()

chrend_orig = copy.deepcopy(chrend)

#############
print '\n...Reading input bim file...'
#############
chroms = []
snps = {}
nbimsnps_valid = 0
bim = open(args.bfile + '.bim', 'r')
for line in bim:
    (chrom, snp_id, cm, bp, a1, a2) = line.split()
    snps[str(snp_id)] = [str(chrom), int(bp)]
    if str(chrom) not in chroms:
        chroms.append(str(chrom))
    if int(chrom) in xrange(1,23):
        nbimsnps_valid += 1
    # prevent later errors
    if int(bp) > chrend[str(chrom)]:
        warnings.warn('SNP %s (chr %s, bp %d) is outside expected chromosome bounds (bp <= %d).' % (str(snp_id), str(chrom), int(bp), int(chrend_orig[str(chrom)])))
        chrend[str(chrom)] = int(bp)
bim.close()
nbimsnps = file_len(args.bfile + '.bim')
print 'Loaded %d autosomal SNPs (of %d total in %s).' % (nbimsnps_valid, nbimsnps, bim.name)



#############
print '\n...Generating genomic chunks...'
#############
chunks = open(outname, 'w')
chunks.write(' '.join(['CHR','START','END','NAME']) + '\n')
idx = 1
nsnps = 0
for ch in xrange(1,23):

    if str(ch) not in chroms:
        continue
    
    if args.ignore_centromeres:
        arms = [0]
    else:
        arms = [0,1]
    
    for arm in arms:
    
        # init region
        if len(arms) == 2:
            if arm == 0:
                reg_start = 1
                reg_end = chrmid[str(ch)]
            else:
                reg_start = chrmid[str(ch)] + 1
                reg_end = chrend[str(ch)]
        else:
            reg_start = 0
            reg_end = chrend[str(ch)]
    
        # init chunks
        curr_snps = {k: v for k, v in snps.items() if (str(v[0]) == str(ch) and int(v[1]) >= int(reg_start) and int(v[1]) <= int(reg_end))}
        first = reg_start
        last = 0
        snps_left = len(curr_snps.keys())
        
        while snps_left > 0 and idx <= args.max_chunks:

            # debug:
            # print 'chr = %s, snps left = %s' % (str(ch),str(snps_left))

            # catch unlikely error mode        
            if snps_left < args.snp_size:
                nth_snp = sorted(curr_snps.keys(), key=lambda x: curr_snps[x][1])[snps_left-1]
                if not args.allow_small_chunks:
                    warnings.warn('Starting with too few SNPs (%d) at chr %d, bp %d. Check for sparse data or misaligned chromosome info?' % (snps_left, int(ch), first))                    
            else:
            # get bp of (minsnp)th SNP in curr_snps
                nth_snp = sorted(curr_snps.keys(), key=lambda x: curr_snps[x][1])[args.snp_size-1]            

            nth_snp_bp = curr_snps[nth_snp][1]
    
            # take later of Mb, #SNP-based location
            last = max(first+bases-1, nth_snp_bp)
            
            # remove SNPs in selected chunk
            curr_snps = {k: v for k, v in curr_snps.items() if v[1] > last}
            
            # check if enough SNPs left
            snps_left_new = len(curr_snps.keys())
            if snps_left_new < args.snp_size:
                
                # if not, then extend current chunk to end
                last = reg_end
                curr_snps = {}
                snps_left_new = 0
            
            # record
            nsnps_curr = snps_left - snps_left_new 
            start_str = str("%03d" % round(first / 1000000.0))
            end_str = str("%03d" % round(last / 1000000.0))
            ch_name = 'chr'+str(ch) +'_'+ start_str +'_'+ end_str
            chunks.write(' '.join([str(ch), str(first), str(last), str(ch_name)]) + '\n')
            
            # prep next iteration
            first = last + 1
            last = 0
            idx += 1
            nsnps += nsnps_curr
            snps_left = snps_left_new
            
chunks.close()

if idx > args.max_chunks:
    warnings.warn('Exceeded maximum number of chunks; stopping early.')
elif nsnps != nbimsnps_valid:
    raise ValueError('Number of chunked SNPs (%d) does not match number of autosomal SNPs in input .bim (%d).' % (nsnps, nbimsnps_valid))


print 'Chunks written to %s.' % outname

print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)
