#! /usr/bin/env python

####################################
# qc_rel.py
# written by Raymond Walters, November 2015
"""
Runs QC for GWAS data
"""
# Overview:
# 1) Input QCed plink bed/bim/fam
# 2) Guess genotyping platform
#     - platform guessing script from ricopili (plague)
#     - parse output here
# 3) Tag IDs in bim file
#     - tag script from ricopili (id_tager - need to modify to refine output structure)
# 4) Pre-filter SNP call rate
# 5) QC individuals
#     - missing
#     - heterozygosity
#     - mendel errors
#     - sex check
# 6) QC SNPs
#     - missings
#     - differential call rate
#     - mendel errors
#     - hwe (case, control)
# 7) Resolve remaining mendel errors
# 8) Generate basic summary info
# 9) Clean up files
#
# TODO: Unify this with strict_qc.py
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
import warnings
# from glob import glob
from args_qc import *
from py_helpers import unbuffer_stdout, read_conf, test_exec, link, file_len
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='qc_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parserqc,parsermendel,parsertag])

args = parser.parse_args()



### print settings in use
print 'Basic Settings:'
print '--bfile '+args.bfile
print '--out '+args.out

print '\nID Tagging Information:'
print '--skip-fid-tags '+str(args.skip_fid_tags)
if not args.skip_fid_tags:
    print '--disname '+str(args.disname)
    print '--popname '+str(args.popname)
    print '--skip-platform '+str(args.skip_platform)

print '\nIndividual QC Thresholds:'
print '--mind-th '+str(args.mind_th)
print '--het-th '+str(args.het_th)
print '--skip-sex-check '+str(args.skip_sex_check)

print '\nSNP QC Thresholds:'
print '--pre-miss '+str(args.pre_miss)
print '--miss-th '+str(args.miss_th)
print '--diff_miss '+str(args.diff_miss)
print '--hwe-th-cas '+str(args.hwe_th_cas)
print '--hwe-th-con '+str(args.hwe_th_cas)
print '--hwe-th-all '+str(args.hwe_th_cas)
if args.maf_th >= 0.0:
    print '--maf-th '+str(args.maf_th)

print '\nMendelian Error Checks:'
print '--mendel '+str(args.mendel)
if args.mendel is not 'none':
    print '--id-mendel-th '+str(args.id_mendel_th)
    print '--snp-mendel-th '+str(args.snp_mendel_th)
    print '--keep-mendel '+str(args.keep_mendel)
print ' '


#############
print '\n...Reading ricopili config file...'
#############

### read plink loc from config
conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

plinkx = configs['p2loc']+"plink"

analyst = configs['init']

if not args.skip_platform:
    # get directory containing current script
    # (hack to get plague script location)
    rp_bin = os.path.dirname(os.path.realpath(__file__))
    plague_ex = rp_bin + '/plague.pl'


#############
print '\n...Checking dependencies...'
# check exists, executable
#############

# verify executables
test_exec(plinkx, 'Plink')
if not args.skip_platform:
    test_exec(plague_ex, 'Platform guessing script')

# verify bfiles are files, not paths
assert '/' not in args.bfile, "--bfile must specify only a file stem, not a path"


print '\n'
print '############'
print 'Begin!'
print '############'

#############
qcdir = 'qc_'+str(args.out)
print '\n...Setting up working directory (./%s)...' % qcdir
#############

wd = os.getcwd()

if os.path.isdir(qcdir):
    raise IOError ('Output directory %s already exists. Stopping to prevent overwriting files.' % qcdir)
else:
    os.makedirs(qcdir)

os.chdir(qcdir)

# link plink files (with verification)
link(str(wd+'/'+args.bfile+'.bed'), str(args.bfile+'.bed'), 'input bed file')
link(str(wd+'/'+args.bfile+'.bim'), str(args.bfile+'.bim'), 'input bim file')
if args.skip_fid_tags:
    link(str(wd+'/'+args.bfile+'.fam'), str(args.bfile+'.fam'), 'input fam file')
else:
    link(str(wd+'/'+args.bfile+'.fam'), str(args.bfile+'.fam.original'), 'input fam file')

# setup SNP, ID QC fail lists
snpout_nam = args.out + '.exclude_snps.txt'
idout_nam = args.out + '.exclude_iids.txt'


#############
if not args.skip_fid_tags:
    print '\n...Preparing tags for FIDs...'
#############

    if args.skip_platform:
        plat = "NONE"
    
    else:
        # run plague script from ricopili
        plague_file = open(str(args.out)+'.plague.out', 'w')
        subprocess.check_call([plague_ex,
                              str(args.bfile+'.bim')],
                              stdout = plague_file)
        
        plague_file.close()

        # initialize
        plat = '?'
        sum_percent_max = 0.0
        
        # parse plague results
        # using same criteria as ricopili (from preimp_dir)
        # - add plague percentages for platform, take highest
        plague = open(str(args.out)+'.plague.out', 'r')
        
        with plague as f:
            # dump first 3 lines
            for _ in xrange(3):
                dump = f.next()
            # loop plague results
            for line in f:
                fields = line.split()
                print 'try ' + str(fields[2])
                sum_percent = float(fields[7]) + float(fields[14])
                print 'sum percent: ' + str(sum_percent)
                if sum_percent > sum_percent_max:
                    plat_full = str(fields[2]).split('_')
                    plat = plat_full[-1]
                    sum_percent_max = sum_percent
        
        plague.close()


    fidtag = 'fam_'+str(args.disname)+'_'+str(args.out)+'_'+str(args.popname)+'_'+str(analyst)+'_'+plat
    print 'Using tag: %s' % fidtag
    
    # tag IDs (easier to do here than modify id_tager_2 output)
    origfam = open(str(args.bfile+'.fam.original'), 'r')
    famout = open(str(args.bfile+'.fam'), 'w')
    for line in origfam:
        outline = str(fidtag)+'*'+line
        famout.write(outline)
        
    origfam.close()
    famout.close()


#############
print '\n...Pre-filtering SNPs on call rate...'
#############

# get SNP call rates
# plus allele freqs for individual QC
prefilter_stats = str(args.out) + ".prefilter_stats"

prefilter_stats_str = [str(plinkx), 
                               "--bfile", args.bfile,
                               "--missing",
                               "--make-founders","require-2-missing",
                               "--freqx",
                               "--silent",
                               "--allow-no-sex",
                               "--out", prefilter_stats]

print 'Running: ' + ' '.join(prefilter_stats_str)
subprocess.check_call(prefilter_stats_str)

# get failing SNPs from plink lmiss file
snp_out = open(snpout_nam, 'w')

prefilter_lmiss_nam = prefilter_stats + '.lmiss'
prefilter_lmiss = open(prefilter_lmiss_nam, 'r')
nex = 0

dumphead = prefilter_lmiss.readline()
for line in prefilter_lmiss:
    (chrom, snp, nmiss, ngeno, fmiss) = line.split()
    
    if float(fmiss) > args.pre_miss:
        snp_out.write(snp + ' prefilter_high_missing\n')
        nex += 1

prefilter_lmiss.close()
snp_out.close()
print 'Excluding %d SNPs failing pre-filter' % nex

# remove SNPs

prefilter_out = args.out+".prefiltered"

prefilter_out_str = [str(plinkx), 
                               "--bfile", args.bfile,
                               "--exclude", snpout_nam,
                               "--silent",
                               "--allow-no-sex",
                               "--make-bed",
                               "--out", prefilter_out]

print ''
print 'Running: ' + ' '.join(prefilter_out_str)
subprocess.check_call(prefilter_out_str)



#############
print '\n...QCing individuals...'
# - get QC metrics
# - find missing rate failures
# - find F_het heterozygosity failures
# - find mendel failures
# -- as proportion of total snps
# - find sex check failures
# - remove all failing individuals
#############

ind_stats = str(args.out) + '.ind_stats'

# setup sex check, if desired
if args.skip_sex_check:
    sexcheck_txt = ''
else:
    # TODO: split PAR region    

    # check number of chr X snps
    bim_prefilter = open(str(prefilter_out)+'.bim', 'r')
    n_chrx = 0
    
    for line in bim_prefilter:
        (chrom, snp, cm, bp, a1, a2) = line.split()
        if chrom == 23 or str(chrom) == 'X':
            n_chrx += 1
    
    bim_prefilter.close()
    
    if n_chrx < args.min_chrx_snps:
        warnings.warn('Insufficent chrX SNPs for sex check (%d). Sex check filter will be omitted.' % n_chrx)
        args.skip_sex_check = True
        sexcheck_txt = ''
    else:
        sexcheck_txt = '--check-sex'


# setup mendel error check, if desired
if args.mendel == 'none':
    mendel_txt = ''
    
elif args.mendel == 'trios':
    mendel_txt = '--mendel'
    
elif args.mendel == 'duos':
    mendel_txt = '--mendel-duos'
    
elif args.mendel == 'multigen':
    mendel_txt = '--mendel-multigen'
    
else:
    # shouldn't be possible from argparse
    raise ValueError ("Invalid argument for '--mendel' (%s)" % str(args.mendel))


# get plink sum stats
ind_stats_str = [str(plinkx), 
                         "--bfile", prefilter_out,
                         "--missing",
                         "--het",
                         sexcheck_txt,
                         mendel_txt,
                         "--read-freq",str(prefilter_stats)+'.frqx',
                         "--silent",
                         "--allow-no-sex",
                         "--out", ind_stats]

# remove empty elements
ind_stats_str = filter(None, ind_stats_str)

print 'Running: ' + ' '.join(ind_stats_str)
subprocess.check_call(ind_stats_str)


# get number of SNPs as denominator for mendel error rate
nsnp = file_len(str(prefilter_out) + '.fam')

# ouput exclusions
id_out = open(idout_nam, 'w')


# filter missingness
imiss_nam = ind_stats + '.imiss'
imiss = open(imiss_nam, 'r')
nex = 0

dumphead = imiss.readline()
for line in imiss:
    (fid, iid, miss_pheno, n_miss, n_geno, f_miss) = line.split()
    
    if float(f_miss) > args.mind_th:
        id_out.write(str(fid) + ' ' + str(iid) + ' high_missing\n')
        nex += 1

imiss.close()
print 'Found %d individuals to exclude for failing missingness rate' % nex


# filter heterozygosity
het_nam = ind_stats + '.het'
het = open(het_nam, 'r')
nex = 0

dumphead = het.readline()
for line in het:
    (fid, iid, o_hom, e_hom, ngeno, fhet) = line.split()
    
    if float(fhet) > args.het_th:
        id_out.write(str(fid) + ' ' + str(iid) + ' high_homozygosity_Fhet\n')
        nex += 1
    elif float(fhet) < (-1.0*args.het_th):
        id_out.write(str(fid) + ' ' + str(iid) + ' low_homozygosity_Fhet\n')
        nex += 1

het.close()
print 'Found %d individuals to exclude for failing Fhet homozygosity rate' % nex


# filter mendel errors
if args.mendel != 'none':

    imendel_nam = ind_stats + '.imendel'
    imendel = open(imendel_nam, 'r')
    nex = 0
    
    dumphead = imendel.readline()
    for line in imendel:
        (fid, iid, nmendel) = line.spilt()
        
        if float(nmendal) / float(nsnp) > args.id_mendal_th:
            id_out.write(str(fid) + ' ' + str(iid) + ' excessive_mendel_errors\n')
            nex += 1

    imendel.close()
    print 'Found %d individuals to exclude for excessive mendel errors' % nex


# filter sex check
if not args.skip_sex_check:
    
    sexcheck_nam = ind_stats + '.sexcheck'
    sexcheck_res = open(sexcheck_nam, 'r')
    sex_warn_out = open(str(args.out)+'.sexcheck_warnings.txt', 'w')
    nex = 0
    nwarn = 0

    dumphead = sexcheck_res.readline()
    for line in sexcheck_res:
        (fid, iid, pedsex, snpsex, status_txt, fsex) = line.split()
        if int(pedsex) == 1:
            if str(fsex) == "NA":
                sex_warn_out.write(fid + ' ' + iid + ' sexcheck_F_NA')
                nwarn += 1
            if float(fsex) < 0.5:
                id_out.write(fid + ' ' + iid + ' sexcheck_pedmale')
                nex += 1
            elif fsex < 0.8:
               sex_warn_out.write(fid + ' ' + iid + ' sexcheck_pedmale')
               nwarn += 1
               
        elif int(pedsex) == 2:
            if str(fsex) == "NA":
                sex_warn_out.write(fid + ' ' + iid + ' sexcheck_F_NA')
                nwarn += 1
            if float(fsex) > 0.5:
                id_out.write(fid + ' ' + iid + ' sexcheck_pedfemale')
                nex += 1
            elif fsex > 0.2:
               sex_warn_out.write(fid + ' ' + iid + ' sexcheck_pedfemale')
               nwarn += 1            
                     
        else:
            sex_warn_out.write(fid + ' ' + iid + ' sexcheck_ped_NA')
            nwarn += 1
    
    sexcheck_res.close()
    sex_warn_out.close()
    print 'Found %d individuals to exclude for sex check errors' % nex
    print 'Sex check warnings for %d additional individuals. See %s' % (nwarn, str(args.out)+'.sexcheck_warnings.txt')


id_out.close()


# remove QC failure IDs
idfilter_out = args.out+".qcind"

idfilter_out_str = [str(plinkx), 
                               "--bfile", prefilter_out,
                               "--remove", idout_nam,
                               "--silent",
                               "--allow-no-sex",
                               "--make-bed",
                               "--out", idfilter_out]
print ''
print 'Running: ' + ' '.join(idfilter_out_str)
subprocess.check_call(idfilter_out_str)



#############
print '\n...QCing SNPs...'
# - get QC metrics
# - find missingness failures
# - find differential missingness failures 
# - find mendel error rate failures
# - find hwe failures
# 
#############

# 6) QC SNPs
#     - missings
#     - differential call rate
#     - mendel errors
#     - hwe (case, control)
# 7) Resolve remaining mendel errors
# 8) Generate basic summary info
# 9) Clean up files


