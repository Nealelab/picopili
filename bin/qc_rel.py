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
from time import strftime
start_time = strftime("%H:%M:%S %d-%B-%Y")
# from glob import glob
from args_qc import *
from py_helpers import unbuffer_stdout, read_conf, test_exec, link, file_len, warn_format, find_exec
unbuffer_stdout()
warnings.formatwarning = warn_format

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
if args.no_cleanup:
    print '--no-cleanup '+str(args.no_cleanup)

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
print '--min-chrx-snps '+str(args.min_chrx_snps)

print '\nSNP QC Thresholds:'
print '--pre-miss '+str(args.pre_miss)
print '--miss-th '+str(args.miss_th)
print '--diff-miss-abs '+str(args.diff_miss_abs)
print '--diff-miss-p '+str(args.diff_miss_p)
print '--hwe-th-cas '+str(args.hwe_th_cas)
print '--hwe-th-con '+str(args.hwe_th_con)
print '--hwe-th-all '+str(args.hwe_th_all)
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
print '\n...Checking dependencies...'
# check exists, executable
#############

### read config
conf_file = os.environ['HOME']+"/picopili.conf"
configs = read_conf(conf_file)
analyst = configs['init']

# find plink
plinkx = find_exec('plink',key='p2loc')

if not args.skip_platform:
    # get directory containing current script
    # (hack to get plague script location)
    rp_bin = os.path.dirname(os.path.realpath(__file__))
    plague_ex = rp_bin + '/plague_pico.pl'
    test_exec(plague_ex, 'Platform guessing script')
# TODO: verify plague works properly across platforms (primary concern is Compress::Zlib loading)

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
                               "--memory", str(2000),
                               "--allow-no-sex",
                               "--out", prefilter_stats]

print 'Running: ' + ' '.join(prefilter_stats_str)
subprocess.check_call(prefilter_stats_str)

# get failing SNPs from plink lmiss file
snp_out = open(snpout_nam, 'w')

prefilter_lmiss_nam = prefilter_stats + '.lmiss'
prefilter_lmiss = open(prefilter_lmiss_nam, 'r')
nex = 0
nsnp_pre = 0 # total number of SNPs pre-QC

dumphead = prefilter_lmiss.readline()
for line in prefilter_lmiss:
    (chrom, snp, nmiss, ngeno, fmiss) = line.split()
    nsnp_pre += 1
    
    if float(fmiss) > args.pre_miss:
        snp_out.write(snp + ' prefilter_high_missing\n')
        nex += 1

prefilter_lmiss.close()
snp_out.close()
print 'Found %d SNPs to exclude for failing pre-filter on missingness rate > %r' % (nex, args.pre_miss)


# remove SNPs
prefilter_out = args.out+".prefiltered"

prefilter_out_str = [str(plinkx), 
                               "--bfile", args.bfile,
                               "--exclude", snpout_nam,
                               "--silent",
                               "--memory", str(2000),
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
        if str(chrom) == '23' or str(chrom) == 'X':
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
    mendel_txt = ['','']
    
elif args.mendel == 'trios':
    mendel_txt = ['--mendel','']
    
elif args.mendel == 'duos':
    mendel_txt = ['--mendel','--mendel-duos']
    
elif args.mendel == 'multigen':
    mendel_txt = ['--mendel','--mendel-multigen']
    
else:
    # shouldn't be possible from argparse
    raise ValueError ("Invalid argument for '--mendel' (%s)" % str(args.mendel))


# get plink sum stats
ind_stats_str = [str(plinkx), 
                         "--bfile", prefilter_out,
                         "--missing",
                         "--het",
                         sexcheck_txt,
                         mendel_txt[0], mendel_txt[1],
                         "--read-freq",str(prefilter_stats)+'.frqx',
                         "--silent",
                         "--memory", str(2000),
                         "--allow-no-sex",
                         "--out", ind_stats]

# remove empty elements
ind_stats_str = filter(None, ind_stats_str)

print 'Running: ' + ' '.join(ind_stats_str)
subprocess.check_call(ind_stats_str)


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
print 'Found %d individuals to exclude for failing missingness rate > %r' % (nex, args.id_mendel_th)


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
print 'Found %d individuals to exclude for failing absolute Fhet homozygosity rate > %r' % (nex, args.het_th)


# filter mendel errors
if args.mendel != 'none':
    
    # get number of SNPs as denominator for mendel error rate
    nsnp = file_len(str(prefilter_out) + '.bim')

    imendel_nam = ind_stats + '.imendel'
    imendel = open(imendel_nam, 'r')
    nex = 0
    
    dumphead = imendel.readline()
    for line in imendel:
        (fid, iid, nmendel) = line.split()
        
        if float(nmendel) / float(nsnp) > args.id_mendel_th:
            id_out.write(str(fid) + ' ' + str(iid) + ' excessive_mendel_errors\n')
            nex += 1

    imendel.close()
    print 'Found %d individuals to exclude for excessive mendel errors > %r of SNPs' % (nex, args.id_mendel_th)


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
                sex_warn_out.write(fid + ' ' + iid + ' sexcheck_F_NA\n')
                nwarn += 1
            if float(fsex) < 0.5:
                id_out.write(fid + ' ' + iid + ' sexcheck_pedmale\n')
                nex += 1
            elif float(fsex) < 0.8:
               sex_warn_out.write(fid + ' ' + iid + ' sexcheck_pedmale\n')
               nwarn += 1
               
        elif int(pedsex) == 2:
            if str(fsex) == "NA":
                sex_warn_out.write(fid + ' ' + iid + ' sexcheck_F_NA\n')
                nwarn += 1
            if float(fsex) > 0.5:
                id_out.write(fid + ' ' + iid + ' sexcheck_pedfemale\n')
                nex += 1
            elif float(fsex) > 0.2:
               sex_warn_out.write(fid + ' ' + iid + ' sexcheck_pedfemale\n')
               nwarn += 1            
                     
        else:
            sex_warn_out.write(fid + ' ' + iid + ' sexcheck_ped_NA\n')
            nwarn += 1
    
    sexcheck_res.close()
    sex_warn_out.close()
    print 'Found %d individuals to exclude for sex check errors' % nex
    print 'Sex check warnings for %d additional individuals. See %s' % (nwarn, str(args.out)+'.sexcheck_warnings.txt')


id_out.close()


# remove QC failure IDs
idfilter_out = args.out+".qc_ind"

idfilter_out_str = [str(plinkx), 
                               "--bfile", prefilter_out,
                               "--remove", idout_nam,
                               "--silent",
                               "--memory", str(2000),
                               "--allow-no-sex",
                               "--make-bed",
                               "--out", idfilter_out]
print ''
print 'Running: ' + ' '.join(idfilter_out_str)
subprocess.check_call(idfilter_out_str)


#############
# Summarizing per-individual QC
#############
# parse initial fam file for summary info
n_pre = 0
n_pre_cas = 0
n_pre_con = 0
n_pre_male = 0
n_pre_female = 0
n_pre_fam = 0
pre_fids = []

if args.skip_fid_tags:
    pre_fam = open(str(args.bfile)+'.fam', 'r')
else:
    pre_fam = open(str(args.bfile)+'.fam.original', 'r')

for line in pre_fam:
    (fid, iid, pat, mat, sex, phen) = line.split()
    
    n_pre += 1

    if str(sex) == '2':
        n_pre_female += 1
    elif str(sex) == '1':
        n_pre_male += 1
    
    if str(phen) == '2':
        n_pre_cas += 1
    elif str(phen) == '1':
        n_pre_con += 1
    
    if fid not in pre_fids:
        pre_fids.append(fid)

pre_fam.close()
n_pre_fam = len(pre_fids)


# parse final fam file for summary info
n_post = 0
n_post_cas = 0
n_post_con = 0
n_post_male = 0
n_post_female = 0
n_post_fam = 0
post_fids = []

post_fam = open(str(idfilter_out)+'.fam', 'r')
for line in post_fam:
    (fid, iid, pat, mat, sex, phen) = line.split()
    
    n_post += 1

    if str(sex) == '2':
        n_post_female += 1
    elif str(sex) == '1':
        n_post_male += 1
    
    if str(phen) == '2':
        n_post_cas += 1
    elif str(phen) == '1':
        n_post_con += 1
    
    if fid not in post_fids:
        post_fids.append(fid)

post_fam.close()
n_post_fam = len(post_fids)



# check if have cases and controls for differential missingness
if not args.skip_diff_miss:
    if n_post_cas == 0 or n_post_con == 0:
        warnings.warn('WARNING: No cases and/or controls remaining after per-ID QC. Will skip differential missingness checks.')
        args.skip_diff_miss = True


#############
print '\n...QCing SNPs...'
# - get QC metrics
# - find missingness failures
# - find differential missingness failures 
# - find mendel error rate failures
# - find hwe failures
# - find maf failures
# 
#############

snp_stats = str(args.out) + '.snp_stats'

if args.skip_diff_miss:
    diff_miss_txt = ''
else:
    diff_miss_txt = '--test-missing'

# get metrics
snp_stats_str = [str(plinkx), 
                         "--bfile", str(idfilter_out),
                         "--missing",
                         diff_miss_txt,
                         mendel_txt[0], mendel_txt[1],
                         "--make-founders","require-2-missing",
                         "--freq",
                         "--hardy",
                         "--silent",
                         "--memory", str(2000),
                         "--allow-no-sex",
                         "--out", snp_stats]

# remove empty elements
snp_stats_str = filter(None, snp_stats_str)

print 'Running: ' + ' '.join(snp_stats_str)
subprocess.check_call(snp_stats_str)


snp_out = open(snpout_nam, 'a')


# filter missingness
lmiss_nam = snp_stats + '.lmiss'
lmiss = open(lmiss_nam, 'r')
dumphead = lmiss.readline()
nex = 0

for line in lmiss:
    (chrom, snp, nmiss, ngeno, fmiss) = line.split()
    
    if float(fmiss) > args.miss_th:
        snp_out.write(snp + ' high_missing\n')
        nex += 1

lmiss.close()
print 'Found %d SNPs to exclude for missingness rate > %r' % (nex, args.miss_th) 


# filter differential missingness
nex_abs = 0
nex_p = 0

if not args.skip_diff_miss:
    diffmiss_nam = str(snp_stats) + '.missing'
    diffmiss = open(diffmiss_nam, 'r')
    dumphead = diffmiss.readline()
    
    
    for line in diffmiss:
        (chrom, snp, miss_cas, miss_con, miss_p) = line.split()
    
        diff = float(miss_cas) - float(miss_con)
    
        if float(miss_p) < args.diff_miss_p:
            snp_out.write(snp + ' differential_missing_pval\n')
            nex_p += 1
            
        if abs(float(diff)) > args.diff_miss_abs:
            snp_out.write(snp + ' absolute_differential_missing\n')
            nex_abs += 1
            
    diffmiss.close()
    print 'Found %d SNPs to exclude for absolute differential missingness > %r' % (nex_abs, args.diff_miss_abs)
    print 'Found %d SNPs to exclude for differential missingness p-value < %r' % (nex_p, args.diff_miss_p)


# filter mendel errors
if args.mendel != 'none':

    # get number of families as denominator for mendel error rate
    nind = file_len(str(snp_stats) + '.fmendel')

    lmendel_nam = str(snp_stats) + '.lmendel'
    lmendel = open(lmendel_nam, 'r')
    nex = 0
    
    dumphead = lmendel.readline()
    for line in lmendel:
        (chrom, snp, nmendel) = line.split()
        
        if float(nmendel) / float(nind) > args.snp_mendel_th:
            snp_out.write(str(snp) + ' excessive_mendel_errors\n')
            nex += 1

    lmendel.close()
    print 'Found %d SNPs to exclude for excessive mendel errors > %r of nuclear families' % (nex, args.snp_mendel_th)


# filter HWE
hwe_nam = str(snp_stats) + '.hwe'
hwe = open(hwe_nam, 'r')
dumphead = hwe.readline()
nex_all = 0
nex_cas = 0
nex_con = 0

for line in hwe:
        (chrom, snp, test, a1, a2, geno, ohet, ehet, hwe_p) = line.split()
        
        if test == 'ALL' or test == 'ALL(QT)' or test == 'ALL(NP)':
            if hwe_p != 'NA' and float(hwe_p) < args.hwe_th_all:
                snp_out.write(str(snp) + ' hardy-weinberg_all\n')
                nex_all += 1
                
        elif test == 'AFF':
            if hwe_p != 'NA' and float(hwe_p) < args.hwe_th_cas:
                snp_out.write(str(snp) + ' hardy-weinberg_cases\n')
                nex_cas += 1
                
        elif test == 'UNAFF':
            if hwe_p != 'NA' and float(hwe_p) < args.hwe_th_con:
                snp_out.write(str(snp) + ' hardy-weinberg_controls\n')
                nex_con += 1
            
        else:
            raise IOError('Failed to parse Hardy-Weinberg results for SNP %s in %s' % (str(snp), str(hwe_nam)))

hwe.close()
print 'Found %d SNPs to exclude for Hardy-Weinberg p-value in founder cases > %r' % (nex_cas, args.hwe_th_cas)
print 'Found %d SNPs to exclude for Hardy-Weinberg p-value in founder controls > %r' % (nex_con, args.hwe_th_con)
print 'Found %d SNPs to exclude for Hardy-Weinberg p-value in all founder IDs > %r' % (nex_all, args.hwe_th_all)


# filter MAF
if float(args.maf_th) >= 0.0:

    maf_nam = str(snp_stats) + '.frq'
    maf = open(maf_nam, 'r')
    dumphead = maf.readline()
    nex_maf = 0
    nex_mis = 0
    
    for line in maf:
        (chrom, snp, a1, a2, a1_freq, nchrobs) = line.split()
    
        if str(a1_freq) is 'NA':
            snp_out.write(str(snp) + ' maf_missing\n')
            nex_mis += 1
            
        else:
            min_freq = min(float(a1_freq), 1.0 - float(a1_freq))
            if min_freq <= float(args.maf_th):
                snp_out.write(str(snp) + ' low_maf\n')
                nex_maf += 1
        
    maf.close()
    print 'Found %d SNPs to exclude for missing MAF in founders' % nex_mis
    print 'Found %d SNPs to exclude for low/invariant MAF <= %r in founders' % (nex_maf, args.maf_th)
    


# remove QC failure SNPs
snp_out.close()
snpfilter_out = args.out+".qc_snp"

snpfilter_out_str = [str(plinkx), 
                               "--bfile", str(idfilter_out),
                               "--exclude", snpout_nam,
                               "--silent",
                               "--memory", str(2000),
                               "--allow-no-sex",
                               "--make-bed",
                               "--out", snpfilter_out]
print ''
print 'Running: ' + ' '.join(snpfilter_out_str)
subprocess.check_call(snpfilter_out_str)



#############
# Resolve remaining mendel errors
#############

if (args.mendel != 'none') and (not args.keep_mendel):
    print '\n...Zeroing out remaining mendelian errors...'

    mendel_out = args.out + '.nomendel'
    
    drop_mendel_str = [str(plinkx),
                           '--bfile', snpfilter_out,
                           mendel_txt[1],
                           '--set-me-missing',
                           '--silent',
                           "--memory", str(2000),
                           '--allow-no-sex',
                           '--make-bed',
                           '--out', mendel_out]
    
    print 'Running: ' + ' '.join(drop_mendel_str)
    subprocess.check_call(drop_mendel_str)


#############
print '\n...Linking final QCed files...'
# - link final file to standardized output name ala ricopili
# - different final file if did mendel error removal  
#############

if args.skip_fid_tags:
    rp_outname = str(args.out)+'-qc'
else:
    rp_outname = str(args.disname)+'_'+str(args.out)+'_'+str(args.popname)+'_'+str(analyst)+'-qc'

os.chdir(wd)

if (args.mendel != 'none') and (not args.keep_mendel):
    link('./'+qcdir+'/'+str(mendel_out)+'.bed', rp_outname+'.bed', 'final QCed .bed file')
    link('./'+qcdir+'/'+str(mendel_out)+'.bim', rp_outname+'.bim', 'final QCed .bim file')
    link('./'+qcdir+'/'+str(mendel_out)+'.fam', rp_outname+'.fam', 'final QCed .fam file')
    
else:
    link('./'+qcdir+'/'+str(snpfilter_out)+'.bed', rp_outname+'.bed', 'final QCed .bed file')
    link('./'+qcdir+'/'+str(snpfilter_out)+'.bim', rp_outname+'.bim', 'final QCed .bim file')
    link('./'+qcdir+'/'+str(snpfilter_out)+'.fam', rp_outname+'.fam', 'final QCed .fam file')  
    
print 'Output files: %s.bed (bim, fam)' % rp_outname

os.chdir(qcdir)


#############
print '\n...Generating summary information...'
# Summary file
# - input file
# - output file
# - time
# - args
# - N, cas/con, m/f, Nfam pre/post
# - num SNPs pre/post
# - Lambda, Lambda1000 pre/post
# QQ plot pre/post 
#############

# TODO: GWAS
# TODO: compute lambda
# TODO: QQ plot




### print file of summary info
sum_file_nam = str(args.out)+'.qc_summary.txt'
sum_file = open(wd+'/'+sum_file_nam, 'w')

sum_file.write('### Files:\n')
sum_file.write('Input bed: %s' % str(wd)+'/'+args.bfile+'.bed\n')
sum_file.write('Output bed: %s' % rp_outname+'.bed\n')
sum_file.write('ID exclusions: %s\n' % idout_nam)
sum_file.write('SNP exclusions: %s\n' % snpout_nam)
sum_file.write('\n')

sum_file.write('Individual QC:\n')
sum_file.write('Pre QC:  %d individuals (%d cases, %d controls, %d male, %d female, %d families)\n' % (n_pre, n_pre_cas, n_pre_con, n_pre_male, n_pre_female, n_pre_fam))
sum_file.write('Post QC: %d individuals (%d cases, %d controls, %d male, %d female, %d families)\n' % (n_post, n_post_cas, n_post_con, n_post_male, n_post_female, n_post_fam))
sum_file.write('\n')

sum_file.write('### SNP QC:\n')
sum_file.write('Pre QC:  %d SNPs\n' % nsnp_pre)
nsnp_post = file_len(wd+'/'+rp_outname+'.bim')
sum_file.write('Post QC: %d SNPs\n' % nsnp_post)
sum_file.write('\n')

# sum_file.write('Genomic Control:\n')
# sum_file.write('Pre QC:  lambda %r, Lambda1000 %r\n' % (lam, lam1k))
# sum_file.write('Post QC: lambda %r, Lambda1000 %r\n' % (lam, lam1k))
# sum_file.write('\n')

sum_file.write('### Arguments:\n')
sum_file.write('Basic Settings:'+'\n')
sum_file.write('--bfile '+args.bfile+'\n')
sum_file.write('--out '+args.out+'\n')
if args.no_cleanup:
    sum_file.write('--no-cleanup '+str(args.no_cleanup)+'\n')

sum_file.write('\nID Tagging Information:'+'\n')
sum_file.write('--skip-fid-tags '+str(args.skip_fid_tags)+'\n')
if not args.skip_fid_tags:
   sum_file.write('--disname '+str(args.disname)+'\n')
   sum_file.write('--popname '+str(args.popname)+'\n')
   sum_file.write('--skip-platform '+str(args.skip_platform)+'\n')

sum_file.write('\nIndividual QC Thresholds:'+'\n')
sum_file.write('--mind-th '+str(args.mind_th)+'\n')
sum_file.write('--het-th '+str(args.het_th)+'\n')
sum_file.write('--skip-sex-check '+str(args.skip_sex_check)+'\n')

sum_file.write('\nSNP QC Thresholds:'+'\n')
sum_file.write('--pre-miss '+str(args.pre_miss)+'\n')
sum_file.write('--miss-th '+str(args.miss_th)+'\n')
sum_file.write('--skip-diff-miss '+str(args.skip_diff_miss)+'\n')
if not args.skip_diff_miss:
    sum_file.write('--diff-miss-abs '+str(args.diff_miss_abs)+'\n')
    sum_file.write('--diff-miss-p '+str(args.diff_miss_p)+'\n')
sum_file.write('--hwe-th-cas '+str(args.hwe_th_cas)+'\n')
sum_file.write('--hwe-th-con '+str(args.hwe_th_con)+'\n')
sum_file.write('--hwe-th-all '+str(args.hwe_th_all)+'\n')
if float(args.maf_th) >= 0.0:
   sum_file.write('--maf-th '+str(args.maf_th)+'\n')

sum_file.write('\nMendelian Error Checks:'+'\n')
sum_file.write('--mendel '+str(args.mendel)+'\n')
if str(args.mendel) != 'none':
    sum_file.write('--id-mendel-th '+str(args.id_mendel_th)+'\n')
    sum_file.write('--snp-mendel-th '+str(args.snp_mendel_th)+'\n')
    sum_file.write('--keep-mendel '+str(args.keep_mendel)+'\n')
sum_file.write('\n')

sum_file.write('### Time:\n')
sum_file.write('Start: %s\n' % start_time)
end_time = strftime("%H:%M:%S %d-%B-%Y")
sum_file.write('Finish: %s\n' % end_time)
sum_file.write('\n')

sum_file.close()
print 'QC summary written to: %s' % sum_file_nam


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
    print 'Zipping interim plink files:'
    #############
    
    # pre-filtered plink files
    subprocess.check_call(["tar", "-zcvf",
                           args.out + '.prefiltered.plink.tar.gz',
                           args.out + '.prefiltered.bed',
                           args.out + '.prefiltered.bim',
                           args.out + '.prefiltered.fam'])   
    # remove files once successfully zipped
    subprocess.check_call(['rm',
                           args.out + '.prefiltered.bed',
                           args.out + '.prefiltered.bim',
                           args.out + '.prefiltered.fam'])

    # ID QC plink files
    subprocess.check_call(["tar", "-zcvf",
                           args.out + '.qc_ind.plink.tar.gz',
                           args.out + '.qc_ind.bed',
                           args.out + '.qc_ind.bim',
                           args.out + '.qc_ind.fam'])
    subprocess.check_call(['rm',
                           args.out + '.qc_ind.bed',
                           args.out + '.qc_ind.bim',
                           args.out + '.qc_ind.fam'])
                           
    if str(args.mendel) != 'none' and (not args.keep_mendel):
        subprocess.check_call(["tar", "-zcvf",
                               args.out + '.qc_snp.plink.tar.gz',
                               args.out + '.qc_snp.bed',
                               args.out + '.qc_snp.bim',
                               args.out + '.qc_snp.fam'])
        subprocess.check_call(['rm',
                               args.out + '.qc_snp.bed',
                               args.out + '.qc_snp.bim',
                               args.out + '.qc_snp.fam'])                               
    

    #############
    print '\nZipping QC metrics:'
    #############
    subprocess.check_call(["tar", "-zcvf",
                           args.out + '.prefilter_stats.tar.gz',
                           args.out + '.prefilter_stats.imiss',
                           args.out + '.prefilter_stats.lmiss',
                           args.out + '.prefilter_stats.frqx'])   

    ind_stats_list = [args.out + '.ind_stats.imiss',
                      args.out + '.ind_stats.lmiss',
                      args.out + '.ind_stats.het']

    if not args.skip_sex_check:
        ind_stats_list.append(args.out + '.ind_stats.sexcheck')

    if str(args.mendel) != 'none':
        ind_stats_list.extend([args.out+'.ind_stats.mendel', 
                               args.out+'.ind_stats.imendel', 
                               args.out+'.ind_stats.lmendel', 
                               args.out+'.ind_stats.fmendel'])
                               
    zip_ind_cmd = ["tar", "-zcvf", args.out+'.ind_stats.tar.gz']
    zip_ind_cmd.extend(ind_stats_list)
    subprocess.check_call(zip_ind_cmd)
    
    
    

    snp_stats_list = [args.out + '.snp_stats.frq',
                      args.out + '.snp_stats.imiss',
                      args.out + '.snp_stats.lmiss',
                      args.out + '.snp_stats.hwe']
    
    if not args.skip_diff_miss:
        snp_stats_list.extend([args.out + '.snp_stats.missing'])
    
    if str(args.mendel) != 'none':    
        snp_stats_list.extend([args.out+'.snp_stats.mendel', 
                               args.out+'.snp_stats.imendel', 
                               args.out+'.snp_stats.lmendel', 
                               args.out+'.snp_stats.fmendel'])

    zip_snp_cmd = ["tar", "-zcvf", args.out+'.snp_stats.tar.gz']
    zip_snp_cmd.extend(snp_stats_list)
    subprocess.check_call(zip_snp_cmd)


    # remove files once zipped
    subprocess.check_call(['rm',
                           args.out + '.prefilter_stats.imiss',
                           args.out + '.prefilter_stats.lmiss',
                           args.out + '.prefilter_stats.frqx'])
    
    rm_ind_cmd = ['rm']
    rm_ind_cmd.extend(ind_stats_list)
    subprocess.check_call(rm_ind_cmd)
    
    rm_snp_cmd = ['rm']
    rm_snp_cmd.extend(snp_stats_list)
    subprocess.check_call(rm_snp_cmd)    
    
    #############
    print '\nCompress exclusion/warning lists:'
    #############
    subprocess.check_call(["gzip", "-v", snpout_nam])
    subprocess.check_call(["gzip", "-v", idout_nam])

    if not args.skip_sex_check:
        subprocess.check_call(["gzip", "-v", str(args.out)+'.sexcheck_warnings.txt'])
    
    
    #############
    print '\nRemove if exist:'
    #############
    # allowing failure, since files may or may not exists
    subprocess.call(["rm", "-v",
                     args.out + '.prefilter_stats.hh',
                     args.out + '.prefiltered.hh',
                     args.out + '.ind_stats.hh',
                     args.out + '.qc_ind.hh',
                     args.out + '.snp_stats.hh',
                     args.out + '.qc_snp.hh',
                     args.out + '.nomendel.hh',
                     ])


print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)
