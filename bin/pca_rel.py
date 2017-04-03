#! /usr/bin/env python

####################################
# pca_rel.py
# written by Raymond Walters, July 2015
"""
Runs PCA for GWAS data with related individuals
"""
# Overview:
# 1) Parse arguments
#    - Get arugments for each task
#    - Handle pass-through args
#    - Record settings (includes applied defaults)
# 2) Submit job to: Run strict QC on the input data
# 3) Submit job to: Run IMUS PCA workflow
#    - Define set of unrelated individuals using PRIMUS
#    - Compute PCA on the unrelated set and projects to remainder
#    - Plot projected PCs
# 4) Submit job to: Confirm completion
#
####################################


import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import argparse
import os
from math import ceil
from args_pca import parserbase, parsergrid, parserqc, parserpca
from args_qc import parsertech
from py_helpers import file_len, unbuffer_stdout
from blueprint import send_job
unbuffer_stdout()


# get directory containing current script
# (to get absolute path for script directory)
rp_bin = os.path.dirname(os.path.realpath(__file__))

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='pca_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parsergrid, parserqc, parserpca, parsertech])
                    
args = parser.parse_args()


# convert bool args to txt
if args.no_cleanup:
    clean_txt = '--no-cleanup'
else:
    clean_txt = ''

if args.keep_mhc:
    mhc_txt = '--keep-mhc'
else:
    mhc_txt = ''

if args.keep_chr8inv:
    chr8inv_txt = '--keep-chr8inv'
else:
    chr8inv_txt = ''

if args.extra_ld_regions == None or args.extra_ld_regions == "None":
    ldregions_txt = ''
else:
    ldregions_txt = str('--extra-ld-regions '+args.extra_ld_regions)

if args.keep_indels:
    indel_txt = '--keep-indels'
else:
    indel_txt = ''

if args.keep_strand_ambiguous:
    strandambi_txt = '--keep-strand-ambiguous'
else:
    strandambi_txt = ''

if args.all_chr:
    allchr_txt = '--all-chr'
else:
    allchr_txt = ''

if args.plot_all:
    plotall_txt = '--plot-all'
else:
    plotall_txt = ''


### print settings in use
print 'Basic settings:'
print '--bfile '+args.bfile
print '--out '+args.out

print '\nQC Thresholds:'
print '--mind-th '+str(args.mind_th)
print '--maf-th '+str(args.maf_th)
print '--hwe-th '+str(args.hwe_th)
print '--miss-th '+str(args.miss_th)

print '\nLD Pruning Parameters:'
print '--ld-th '+str(args.ld_th)
print '--ld_wind '+str(args.ld_wind)
print '--keep-mhc '+str(args.keep_mhc)
print '--keep-chr8inv '+str(args.keep_chr8inv)
print '--extra-ld-regions '+str(args.extra_ld_regions)

print '\nAdditional SNP Criteria:'
print '--keep-indels '+str(args.keep_indels)
print '--keep-strand-ambiguous '+str(args.keep_strand_ambiguous)
print '--all_chr '+str(args.all_chr)

print '\nUnrelated Set (IMUS) Criteria:'
print '--rel-th '+str(args.rel_th)

print '\nPrincipal Components (PCA):'
print '--npcs '+str(args.npcs) 
print '--plot-all '+str(args.plot_all)



#####
# check imus memory requirements
# = 6 GB + 400MB*(n/1000)^2, rounded up to nearest 4GB
# based on previous runs of PRIMUS
#####

warn_mem = False

nsamp = float(file_len(str(args.bfile)+'.fam'))

imus_mem = int(ceil( (6000.0+400.0*( (nsamp/1000.0)**2) )/4000.0 ) * 4)

if imus_mem > 16 and not args.large_mem_ok:
    warn_mem = True
    args.test_sub = True


#####
# submit strict qc
print '\n...Submitting Strict QC job...'
#####
strictqc_call = ' '.join(['strict_qc.py', 
                         '--bfile', args.bfile,
                         '--out', args.out,
                         clean_txt,
                         '--mind-th', str(args.mind_th),
                         '--maf-th', str(args.maf_th),
                         '--hwe-th', str(args.hwe_th),
                         '--miss-th', str(args.miss_th),
                         '--ld-th', str(args.ld_th),
                         '--ld-wind', str(args.ld_wind),
			 '--plink-mem',str(args.plink_mem),
                         mhc_txt,
                         chr8inv_txt,
                         ldregions_txt,
                         indel_txt,
                         strandambi_txt,
                         allchr_txt])

jobres = send_job(jobname=str('strictqc_'+args.out),
                  arrayfile=None,
                  cmd=str(strictqc_call),
                  logname=str('strictqc_'+args.out+'.sub.log'),
                  mem=int(args.plink_mem),
                  walltime=2,
                  sleep=0,
                  testonly=args.test_sub)


#####
# submit imus pca
print '\n...Submitting IMUS PCA job...'
#####

imuspca_call = ' '.join(['imus_pca.py',
                         '--bfile', str(args.out+'.strictqc.pruned'),
                         '--out', args.out,
                         clean_txt,
                         '--rel-th', str(args.rel_th),
                         '--npcs', str(args.npcs),
                         plotall_txt,
                         '--pcadir', str(args.pcadir),
                         '--rscript-ex', str(args.rscript_ex),
                         '--primus-ex', str(args.primus_ex)
                         ])

jobres2 = send_job(jobname=str('imuspca_'+args.out),
                   cmd=str(imuspca_call),
                   logname=str('imuspca_'+args.out+'.sub.log'),
                   mem=int(imus_mem)*1000,
                   walltime=30,
                   wait_name=str('strictqc_'+args.out),
                   wait_num=str(jobres).strip(),
                   sleep=args.sleep,
                   testonly=args.test_sub)


#####
# submitting final file check
print '\n...Submitting final file checker/notifier...'
#####

wd = os.getcwd()
if args.pcadir == None or args.pcadir == "None":
    pcaout = str(args.out + '_imus_pca')
else:
    pcaout = str(args.pcadir)


final_call = ' '.join(['final_file_check.py',
                       '--filename', str(wd+'/'+pcaout+'/plots/'+args.out+'.pca.pairs.png'),
                       '--taskname', str('pca_rel_'+args.out)])

send_job(jobname=str('checkfinal_'+args.out),
         arrayfile=None,
         cmd=str(final_call),
         logname=str('checkfinal_'+args.out+'.sub.log'),
         mem=100,
         walltime=1,
         wait_name=str('imuspca_'+args.out),
         wait_num=str(jobres2).strip(),
         sleep=str(args.sleep),
         testonly=args.test_sub)

#######
# Print completion message
# - flags if need rerun for large mem
# - depends on if actually submitted or just test run
#######
print '\n############\n'
if warn_mem:
    print 'WARNING: Commands not submitted!'
    print 'IMUS PCA expected to require > 16 GB of RAM (estimate %d GB based on sample size)' % imus_mem
    print 'Please resubmit with \'--large-mem-ok\' if this is acceptable.\n'
    exit(1)
elif args.test_sub:
    print 'Completed script. See dummy submit commands above.\n' 
else:
    print 'Finished submitting jobs.'
    print 'See qstat to track job status.'
    print 'Email will be sent when workflow completes.\n'

exit(0)
