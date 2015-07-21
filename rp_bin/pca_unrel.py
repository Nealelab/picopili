#! /usr/bin/env python

####################################
# pca_unrel.py
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
import subprocess
from args_pca import *
from py_helpers import file_len


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='pca_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parsergrid, parserqc, parserpca])
                    
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
print '\n'

print 'QC Thresholds:'
print '--mind-th '+str(args.mind_th)
print '--maf-th '+str(args.maf_th)
print '--hwe-th '+str(args.hwe_th)
print '--miss-th '+str(args.miss_th)
print '\n'

print 'LD Pruning Parameters:'
print '--ld-th '+str(args.ld_th)
print '--ld_wind '+str(args.ld_wind)
print '--keep-mhc '+str(args.keep_mhc)
print '--keep-chr8inv '+str(args.keep_chr8inv)
print '--extra-ld-regions '+str(args.extra_ld_regions)
print '\n'

print 'Additional SNP Criteria:'
print '--keep-indels '+str(args.keep_indels)
print '--keep-strand-ambiguous '+str(args.keep_strand_ambiguous)
print '--all_chr '+str(args.all_chr)
print '\n'

print 'Unrelated Set (IMUS) Criteria:'
print '--rel-deg '+str(args.rel_deg)
print '\n'

print 'Principal Components (PCA):'
print '--npcs '+str(args.npcs) 
print '--plot-all '+str(args.plot_all)
print '\n'



#####
# check imus memory requirements
# = 6 GB + 400MB*(n/1000)^2, rounded up to nearest 4GB
# based on previous runs of PRIMUS
#####

warn_mem = False

nsamp = float(file_len(str(args.bfile)+'.fam'))

imus_mem = int(ceil( (6000.0+400.0*( (nsamp/1000.0)**2) )/4000.0 ) * 4000)

if imus_mem > 16000 and not args.large_mem_ok:
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
                         mhc_txt,
                         chr8inv_txt,
                         ldregions_txt,
                         indel_txt,
                         strandambi_txt,
                         allchr_txt])

strictqc_lsf = ' '.join(["bsub",
                         "-q", 'hour',
                         "-R", str('\"rusage[mem=2]\"'),
                         "-J", str('strictqc_'+args.out),
                         "-P", str('pico_'+args.out),
                         "-o", str('strictqc_'+args.out+'.bsub.log'),
                         "-r",
                         str('\"'+strictqc_call+'\"')])

print strictqc_lsf
if not args.test_sub:
    subprocess.check_call(strictqc_lsf, shell=True)


#####
# submit imus pca
print '\n...Submitting IMUS PCA job...'
#####

imuspca_call = ' '.join(['imus_pca.py',
                         '--bfile', str(args.out+'.strictqc.pruned'),
                         '--out', args.out,
                         clean_txt,
                         '--rel-deg', str(args.rel_deg),
                         '--npcs', str(args.npcs),
                         plotall_txt,
                         '--pcadir', str(args.pcadir),
                         '--rscript-ex', str(args.rscript_ex),
                         '--primus-ex', str(args.primus_ex)
                         ])

imuspca_lsf = ' '.join(["bsub",
                        "-w", str('\'ended(\"'+str('strictqc_'+args.out)+'\")\''),
                        "-E", str('\"sleep '+str(args.sleep)+'\"'),
                        "-q", 'week',
                        "-R", str('\"rusage[mem='+str(imus_mem)+']\"'),
                        "-J", str('imuspca_'+args.out),
                        "-P", str('pico_'+args.out),
                        "-o", str('imuspca_'+args.out+'.bsub.log'),
                        "-r",
                        str('\"'+imuspca_call+'\"')])

print imuspca_lsf
if not args.test_sub:
    subprocess.check_call(imuspca_lsf, shell=True)


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

final_lsf = ' '.join(["bsub",
                        "-w", str('\'ended(\"'+str('imuspca_'+args.out)+'\")\''),
                        "-E", str('\"sleep '+str(args.sleep)+'\"'),
                        "-q", 'hour',
                        "-J", str('checkfinal_'+args.out),
                        "-P", str('pico_'+args.out),
                        "-o", str('checkfinal_'+args.out+'.bsub.log'),
                        "-r",
                        str('\"'+final_call+'\"')])

print final_lsf
if not args.test_sub:
    subprocess.check_call(final_lsf, shell=True)


#######
# Print completion message
# - flags if need rerun for large mem
# - depends on if actually submitted or just test run
#######
print '\n############\n'
if warn_mem:
    print 'WARNING: Commands not submitted!'
    print 'IMUS PCA expected to require > 16 GB of RAM (estimate %d MB based on sample size)' % imus_mem
    print 'Please resubmit with \'--large-mem-ok\' if this is acceptable.\n'
    exit(1)
elif args.test_sub:
    print 'Completed script. See dummy submit commands above.\n' 
else:
    print 'Finished submitting jobs.'
    print 'See bjobs to track job status.'
    print 'Email will be sent when workflow completes.\n'

exit(0)
