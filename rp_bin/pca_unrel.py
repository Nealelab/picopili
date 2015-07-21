#! /usr/bin/env python

####################################
# pca_unrel.py
# written by Raymond Walters, July 2015
"""
Runs PCA for GWAS data with related individuals
"""
# Overview:
# 1) Input plink bed/bim/fam
# 2) Run strict QC on the input data
# 2) Define set of unrelated individuals using PRIMUS
# 3) Compute PCA on the unrelated set and projects to remainder
# 4) Plot projected PCs
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


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='pca_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parsergrid, parserqc, parserpca])
                    
# parser.set_defaults(x=y, z=2, etc)
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

if args.extra_ld_regions == None:
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
                        "-R", str('\"rusage[mem=4]\"'),
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


print '\n############\n'
if args.test_sub:
    print 'Completed script. See dummy submit commands above.\n\n' 
else:
    print 'Finished submitting jobs.\n'
    print 'See bjobs to track job status.\n'
    print 'Email will be sent when workflow completes.\n\n'

exit(0)
