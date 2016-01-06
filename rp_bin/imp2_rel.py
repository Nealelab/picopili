#! /usr/bin/env python

####################################
# imp2_rel.py
# written by Raymond Walters, January 2016
"""
Runs IMPUTE2 for GWAS data with related individuals
"""
# Overview:
# 1) Parse arguments
#    - check dependencies, print args, etc
# 2) Define genomic chunks
# 3) Run impute2
#
# TODO: better reference specification options
# TODO: chrx imputation
# 
####################################

import os
import subprocess
from args_impute import *
from py_helpers import unbuffer_stdout, file_len, link, read_conf
unbuffer_stdout()


parser = argparse.ArgumentParser(prog='imp2_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parserimpute, parserref, parserchunk, parsercluster])
                    
args = parser.parse_args()

# TODO: allow options
args.refstem = '/humgen/atgu1/fs03/shared_resources/1kG/shapeit/1000GP_Phase3'



# get useful modified args
if args.addout is not None and str(args.addout) != '':
    addout_txt = ['--addout',str(args.addout)]
    outdot = str(args.out)+'.'+str(args.addout)
else:
    addout_txt = ['','']
    outdot = str(args.out)
    
if args.seed is not None and str(args.seed) != '' and int(args.seed) > 0:
    seedtxt = '-seed '+str(args.seed)
else:
    seedtxt = ''


# report settings in use
print '\nBasic settings:'
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
if args.addout is not None:
    print '--addout '+str(args.addout)

print '\nIMPUTE2 arguments:'
print '--ne '+str(args.ne)
print '--buffer '+str(args.buffer)
if args.seed is not None and str(args.seed) != '' and int(args.seed) > 0:
    print '--seed '+str(args.seed)

print '\nImputation reference:'
print '--refstem '+str(args.refstem)
print '--map-dir '+str(args.map_dir)

print '\nGenomic chunks:'
print '--Mb-size '+str(args.Mb_size)
print '--snp_size '+str(args.snp_size)
print '--chr_info_file '+str(args.chr_info_file)

print '\nCluster settings:'
print '--sleep '+str(args.sleep)



#############
print '\n...Reading ricopili config file...'
#############

### read plink loc from config
conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

impute_ex = configs['i2loc']+"impute2"

# get directory containing current script
# (to get absolute path for scripts)
rp_bin = os.path.dirname(os.path.realpath(__file__))
chunker_ex = rp_bin+'/chunk_snps.py'



#############
print '\n...Checking dependencies...'
#############



# TODO: here
# .hg19.ch.fl.bim for chunking
# imp. references
# executables




print '\n'
print '############'
print 'Begin!'
print '############'


######################
print '\n...Setting up working directories...'
######################

# working directories
wd = os.getcwd()
shape_dir = wd + '/phase_chr'
chunk_dir = wd + '/chunks_for_imp'
imp_dir = wd + '/imp_sub'
os.mkdir(chunk_dir)
os.mkdir(imp_dir)
print 'Genomic chunks: %s' % chunk_dir
print 'Imputation: %s' % imp_dir



######################
print '\n...Creating genomic chunks for imputation...'
######################

os.chdir(chunk_dir)
link(str(shape_dir)+'/'+str(args.bfile)+'.hg19.ch.fl.bim', str(args.bfile)+'.hg19.ch.fl.bim', 'reference-aligned bim file')

# create chunks
chunk_call = [chunker_ex,
              '--bfile',str(args.bfile)+'.hg19.ch.fl',
              '--out',str(args.out),
              addout_txt[0],addout_txt[1],
              '--Mb-size',str(args.Mb_size),
              '--snp-size',str(args.snp_size),
              '--chr-info-file',str(args.chr_info_file)]
chunk_call = filter(None,chunk_call)

chunk_log = open('chunk.'+str(outdot)+'.log', 'w')
print ' '.join(chunk_call) + '\n'
subprocess.check_call(chunk_call, stderr=subprocess.STDOUT, stdout=chunk_log)
chunk_log.close()



######################
print '\n...Submitting IMPUTE2 job array...'
# - require "allow_large_regions" to handle centromere boundaries
# - write submit script to include chunk name parsing
# TODO: consider making queue/resources flexible
#
# impute2 
# -use_prephased_g 
# -known_haps_g chr20.phased.haps 
# -h chr20.reference.hap.gz 
# -l chr20.reference.legend.gz 
# -m chr20.gmap.gz 
# -int 35e6 36e6 
# -Ne 20000
# -allow_large_regions
# -o_gz
# -o chr20.imputed
# -seed 12345
#
######################

os.chdir(imp_dir)
link(str(chunk_dir)+'/'+str(outdot)+'.chunks.txt', str(outdot)+'.chunks.txt', 'genomic chunk results')

uger_imp_template = """#!/usr/bin/env sh
#$ -j y
#$ -cwd
#$ -V
#$ -N {jname}
#$ -q short
#$ -l m_mem_free=8g
#$ -t 1-{nchunk}
#$ -o {outlog}

source /broad/software/scripts/useuse
reuse -q Anaconda
sleep {sleep}

cchr=`awk -v a=${{SGE_TASK_ID}} 'NR==a+1{{print $1}}' {cfile}`
cstart=`awk -v a=${{SGE_TASK_ID}} 'NR==a+1{{print $2}}' {cfile}`
cend=`awk -v a=${{SGE_TASK_ID}} 'NR==a+1{{print $3}}' {cfile}`
cname=`awk -v a=${{SGE_TASK_ID}} 'NR==a+1{{print $4}}' {cfile}`

{impute_ex} -use_prephased_g -known_haps_g {in_haps} -h {ref_haps} -l {ref_leg} -m {map} -int ${{cstart}} ${{cend}} -buffer {buffer} -Ne {Ne} -allow_large_regions -o_gz -o {out} {seedtxt}

# eof
"""
    
# get number of chunks (-1 is for header)
nchunks = file_len(outdot+'.chunks.txt')-1

# fill in template
jobdict = {"jname": 'imp.chunks.'+str(outdot),
           "nchunk": str(nchunks),
           "outlog": str('imp.chunks.'+str(outdot)+'.$TASK_ID.qsub.log'),
           "sleep": str(args.sleep),
           "cfile": str(outdot)+'.chunks.txt',
           "impute_ex": str(impute_ex),
           "in_haps": str(shape_dir)+'/'+str(outdot)+'.chr${cchr}.phased.haps',
           "ref_haps": str(args.refstem)+'_chr${cchr}.hap.gz',
           "ref_leg": str(args.refstem)+'_chr${cchr}.legend.gz',
           "map": str(args.map_dir)+'/genetic_map_chr${cchr}_combined_b37.txt',
           "Ne": str(args.ne),
           "buffer": str(args.buffer),
           "out": str(outdot)+'.imp.${cname}',
           "seedtxt": str(seedtxt)
           }

uger_imp = open(str(outdot)+'.imp_chunks.sub.sh', 'w')
uger_imp.write(uger_imp_template.format(**jobdict))
uger_imp.close()

# submit
print ' '.join(['qsub',uger_imp.name]) + '\n'
# subprocess.check_call(' '.join(['qsub',uger_imp.name]), shell=True)
print 'Imputation jobs submitted for %d chunks.\n' % nchunks




print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)

# TODO:
# agg script
# - resub failed chunks in long
# - best-guess
# - re-combine chunks
# - mendel, info qc?
