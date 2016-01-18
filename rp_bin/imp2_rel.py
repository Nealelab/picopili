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

import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import os
import subprocess
from args_impute import *
from py_helpers import unbuffer_stdout, file_len, link, read_conf
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############
parser = argparse.ArgumentParser(prog='imp2_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parserimpute, parserref, parserchunk, parsercluster])
                    
args, extra_args = parser.parse_known_args()

# TODO: allow options
args.refstem = '/humgen/atgu1/fs03/shared_resources/1kG/shapeit/1000GP_Phase3'



# get useful modified args
if args.addout is not None and str(args.addout) != '':
    addout_txt = ['--addout',str(args.addout)]
    outdot = str(args.out)+'.'+str(args.addout)
else:
    addout_txt = ['','']
    outdot = str(args.out)
    
if args.imp_seed is not None and str(args.imp_seed) != '' and int(args.imp_seed) > 0:
    seedtxt = '-seed '+str(args.imp_seed)
else:
    seedtxt = ''


# report settings in use
print '\nBasic settings:'
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
if args.addout is not None:
    print '--addout '+str(args.addout)

print '\nIMPUTE2 arguments:'
print '--Ne '+str(args.Ne)
print '--buffer '+str(args.buffer)
if args.imp_seed is not None and str(args.imp_seed) != '' and int(args.imp_seed) > 0:
    print '--seed '+str(args.imp_seed)

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
shapeit_ex = configs['shloc'] + '/bin/shapeit'

# get directory containing current script
# (to get absolute path for scripts)
rp_bin = os.path.dirname(os.path.realpath(__file__))
chunker_ex = rp_bin+'/chunk_snps.py'



# directories
wd = os.getcwd()
shape_dir = wd + '/phase_chr'





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
print '\n...Verifying pre-phasing was successful...'
######################

bad_chr = []

for chrom in xrange(1,22):
    haps_out = str(shape_dir)+'/'+str(outdot)+'.chr'+str(chrom)+'.phased.haps'
    samp_out = str(shape_dir)+'/'+str(outdot)+'.chr'+str(chrom)+'.phased.sample'


    if not os.path.isfile(haps_out):
        bad_chr.append(chrom)
    elif not os.path.isfile(samp_out):
        bad_chr.append(chrom)
        

# TODO: resub shapeit if failed
# TODO: re-queue this job
if bad_chr:    
    num_chr = len(bad_chr)
    print 'Missing pre-phasing results for %d chromosomes. Preparing to resubmit...' % num_chr

    os.chdir(shape_dir)
        
    # verify haven't already tried this resub
    uger_phase_name = str(outdot)+'.shape.resub_'+str(num_chr)+'_chr.sub.sh'
    if os.path.isfile(uger_phase_name):
        print '\n####################'
        print 'ERROR:'
        print 'Found previous attempt to resubmit %d failed chromosomes.' % int(num_chr)
        print 'Pre-phasing is likely stuck.'
        print 'Problem chromosomes: %s' % (','.join(bad_chr))
        print 'Exiting...\n'
        exit(1)

    # make submit script
    # using this structure to get adaptive chromosome list
    uger_phase_template = """#!/usr/bin/env sh
    #$ -j y
    #$ -cwd
    #$ -V
    #$ -N {jname}
    #$ -q long
    #$ -l m_mem_free={mem}g
    #$ -pe smp {threads}
    #$ -t 1-{nchr}
    #$ -o {outlog}
    
    source /broad/software/scripts/useuse
    reuse -q Anaconda
    sleep {sleep}
    
    chrs=({chr_list})
    chrom=${{chrs[${{SGE_TASK_ID}}-1]}}

    {shape_ex} {bed} {map} {ref} {window} {thread_str} {seed_str} {outmax} {shapelog} 
    
    # eof
    """

#    shape_call = [shapeit_ex,
#                  '--input-bed', chrstem+'.bed', chrstem+'.bim', chrstem+'.fam',
#                  '--input-map', args.map_dir+'/genetic_map_chr\$tasknum_combined_b37.txt',
#                  '--input-ref', refstem+'_chr\$tasknum.hap.gz', refstem+'_chr\$tasknum.legend.gz', refstem+'.sample',
#                  '--window', str(args.window),
#                  '--duohmm',
#                  '--thread', str(args.threads),
#                  '--seed', str(args.shape_seed),
#                  '--output-max', outstem+'.phased.haps', outstem+'.phased.sample',
#                  '--output-log', outstem+'.shape.log']

    
    # fill in template
    chrstem = str(args.bfile)+'.hg19.ch.fl.chr${chrom}'
    outstem = str(outdot)+'.chr${chrom}'    
    jobdict = {"jname": 'shape.'+str(outdot)+'.resub_'+str(num_chr),
               "mem": str(args.mem_req),
               "threads": str(args.threads),
               "nchr": str(num_chr),
               "outlog": 'shape.'+str(outdot)+'.resub_'+str(num_chr)+'.qsub.$TASK_ID.log',
               "sleep": str(args.sleep),
               "chr_list": ' '.join(bad_chr),
               "shape_ex": str(shape_ex),
               "bed": '--input-bed '+str(chrstem)+'.bed '+str(chrstem)+'.bim '+str(chrstem)+'.fam',
               "map": '--input-map '+str(args.map_dir)+'/genetic_map_chr${chrom}_combined_b37.txt',
               "ref": '--input-ref '+str(args.refstem)+'_chr${chrom}.hap.gz '+str(args.refstem)+'_chr${chrom}.legend.gz '+str(args.refstem)+'.sample',
               "window": '--window '+str(args.window),
               "thread_str": '--thread '+str(args.threads),
               "seed_str": '--seed '+str(args.shape_seed),
               "outmax": '--output-max '+str(outstem)+'.phased.haps '+str(outstem)+'.phased.sample',
               "shapelog": str(outstem)+'.shape.resub_'+str(num_chr)+'.log',
               }
    
    uger_phase = open(uger_phase_name, 'w')
    uger_phase.write(uger_phase_template.format(**jobdict))
    uger_phase.close()
    
    # submit
    print ' '.join(['qsub',uger_phase_name]) + '\n'
    subprocess.check_call(' '.join(['qsub',uger_phase_name]), shell=True)
    print 'Pre-phasing jobs re-submitted for %d chromosomes.\n' % num_chr



    # put this job back in the queue
    print '\n...Replacing this imputation job in the queue...'
    
    os.chdir(wd)
    imp_log = 'imp_chunks.'+str(outdot)+'.qsub.log'
    uger_imp = ' '.join(['qsub',
                            '-hold_jid','imp.chunks.'+str(outdot),
                            '-q', 'short',
                            '-l', 'm_mem_free=8g',
                            '-N', 'imp.chunks.'+str(outdot),
                            '-o', imp_log,
                            str(rp_bin)+'/uger.sub.sh',
                            str(args.sleep),
                            ' '.join(sys.argv[:])])
    
    print uger_imp + '\n'
    subprocess.check_call(uger_imp, shell=True)

    print '\n############'
    print '\n'
    print 'All jobs submitted.\n'
    exit(0)




######################
print '\n...Setting up working directories...'
######################

# working directories
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
           "Ne": str(args.Ne),
           "buffer": str(args.buffer),
           "out": str(outdot)+'.imp.${cname}',
           "seedtxt": str(seedtxt)
           }

uger_imp = open(str(outdot)+'.imp_chunks.sub.sh', 'w')
uger_imp.write(uger_imp_template.format(**jobdict))
uger_imp.close()

# submit
print ' '.join(['qsub',uger_imp.name]) + '\n'
subprocess.check_call(' '.join(['qsub',uger_imp.name]), shell=True)
print 'Imputation jobs submitted for %d chunks.\n' % nchunks




###
# submit next imputation task
###
if args.full_pipe:
    ######################
    print '\n...Queuing best-guess script...'
    ######################
    
    os.chdir(wd)
    next_call = str(rp_bin) + '/bg_imp.py '+' '.join(sys.argv[1:])

    # TODO: consider queue/mem for agg
    bg_log = 'bg_imp.'+str(outdot)+'.qsub.log'
    uger_bg = ' '.join(['qsub',
                            '-hold_jid','imp.chunks.'+str(outdot),
                            '-q', 'short',
                            '-l', 'm_mem_free=4g',
                            '-N', 'bg.chunks.'+str(outdot),
                            '-o', bg_log,
                            str(rp_bin)+'/uger.sub.sh',
                            str(args.sleep),
                            next_call])
    
    print uger_bg + '\n'
    subprocess.check_call(uger_bg, shell=True)






print '\n############'
print '\n'
print 'SUCCESS!\n'
print 'All jobs submitted.\n'
exit(0)
# eof