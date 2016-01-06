#! /usr/bin/env python

####################################
# agg_imp.py
# written by Raymond Walters, January 2016
"""
Process IMPUTE2 output for GWAS data with related individuals
"""
# Overview:
# 1) Parse arguments
#    - check dependencies, print args, etc
# 2) Get genomic chunk definitions
# 3) Check all chunks imputed successfully
#    - if not, resubmit (tag resub with file)
#    - if previously resubbed, fail
# 4) Setup working directories for remaining tasks 
# 5) Convert chunks to plink files
#    - farm 5 & 6 out to cluster
# 6) Fix bim files (extract rsids, etc)
# 67) Merge and filter
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
from py_helpers import unbuffer_stdout, read_conf, file_tail #, file_len, link, read_conf
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############
parser = argparse.ArgumentParser(prog='agg_imp.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parsercluster])
                    
args = parser.parse_args()


# get useful modified args
if args.addout is not None and str(args.addout) != '':
    addout_txt = ['--addout',str(args.addout)]
    outdot = str(args.out)+'.'+str(args.addout)
else:
    addout_txt = ['','']
    outdot = str(args.out)

# directories
wd = os.getcwd()
shape_dir = wd + '/phase_chr'
chunk_dir = wd + '/chunks_for_imp'
imp_dir = wd + '/imp_sub'


# report settings in use
print '\nBasic settings:'
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
if args.addout is not None:
    print '--addout '+str(args.addout)

# TODO: Rest of here




#############
print '\n...Reading ricopili config file...'
#############

### read plink loc from config
conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

plink_ex = configs['p2loc']+"plink"

# get directory containing current script
# (to get absolute path for scripts)
rp_bin = os.path.dirname(os.path.realpath(__file__))


#############
print '\n...Checking dependencies...'
#############



# TODO: here




print '\n'
print '############'
print 'Begin!'
print '############'

######################
print '\n...Loading genomic chunks definitions...'
######################

chunk_file_name = str(outdot)+'.chunks.txt'
chunks_in = open(str(chunk_dir) +'/'+ chunk_file_name, 'r')


######################
print '\n...Checking imputation was successful...'
######################

chunks = {}
mis_chunks = {}

dumphead = chunks_in.readline()
for line in chunks_in:
    (chrom, start, end, chname) = line.split()
    chunks[str(chname)] = [str(chrom), int(start), int(end)]

    # verify output file exists
    ch_imp = imp_dir + '/' + str(outdot) + '.imp.' + str(chname) + '.gz'
    
    # verify completed successfully
    # - based on expected output of concordance table on last line
    ch_sum = imp_dir + '/' + str(outdot) + '.imp.' + str(chname) + '_summary'
    
    # record failed chunks
    if not os.path.isfile(ch_imp):
        mis_chunks[str(chname)] = [str(chrom), int(start), int(end)]
    elif '[0.9-1.0]' not in file_tail(ch_sum, n=1):
        mis_chunks[str(chname)] = [str(chrom), int(start), int(end)]

chunks_in.close()



###############
# if there are missing chunks, restart their imputation and resub agg script
###############
if len(mis_chunks) > 0:
    nummiss = len(mis_chunks)
    print 'Missing results for %d imputation chunks. Preparing to resubmit...' % nummiss
    
    # check if already tried with this number of missing chunks
    os.chdir(imp_dir)
    tmp_chunk_file_name = 'tmp_missing_'+str(nummiss)+'_chunks.'+str(outdot)+'.txt'
    if os.path.isfile(tmp_chunk_file_name):
        print '\n####################'
        print 'ERROR:'
        print 'Found previous attempt to resubmit %d failed chunks.' % int(nummiss)
        print 'Imputation is likely stuck.'
        print 'See %s/%s for failing chunks.' % (imp_dir, tmp_chunk_file_name)
        print 'Exiting...\n'
        exit(1)
    
    
    # else write file of just missing chunks for task array
    tmp_chunk_file = open(tmp_chunk_file_name, 'w')
    tmp_chunk_file.write(' '.join(['CHR','START','END','NAME']) + '\n')

    for ch in mis_chunks.keys():
        tmp_chunk_file.write(' '.join([str(mis_chunks[ch][0]), str(mis_chunks[ch][1]), str(mis_chunks[ch][2]), str(ch)]) + '\n')    
    tmp_chunk_file.close()
    
    print 'List of missing chunks: %s' % tmp_chunk_file.name
    
    # copy original submit script
    # replace chunk list, name, number of tasks
    orig_uger_file = open(str(outdot)+'.imp_chunks.sub.sh', 'r')
    new_uger_file = open(str(outdot)+'.imp_chunks.resub_'+ str(nummiss)+'_chunks.sub.sh', 'w')
    
    for line in orig_uger_file:
        if '#$ -t ' in line:
            new_uger_file.write('#$ -t 1-'+str(nummiss)+'\n')
        elif '#$ -l m_mem_free' in line:
	    new_uger_file.write('#$ -l m_mem_free=12g \n')     
        elif '#$ -q short' in line:
	    new_uger_file.write('#$ -q long \n')
        else:
            line=line.replace(chunk_file_name, tmp_chunk_file_name)
            line=line.replace('.$TASK_ID.','.tmp'+str(nummiss)+'.$TASK_ID.')
            line=line.replace('#$ -N imp.chunks.'+str(outdot), '#$ -N imp.chunks.'+str(outdot)+'.resub_'+str(nummiss))
            new_uger_file.write(line)
            
    orig_uger_file.close()
    new_uger_file.close()

    print ' '.join(['qsub',new_uger_file.name]) + '\n'
#    subprocess.check_call(' '.join(['qsub',new_uger_file.name]), shell=True)
    print 'GWAS jobs resubmitted for %d chunks.\n' % nummiss
    
    
    print '\n...Replacing this agg job in the queue...'

# TODO: consider queue/mem for agg
    agg_log = 'agg.'+str(outdot)+'.resub_'+str(nummiss)+'.qsub.log'
    uger_agg = ' '.join(['qsub',
                            '-hold_jid','imp.chunks.'+str(outdot)+'.resub_'+str(nummiss),
                            '-q', 'long',
                            '-l', 'm_mem_free=4g',
                            '-N', 'agg_'+str(outdot),
                            '-o', agg_log,
                            str(rp_bin)+'/uger.sub.sh',
                            str(args.sleep),
                            ' '.join(sys.argv[:])])
    
    print uger_agg + '\n'
#    subprocess.check_call(uger_agg, shell=True)

    print '\n############'
    print '\n'
    print 'All jobs submitted.\n'
    exit(0)



######################
print '\n...Setting up working directories...'
######################

# working directories
proc_dir = wd + '/imp_postproc'
bg_dir = wd + '/imp_results'
os.mkdir(proc_dir)
os.mkdir(bg_dir)
print 'Impute2 output processing: %s' % proc_dir
print 'Final imputation results: %s' % bg_dir


######################
print '\n...Generating best-guess genotypes...'
######################

# ~/plink2/plink 
# --gen tester.addertxt.imp.chr1_117_123.gz.gz 
# --sample tester.addertxt.chr1.phased.sample 
# --oxford-single-chr 1 
# --oxford-pheno-name plink_pheno 
# --hard-call-threshold .2 
# --missing-code '-9','NA','na' 
# --out test_parse_imp2

# format bim files

# write a stand-alone task for this to submit, then check here for completion



######################
print '\n...Merging and filtering imputation results..'
######################

# info
# mendel errors
# call rate
# maf 




# finish
print '\n############'
print '\n'
print 'SUCCESS!'
print 'All jobs submitted.\n'
exit(0)

# eof