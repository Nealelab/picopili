#! /usr/bin/env python

####################################
# agg_imp.py
# written by Raymond Walters, January 2016
"""
Aggregate best-guess calls and info scores
"""
# Overview:
# 1) Parse arguments
#    - check dependencies, print args, etc
# 2) Get genomic chunk definitions
# 3) Check files for all chunks present
#    - if missing, resub (fail if resubbing same)
# 4) Combine best-guess files
# 
####################################

# TODO: enable failed chunk check



import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import os
import subprocess
from args_impute import *
from py_helpers import unbuffer_stdout, read_conf, file_len #, file_tail, link, warn_format
unbuffer_stdout()
# warnings.formatwarning = warn_format

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############
parser = argparse.ArgumentParser(prog='agg_imp.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parsercluster])
                    
args, extra_args = parser.parse_known_args()



# TODO: arg print



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
bg_dir = wd + '/imp_postproc'
out_dir = wd + '/cobg_dir_'+str(args.out)+'_'+str(args.addout)







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
print '\n...Checking best-guess calling/filtering was successful...'
######################

chunks = {}
mis_chunks = {}


dumphead = chunks_in.readline()
for line in chunks_in:
    (chrom, start, end, chname) = line.split()
    chunks[str(chname)] = [str(chrom), int(start), int(end)]

    # verify output files exists
    ch_bed = bg_dir + '/' + str(outdot) + '.bg.filtered.' + str(chname) + '.bed'
    ch_bim = bg_dir + '/' + str(outdot) + '.bg.filtered.' + str(chname) + '.bim_rsids'
    ch_fam = bg_dir + '/' + str(outdot) + '.bg.filtered.' + str(chname) + '.fam_trans'
    ch_inf = bg_dir + '/' + str(outdot) + '.bg.filtered.' + str(chname) + '.info'
   
    # record failed chunks
    if not os.path.isfile(ch_bed):
        mis_chunks[str(chname)] = [str(chrom), int(start), int(end)]
    elif not os.path.isfile(ch_bim):
        mis_chunks[str(chname)] = [str(chrom), int(start), int(end)]
    elif not os.path.isfile(ch_fam):
        mis_chunks[str(chname)] = [str(chrom), int(start), int(end)]
    elif not os.path.isfile(ch_inf):
        mis_chunks[str(chname)] = [str(chrom), int(start), int(end)]

chunks_in.close()



###############
# if there are missing chunks, restart and resub agg script
###############

if len(mis_chunks) > 0:
    nummiss = len(mis_chunks)
    print 'Missing filtered best-guess for %d chunks. Preparing to resubmit...' % nummiss
    
    # check if already tried with this number of missing chunks
    os.chdir(bg_dir)
    tmp_chunk_file_name = 'tmp_missing_'+str(nummiss)+'_chunks.'+str(outdot)+'.txt'
    if os.path.isfile(tmp_chunk_file_name):
        print '\n####################'
        print 'ERROR:'
        print 'Found previous attempt to resubmit %d failed chunks.' % int(nummiss)
        print 'Imputation is likely stuck.'
        print 'See %s/%s for failing chunks.' % (bg_dir, tmp_chunk_file_name)
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
    orig_uger_file = open(str(outdot)+'.bg_chunks.sub.sh', 'r')
    new_uger_file = open(str(outdot)+'.bg_chunks.resub_'+ str(nummiss)+'_chunks.sub.sh', 'w')
    
    for line in orig_uger_file:
        if '#$ -t ' in line:
            new_uger_file.write('#$ -t 1-'+str(nummiss)+'\n')
        elif '#$ -l m_mem_free' in line:
	    new_uger_file.write('#$ -l m_mem_free=8g \n')     
        elif '#$ -q short' in line:
	    new_uger_file.write('#$ -q long \n')
        else:
            line=line.replace(chunk_file_name, tmp_chunk_file_name)
            line=line.replace('.$TASK_ID.','.tmp'+str(nummiss)+'.$TASK_ID.')
            line=line.replace('#$ -N bg.chunks.'+str(outdot), '#$ -N bg.chunks.'+str(outdot)+'.resub_'+str(nummiss))
            new_uger_file.write(line)
            
    orig_uger_file.close()
    new_uger_file.close()

    print ' '.join(['qsub',new_uger_file.name]) + '\n'
    subprocess.check_call(' '.join(['qsub',new_uger_file.name]), shell=True)
    print 'Best-guess jobs resubmitted for %d chunks.\n' % nummiss
    
    
    print '\n...Replacing this aggregation job in the queue...'

    # TODO: consider queue/mem for agg
    agg_log = 'agg_imp.'+str(outdot)+'.resub_'+str(nummiss)+'.qsub.log'
    uger_agg = ' '.join(['qsub',
                            '-hold_jid','bg.chunks.'+str(outdot)+'.resub_'+str(nummiss),
                            '-q', 'long',
                            '-l', 'm_mem_free=8g',
                            '-N', 'agg.imp.'+str(outdot),
                            '-o', agg_log,
                            str(rp_bin)+'/uger.sub.sh',
                            str(args.sleep),
                            ' '.join(sys.argv[:])])
    
    print uger_agg + '\n'
    subprocess.check_call(uger_agg, shell=True)

    print '\n############'
    print '\n'
    print 'All jobs submitted.\n'
    exit(0)




######################
print '\n...Merging imputation info scores...'
######################

# re-read from chunk file to get chunks in order
# also merge info file

os.mkdir(out_dir)
os.chdir(out_dir)
merge_list = open(str(outdot)+'.chunk_merge_list.txt', 'w')
info_file = open(str(outdot)+'.cobg.filtered.info', 'a')

chunks_in = open(str(chunk_dir) +'/'+ chunk_file_name, 'r')
dumphead = chunks_in.readline()

# track first for info header
first = True

for line in chunks_in:
    (chrom, start, end, chname) = line.split()
    
    # double check in case resub was skipped
    if str(chname) in mis_chunks.keys():
        continue

    # file names
    ch_bed = bg_dir + '/' + str(outdot) + '.bg.filtered.' + str(chname) + '.bed'
    ch_bim = bg_dir + '/' + str(outdot) + '.bg.filtered.' + str(chname) + '.bim_rsids'
    ch_fam = bg_dir + '/' + str(outdot) + '.bg.filtered.' + str(chname) + '.fam_trans'
    ch_inf = bg_dir + '/' + str(outdot) + '.bg.filtered.' + str(chname) + '.info'
    
    # put chunk in merge list
    merge_list.write(' '.join([ch_bed, ch_bim, ch_fam]) + '\n')
    
    # add info to full file
    if first:
        subprocess.check_call(['cat',ch_inf],
                              stdout = info_file)
        first = False
    else:
        subprocess.check_call(['tail','-n','+2',ch_inf],
                              stdout = info_file)

chunks_in.close()
merge_list.close()
info_file.close()



######################
print '\n...Merging best-guess genotypes...'
######################

merge_log = open(str(outdot)+'.cobg.filtered.merge.log', 'w')
subprocess.check_call([plink_ex,
                       '--merge-list', str(merge_list.name),
                       '--merge-mode',str(4),
                       '--make-bed',
                       '--out', str(outdot)+'.cobg.filtered'],
                       stderr=subprocess.STDOUT, 
                       stdout=merge_log) 
merge_log.close()



######################
print '\n...Verifying output...'
######################

assert os.path.isfile(str(outdot)+'.cobg.filtered.bed')
assert os.path.isfile(str(outdot)+'.cobg.filtered.bim')
assert os.path.isfile(str(outdot)+'.cobg.filtered.fam')
assert file_len(str(outdot)+'.cobg.filtered.bim')+1 == file_len(str(outdot)+'.cobg.filtered.info')
# TODO: here




# finish
print '\n############'
print '\n'
print 'SUCCESS!'
exit(0)

# eof