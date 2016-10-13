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


import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import os
import subprocess
import argparse
from args_impute import parserbase, parsercluster, parserjob
from py_helpers import unbuffer_stdout, find_exec, file_len
from blueprint import send_job, load_job, save_job, read_clust_conf
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############
parser = argparse.ArgumentParser(prog='agg_imp.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parsercluster, parserjob])
                    
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
if args.addout is not None and str(args.addout) != '':
    out_dir = wd + '/cobg_dir_'+str(args.out)+'_'+str(args.addout)
else:
    out_dir = wd + '/cobg_dir_'+str(args.out)



#############
print '\n...Checking dependencies...'
#############

plink_ex = find_exec('plink', key='p2loc')

# get directory containing current script
# (to get absolute path for scripts)
rp_bin = os.path.dirname(os.path.realpath(__file__))
clust_conf = read_clust_conf()

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
    
    ###
    # copy original submit script
    # replace chunk list, name, number of tasks, memory spec
    # resubmit
    ###

    # load pickle of job info
    orig_job_conf = 'bg.chunks.'+str(outdot)+'.pkl'
    
    if not os.path.isfile(orig_job_conf):
        orig_job_file = str(outdot)+'.bg_chunks.sub.sh'
        raise IOError("Unable to find previous job configuration pickle %s.\
            \nRefer to previous submit script %s to modify/resubmit.\n" % (str(orig_job_conf),str(orig_job_file)))

    
    cmd_templ, job_dict, sendjob_dict = load_job(orig_job_conf)

    # rename resub
    sendjob_dict['jobname'] = 'bg.chunks.'+str(outdot)+'.resub_'+str(nummiss)
    
    sendjob_dict['logname'] = str('bg.chunks.'+str(outdot)+'.resub_'+str(nummiss)+'.'+str(clust_conf['log_task_id'])+'.sub.log')

    # increase memory and walltime
    # TODO: consider how to scale mem/time here
    oldmem = sendjob_dict['mem']
    sendjob_dict['mem'] = int(oldmem) + 4000

    oldtime = sendjob_dict['walltime']
    sendjob_dict['walltime'] = int(oldtime)*4
    
    # replace chunk file and set number of new jobs
    sendjob_dict['njobs'] = int(nummiss)

    job_dict['cfile'] = tmp_chunk_file_name
    
    
    # re-save new settings (primarily to track updating mem and walltime)
    save_job(jfile=orig_job_conf, cmd_templ=cmd_templ, job_dict=job_dict, sendjob_dict=sendjob_dict)

    
    # submit
    bg_cmd = cmd_templ.format(**job_dict)

    jobres = send_job(jobname=sendjob_dict['jobname'],
	              cmd=bg_cmd,
	              logname=sendjob_dict['logname'],
	              mem=sendjob_dict['mem'],
	              walltime=sendjob_dict['walltime'],
	              njobs=sendjob_dict['njobs'],
	              sleep=sendjob_dict['sleep'])
        
    print 'Best-guess jobs resubmitted for %d chunks.\n' % nummiss
    
    
    
    print '\n...Replacing this aggregation job in the queue...'

    os.chdir(wd)
    agg_log = 'agg_imp.'+str(outdot)+'.resub_'+str(nummiss)+'.sub.log'

    # TODO: consider queue/mem for agg
    send_job(jobname='agg.imp.'+str(outdot),
             cmd=' '.join(sys.argv[:]),
             logname=agg_log,
             mem=8000,
             walltime=30,
             wait_name='bg.chunks.'+str(outdot)+'.resub_'+str(nummiss),
             wait_num=str(jobres).strip(),
             sleep=args.sleep)

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
