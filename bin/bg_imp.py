#! /usr/bin/env python

####################################
# bg_imp.py
# written by Raymond Walters, January 2016
"""
Generate best-guess calls from IMPUTE output
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
# 6) Fix bim files (extract rsids, etc)
#    - farm 5 & 6 out to cluster
# 
####################################

# TODO: turn on subprocess submissions
# enable missing chunk check


import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import os
import warnings
import argparse
from warnings import warn
from textwrap import dedent
from args_impute import parserbase, parserbg, parsercluster, parserjob
from py_helpers import unbuffer_stdout, find_exec, file_tail, link, warn_format, read_conf
from blueprint import send_job, init_sendjob_dict, save_job, load_job, read_clust_conf
unbuffer_stdout()
warnings.formatwarning = warn_format

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############
parser = argparse.ArgumentParser(prog='bg_imp.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parserbg, parsercluster, parserjob])
                   
args, extra_args = parser.parse_known_args()


# process call threshold
if args.hard_call_th is None:
    hard_call_th = 1.0 - float(args.bg_th)
else:
    hard_call_th = float(args.hard_call_th)
    if args.bg_th > 1.0-hard_call_th:
        if args.bg_th == 0.8:
            warn("Default value of --bg-th (0.8) overridden by less strict value for --hard-call-th (%f)." % hard_call_th)
        else:
            hard_call_th = 1.0 - float(args.bg_th)
            warn("Both --hard-call-th and --bg_th specified. Using stricter value (==> hard-call-th %f)" % hard_call_th)
    elif args.bg_th < 1.0-hard_call_th:
        if args.bg_th != 0.8:
	    warn("Both --hard-call-th and --bg_th specified. Using stricter value (==> hard-call-th %f)" % hard_call_th)


assert float(hard_call_th) > 0
assert float(hard_call_th) < 1



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


######
# setup text for plink filtering
######

# mendel errors
if args.keep_mendel:
    mendel_txt = ''
elif args.mendel == 'none':
    args.keep_mendel = True
    mendel_txt = ''
elif args.mendel == 'multigen':
    mendel_txt = '--set-me-missing --mendel-multigen'
elif args.mendel == 'duos':
    mendel_txt = '--set-me-missing --mendel-duos'
elif args.mendel == 'trios':
    mendel_txt = '--set-me-missing'
else:
    raise ValueError("Invalid argument for --mendel. Something has failed.")

# info scores
# based on impute2 output format for _info files
if args.info_th is None and args.max_info_th is None:
    info_txt = ''
else:
    # init, then add thresholds
    if args.imp_version==2:
        info_txt = '--qual-scores '+str(imp_dir)+'/'+str(outdot)+'.imp.${{cname}}_info' +' 5 2 1'
    elif args.imp_version==4:
        info_txt = '--qual-scores '+str(imp_dir)+'/'+str(outdot)+'.imp.${{cname}}.qctool_info.txt' +' 17 2 \'#\''
    # minimum info
    if args.info_th >= 0.0 and args.info_th <= 1.0:
        info_txt = info_txt + ' --qual-threshold '+str(args.info_th)
    elif args.info_th is not None:
        raise ValueError("Invalid argument for --info_th. Must be between 0 and 1.")
    # maximum info
    if args.max_info_th >= 0.0:
        info_txt = info_txt + ' --qual-max-threshold '+str(args.max_info_th)
    elif args.max_info_th is not None:
        raise ValueError("Invalid argument for --max_info_th. Must be between 0 and 1.")    

# maf
if args.maf_th is None:
   maf_txt = '' 
elif args.maf_th >= 0.0 and args.maf_th <= 0.5:
    maf_txt = '--maf '+str(args.maf_th)
else:
    raise ValueError("Invalid argument for --maf_th. Must be between 0 and 0.5.")

# mac
if args.mac_th is None:
    mac_txt = ''
elif int(args.mac_th) > 0:
    mac_txt = '--mac '+str(int(args.mac_th))
else:
    raise ValueError("Invalid argument for --mac_th. Must be greater than 0.")

# call rate
if args.miss_th is None:
    geno_txt = ''
elif args.miss_th > 0.0 and args.miss_th < 1.0:
    geno_txt = '--geno '+str(args.miss_th)
else:
    raise ValueError("Invalid argument for --miss_th. Must be between 0 and 1.")



# report settings in use
print '\nBasic settings:'
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
if args.addout is not None:
    print '--addout '+str(args.addout)

print '\nBest-guess genotype calling:'
if args.hard_call_th is None:
    print '--bg-th '+str(args.bg_th)
else:
    print '--hard-call-th '+str(hard_call_th)
print '--info-th '+str(args.info_th)
print '--max-info-th '+str(args.max_info_th)
if args.keep_mendel:
    print '--keep-mendel'
else:
    print '--mendel '+str(args.mendel)
print '--maf-th '+str(args.maf_th)
if args.mac_th is not None:
    print '--mac-th '+str(args.mac_th)
print '--miss-th '+str(args.miss_th)

print '\nCluster settings:'
print '--sleep '+str(args.sleep)
if args.full_pipe:
    print '--full-pipe'


#############
print '\n...Checking dependencies...'
#############

plink_ex = find_exec('plink',key='p2loc')

# get directory containing current script
# (to get absolute path for scripts)
rp_bin = os.path.dirname(os.path.realpath(__file__))
rs_ex = str(rp_bin)+'/rs_trans.py'

# get cluster configuration
# needed for specifying logfile names with clust_conf['log_task_id']
conf_file = os.environ['HOME']+"/picopili.conf"
configs = read_conf(conf_file)
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
print '\n...Checking imputation was successful...'
######################

chunks = {}
mis_chunks = {}

dumphead = chunks_in.readline()
for line in chunks_in:
    (chrom, start, end, chname) = line.split()
    chunks[str(chname)] = [str(chrom), int(start), int(end)]

    # verify output file exists
    if args.imp_version==2:
        ch_imp = imp_dir + '/' + str(outdot) + '.imp.' + str(chname) + '.gz'
    elif args.imp_version==4:
        ch_imp = imp_dir + '/' + str(outdot) + '.imp.' + str(chname) + '.gen.gz'
    
    # verify completed successfully
    # make sure use most recent log for resubs
    if args.imp_version==2:
        ch_sum = imp_dir + '/' + str(outdot) + '.imp.' + str(chname) + '_summary'
        fin_string = '[0.9-1.0]'
    elif args.imp_version==4:
        ch_sum = imp_dir + '/' + str(outdot) + '.imp.' + str(chname) + '.qctool_info.log'
        fin_string = 'Thank you for using qctool'
    
    # record failed chunks
    if not os.path.isfile(ch_imp):
        mis_chunks[str(chname)] = [str(chrom), int(start), int(end)]
    elif fin_string not in file_tail(ch_sum, n=1):
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
    
    ###
    # copy original submit script
    # replace chunk list, name, number of tasks, memory spec
    # resubmit
    ###

    # load pickle of job info
    orig_job_conf = 'imp.chunks.'+str(outdot)+'.pkl'
    
    if not os.path.isfile(orig_job_conf):
        orig_job_file = str(outdot)+'.imp_chunks.sub.sh'
        raise IOError("Unable to find previous job configuration pickle %s.\
            \nRefer to previous submit script %s to modify/resubmit.\n" % (str(orig_job_conf),str(orig_job_file)))

    
    cmd_templ, job_dict, sendjob_dict = load_job(orig_job_conf)

    # rename resub
    sendjob_dict['jobname'] = 'imp.chunks.'+str(outdot)+'.resub_'+str(nummiss)
    
    sendjob_dict['logname'] = str('imp.chunks.'+str(outdot)+'.resub_'+str(nummiss)+'.'+str(clust_conf['log_task_id'])+'.sub.log')

    # increase memory and walltime
    # TODO: consider how to scale mem/time here
    oldmem = sendjob_dict['mem']
    sendjob_dict['mem'] = int(oldmem)*2

    oldtime = sendjob_dict['walltime']
    sendjob_dict['walltime'] = int(oldtime)*4
    
    # replace chunk file and set number of new jobs
    sendjob_dict['njobs'] = int(nummiss)

    job_dict['cfile'] = tmp_chunk_file_name
    
    
    # re-save new settings (primarily to track updating mem and walltime)
    save_job(jfile=orig_job_conf, cmd_templ=cmd_templ, job_dict=job_dict, sendjob_dict=sendjob_dict)

    
    # submit
    imp_cmd = cmd_templ.format(**job_dict)

    jobres = send_job(jobname=sendjob_dict['jobname'],
             	      cmd=imp_cmd,
             	      logname=sendjob_dict['logname'],
	              mem=sendjob_dict['mem'],
	              walltime=sendjob_dict['walltime'],
	              njobs=sendjob_dict['njobs'],
	              sleep=sendjob_dict['sleep'],
		      forcearray=True)
        
    print 'GWAS jobs resubmitted for %d chunks.\n' % nummiss

    
    
    print '\n...Replacing this best-guess job in the queue...'

    os.chdir(wd)
    bg_log = 'bg.'+str(outdot)+'.resub_'+str(nummiss)+'.sub.log'

    # TODO: consider queue/mem for agg
    send_job(jobname='bg.chunks.'+str(outdot),
             cmd=' '.join(sys.argv[:]),
             logname=bg_log,
             mem=8000,
             walltime=2, # week
             wait_name='imp.chunks.'+str(outdot)+'.resub_'+str(nummiss),
	     wait_num=str(jobres).strip(),
             sleep=args.sleep)

    print '\n############'
    print '\n'
    print 'All jobs submitted.\n'
    exit(0)



######################
print '\n...Setting up working directory...'
######################

# working directories
proc_dir = wd + '/imp_postproc'
os.mkdir(proc_dir)
print 'Impute output processing: %s' % proc_dir
os.chdir(proc_dir)
link(str(chunk_dir)+'/'+str(outdot)+'.chunks.txt', str(outdot)+'.chunks.txt', 'genomic chunk results')


######################
print '\n...Generating best-guess genotypes...'
######################

# best-guess job script for each chunk
bg_templ = dedent("""\
    cname=`awk -v a={task} 'NR==a+1{cbopen}print $4{cbclose}' {cfile}`
    cchr=`awk -v a={task} 'NR==a+1{cbopen}print $1{cbclose}' {cfile}`
    
    {plink_ex} --gen {gen_in} --sample {samp_in} --oxford-single-chr ${cbopen}cchr{cbclose} --oxford-pheno-name plink_pheno --hard-call-threshold {hard_call_th} --missing-code -9,NA,na --allow-no-sex --silent --memory 4000 --out {out_str} 
    
    sleep {sleep}
    # note: Mendel errors checked after --update-parents, see https://www.cog-genomics.org/plink2/order
    {plink_ex} --bfile {out_str} {mendel_txt} --pheno {idnum} --mpheno 4 --update-parents {idnum} --allow-no-sex --make-bed --silent --memory 2000 --out {out_str2}
    rm {out_str}.bed
    rm {out_str}.bim
    rm {out_str}.fam
    
    sleep {sleep}
    {plink_ex} --bfile {out_str2} {maf_txt} {mac_txt} {geno_txt} {info_txt} --allow-no-sex --make-bed --silent --memory 2000 --out {out_str_filt}
    rm {out_str2}.bed
    rm {out_str2}.bim
    rm {out_str2}.fam
    
    {rs_ex} --chunk ${cbopen}cname{cbclose} --name {outdot} --imp-dir {imp_dir} --fam-trans {trans}
""")

# get number of chunks
nchunks = len(chunks)

# info to fill in job template

if args.imp_version==2:
    gen_in = str(imp_dir)+'/'+str(outdot)+'.imp.${{cname}}.gz'
elif args.imp_version==4:
    gen_in = str(imp_dir)+'/'+str(outdot)+'.imp.${{cname}}.gen.gz'


jobdict = {"task": "{task}",
           "sleep": str(args.sleep),
           "cfile": str(outdot)+'.chunks.txt',
           "plink_ex": str(plink_ex),
           "gen_in": str(gen_in),
           "samp_in": str(shape_dir)+'/'+str(outdot)+'.chr${{cchr}}.phased.sample',
           "hard_call_th": str(hard_call_th),
           "out_str": str(outdot)+'.bg.${{cname}}',
           "mendel_txt": str(mendel_txt),
           "info_txt": str(info_txt),
           "out_str2": str(outdot)+'.bg.tmp.${{cname}}',
           "maf_txt": str(maf_txt),
           "mac_txt": str(mac_txt),
           "geno_txt": str(geno_txt),
           "out_str_filt": str(outdot)+'.bg.filtered.${{cname}}',
           "rs_ex": str(rs_ex),
           "outdot": str(outdot),
           "imp_dir": str(imp_dir),
           "idnum": str(shape_dir)+'/'+str(args.bfile)+'.hg19.ch.fl.fam',
           "trans": str(shape_dir)+'/'+str(args.bfile)+'.hg19.ch.fl.fam.transl',
	   "cbopen":'{{',
	   "cbclose":'}}',
           }


# store job information for possible resubs
job_store_file = 'bg.chunks.'+str(outdot)+'.pkl'

clust_dict = init_sendjob_dict()
clust_dict['jobname'] = 'bg.chunks.'+str(outdot)
clust_dict['logname'] = str('bg.chunks.'+str(outdot)+'.'+str(clust_conf['log_task_id'])+'.sub.log')
clust_dict['mem'] = 8000
clust_dict['walltime'] = 2
clust_dict['njobs'] = int(nchunks)
clust_dict['sleep'] = args.sleep

save_job(jfile=job_store_file, cmd_templ=bg_templ, job_dict=jobdict, sendjob_dict=clust_dict)


# submit
# TODO: flex queue/mem reqs
bg_cmd = bg_templ.format(**jobdict)

jobres2 = send_job(jobname='bg.chunks.'+str(outdot),
	           cmd=bg_cmd,
	           logname=str('bg.chunks.'+str(outdot)+'.'+str(clust_conf['log_task_id'])+'.sub.log'),
	           mem=8000,
	           walltime=2,
	           njobs=int(nchunks),
	           sleep=args.sleep)

print 'Best-guess jobs submitted for %d chunks.\n' % nchunks


###
# submit next imputation task
###
if args.full_pipe:
    ######################
    print '\n...Queuing chunk aggregation script...'
    ######################
    
    os.chdir(wd)
    next_call = str(rp_bin) + '/agg_imp.py '+' '.join(sys.argv[1:])

    agg_log = 'agg_imp.'+str(outdot)+'.sub.log'

    # some dynamic adjustment of mem based on sample size population
    # (empirically, seem to get ~2x sites from afr vs eur)
    fam_n = file_len(str(shape_dir)+'/'+str(args.bfile)+'.hg19.ch.fl.fam')
    if fam_n > 3000:
        agg_mem = 32000
    elif fam_n > 1000:
        agg_mem = 16000
    else:
        agg_mem = 8000

    if "afr" in sys.argv[1:]:
        agg_mem = 2*agg_mem

    # TODO: consider queue/mem for agg
    send_job(jobname='agg.imp.'+str(outdot),
             cmd=next_call,
             logname=agg_log,
             mem=int(agg_mem),
             walltime=30,
             wait_name='bg.chunks.'+str(outdot),
	     wait_num=str(jobres2).strip(),
             sleep=args.sleep)

# finish
print '\n############'
print '\n'
print 'SUCCESS!'
print 'All jobs submitted.\n'
exit(0)


# eof
