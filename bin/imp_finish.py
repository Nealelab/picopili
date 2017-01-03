#! /usr/bin/env python

####################################
# final_imp.py
# written by Raymond Walters, December 2016
"""
Verifies imputation completed successfully
"""
#
####################################


import os
import sys
from py_helpers import file_check_email, file_len
from blueprint import load_job, save_job, read_clust_conf, send_job

outdot = str(sys.argv[1])


# get cluster configuration
# needed for specifying logfile names with clust_conf['log_task_id']
clust_conf = read_clust_conf()



if not os.path.isfile(str(outdot)+'.cobg.filtered.bed'):
    filt_done = False    
elif not os.path.isfile(str(outdot)+'.cobg.filtered.bim'):
    filt_done = False
elif not os.path.isfile(str(outdot)+'.cobg.filtered.fam'):
    filt_done = False
elif not file_len(str(outdot)+'.cobg.filtered.bim')+1 == file_len(str(outdot)+'.cobg.filtered.info'):
    filt_done = False
else:
    filt_done = True

if not os.path.isfile(str(outdot)+'.cobg.bed'):
    bg_done = False    
elif not os.path.isfile(str(outdot)+'.cobg.bim'):
    bg_done = False
elif not os.path.isfile(str(outdot)+'.cobg.fam'):
    bg_done = False
elif not file_len(str(outdot)+'.cobg.bim')+1 == file_len(str(outdot)+'.cobg.info'):
    bg_done = False
else:
    bg_done = True

if not bg_done:
    # resub
    j1_conf = 'merge.bg.'+str(outdot)+'.pkl'
    cmd_templ, job_dict, sendjob_dict = load_job(j1_conf)

    # rename resub
    sendjob_dict['jobname'] = 'merge.bg.'+str(outdot)+'.resub'
    sendjob_dict['logname'] = str('merge.bg.'+str(outdot)+'.resub.'+str(clust_conf['log_task_id'])+'.sub.log')

    # increase memory and walltime
    # TODO: consider how to scale mem/time here
    oldmem = sendjob_dict['mem']
    sendjob_dict['mem'] = int(oldmem)*2

    oldtime = sendjob_dict['walltime']
    sendjob_dict['walltime'] = int(oldtime)*2
    
    # re-save new settings (primarily to track updating mem and walltime)
    save_job(jfile=j1_conf, cmd_templ=cmd_templ, job_dict=job_dict, sendjob_dict=sendjob_dict)

    # submit
    merge_cmd1 = cmd_templ.format(**job_dict)

    jobres = send_job(jobname=sendjob_dict['jobname'],
	              cmd=merge_cmd1,
	              logname=sendjob_dict['logname'],
	              mem=sendjob_dict['mem'],
	              walltime=sendjob_dict['walltime'],
	              sleep=sendjob_dict['sleep'])    


if not filt_done:
    # resub
    j2_conf = 'merge.bg_filt.'+str(outdot)+'.pkl'
    cmd_templ, job_dict, job2_dict = load_job(j2_conf)

    # rename resub
    job2_dict['jobname'] = 'merge.bg_filt.'+str(outdot)+'.resub'
    job2_dict['logname'] = str('merge.bg_filt.'+str(outdot)+'.resub.'+str(clust_conf['log_task_id'])+'.sub.log')

    # increase memory and walltime
    # TODO: consider how to scale mem/time here
    oldmem = sendjob_dict['mem']
    sendjob_dict['mem'] = int(oldmem)*2

    oldtime = sendjob_dict['walltime']
    sendjob_dict['walltime'] = int(oldtime)*2
    
    # re-save new settings (primarily to track updating mem and walltime)
    save_job(jfile=j2_conf, cmd_templ=cmd_templ, job_dict=job_dict, sendjob_dict=job2_dict)

    # submit
    merge_cmd2 = cmd_templ.format(**job_dict)

    jobres2 = send_job(jobname=job2_dict['jobname'],
	              cmd=merge_cmd2,
	              logname=job2_dict['logname'],
	              mem=job2_dict['mem'],
	              walltime=job2_dict['walltime'],
	              sleep=job2_dict['sleep'])



if bg_done and filt_done:
    file_check_email(str(outdot)+'.cobg.filtered.bed', 'merge.bg_filt.'+str(outdot))

    # TODO: add cleanup here

    # finish
    print '\n############'
    print '\n'
    print 'SUCCESS!'
    return(0)

else:
    # resub this
    fin_conf = 'merge.bg_filt.'+str(outdot)+'.pkl'
    final_call, job_dict, final_dict = load_job(fin_conf)
    
    final_dict['logname'] = 'imp.check_fin.'+str(outdot)+'.resub.log'
    
    # which to wait on
    if not bg_done:
        final_dict['wait_name'] = sendjob_dict['jobname']
        final_dict['wait_num'] = str(jobres).strip()
    else:
        final_dict['wait_name'] = job2_dict['jobname']
        final_dict['wait_num'] = str(jobres2).strip()

    
    send_job(jobname=final_dict['jobname'],
         cmd=str(final_call),
         logname=final_dict['logname'],
         mem=final_dict['mem'],
         walltime=final_dict['walltime'],
         wait_name=final_dict['wait_name'],
         wait_num=final_dict['wait_num'],
         sleep=final_dict['sleep'])
    
    return(1)
    