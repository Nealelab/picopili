#! /usr/bin/env python

####################################
# blueprint.py
# written by Raymond Walters, September 2016
"""
manages job submission on different cluster architectures
"""
#
####################################

import os
import subprocess
from textwrap import dedent
from py_helpers import read_conf, file_len

def send_job(jobname,
             arrayfile=None,
             cmd=None,
             logname=None,
             logloc=None,
             mem=None,
             walltime=None,
#             week=None,
             njobs=None,
             maxpar=10000,
             threads=None,
             wait_file=None,
             wait_name=None,
             cluster=None,
             sleep=30,
             testonly=False):
    
    # validate args
    if arrayfile is None and cmd is None:
        raise ValueError("Require either array file or command.")
    
    elif arrayfile is not None and cmd is not None:
        raise ValueError("Require either array file or command, not both.")


    if logloc is None:
        logloc = os.getcwd()
    
    if not os.path.isdir(logloc):
        os.mkdir(logloc)
        
    if maxpar < 1:
        maxpar = 10000

    # get cluster queue name
    if cluster is None:
        conf_file = os.environ['HOME']+"/picopili.conf"
        configs = read_conf(conf_file)
        cluster = configs['cluster']

    # get queue template
    pico_bin = os.path.dirname(os.path.realpath(__file__))
    clust_dir = os.path.dirname(pico_bin) + '/cluster_templates'
    
    assert os.path.isdir(clust_dir), "Unable to find cluster job submission template directory %s" % str(clust_dir)

    # load queue configuration info
    # - submission syntax, queue names, job holds
    clust_conf = read_conf(str(clust_dir)+'/'+str(cluster)+'.conf')

    # basic template
    with open(str(clust_dir)+'/'+str(cluster)+'.sub.sh','r') as single_templ:
        templ = single_templ.read()

    # setup memory args
    if mem is None:
        mem = 2000
    mem_mb = str(int(mem))
    if int(mem) > 1000:
        mem_gb = str(int(mem)/1000)
    else:
        mem_gb = str(1)

    # multithreading arguments
    if threads is None:
        threads = 1

    # queue picking from job length
    if walltime is None:
        walltime = 1
        queue_name = clust_conf['hour_q']
    elif walltime <= 1.0:
        queue_name = clust_conf['hour_q']
    elif walltime <= 2.0:
        queue_name = clust_conf['hour2_q']
    elif walltime <= 4.0:
        queue_name = clust_conf['hour4_q']
    elif walltime <= 24.0:
        queue_name = clust_conf['day_q']
    else:
        queue_name = clust_conf['long_q']
    
    # job dependencies
    if wait_name is not None:
        hold_str = clust_conf['hold_flag'] + ' ' + str(wait_name)
        
    elif wait_file is not None:
        with open(wait_file, 'r') as wait_fi:
            wait_name = wait_fi.readline()
            hold_str = clust_conf['hold_flag'] + ' ' + str(wait_name)

    else:
        hold_str = ""
        

    # load base template
    


    # for single jobs
    if cmd is not None and (njobs is None or njobs <= 1):
                    
        njobs = 1
        tot_threads = int(threads)
        
        # log name
        if logname is None:
            logname = str(jobname)+'.sub.log'
            
        # command line
        cmd_str = cmd

        # dummy task array args for dict
        array_jobs = njobs
        j_per_core = 1


    # for array jobs
    else:

        # setup indexing tasks
        j_per_core = int(clust_conf['j_per_node'])
        if j_per_core == 1:
            task_index = str(clust_conf['task_id'])
        else:
            task_index = "${tid}"

        # cmd or array file spec
        if cmd is not None:
            cmd_line = cmd.format(task=task_index)
            tot_threads = int(njobs)*int(threads)
        
        else:
            assert os.path.isfile(arrayfile), "Job array file %s not found." % str(arrayfile)
            
            njobs = file_len(arrayfile)
            tot_threads = int(njobs)*int(threads)
            
            cmd_tmp = dedent("""\
                cline=`head -n {task} {fi} | tail -n 1`
                echo $cline
                $cline
            """)
            cmd_line = cmd_tmp.format(task=task_index, fi=arrayfile)

        # parallelization of array jobs on a node
        if j_per_core > 1:
            
            from math import floor, ceil
            
            # max simul tasks with memory limit
            node_mem = float(clust_conf['array_mem_mb'])
            task_mem_lim = floor((node_mem-1.0)/float(mem))
            
            # max simul tasks with threading
            if task_mem_lim > floor(int(j_per_core)/int(threads)):
                task_mem_lim = floor(int(j_per_core)/int(threads))
            
            if task_mem_lim < 1:
                task_mem_lim=1            
            
            # number of jobs to cover all tasks
            array_jobs = ceil(float(njobs)/float(task_mem_lim))
            
            # setup to do task_mem_lim jobs on each node
            # note: specified above that cmd_line uses ${tid} as task index 
            par_tmp = dedent("""\
                # array index for this job            
                jj={job_index}
                
                # number of jobs to run on node
                nodej={nodej}
                
                # total number of jobs to run in task array
                maxj={njobs}
                
                # task index of first task on this node
                tid=$(($nodej * ($jj - 1) + 1))
                
                # find index of last task for this node
                # - from either node task limit (nodej)
                #   or total numebr of tasks (maxj)
                if [$tid -le $(($maxj - $nodej + 1))]; then
                    last_task = $(($tid + $nodej - 1))
                else
                    last_task = $(($maxj))
                
                # start the tasks
                while [ $tid -le max_task ]; do
                    {cmd_line} &
                    tid=$(($tid+1))
                done
                
                # let all tasks finish
                wait
            """)
            
            cmd_str = par_tmp.format(njobs=str(njobs),
                                     nodej=str(task_mem_lim),
                                     job_index=str(clust_conf['task_id']),
                                     cmd_line=cmd_line)
            
            
        else:
            array_jobs = njobs
            cmd_str = cmd_line
            
            
        # log name
        if logname is None:
            logname = str(jobname)+'.sub.'+str(clust_conf['log_task_id'])+'.log'        



    # fill in template
    jobdict = {"job_name": str(jobname),
               "cmd_string": cmd_str, # formatted elsewhere
               "log_name": str(logloc)+'/'+str(logname),
               "mem_in_mb": str(mem_mb),
               "mem_in_gb": str(mem_gb),
               "threads": str(threads),
               "total_threads": str(tot_threads),
               "wall_hours": str(walltime),
               "njobs": str(njobs),
               "array_jobs": str(array_jobs),
               "array_max": str(maxpar),
               "core_par": str(j_per_core),
               "task_id": str(clust_conf['task_id']),
               "log_task_id": str(clust_conf['log_task_id']),
               "queue_name": str(queue_name),
               "sleep_time": str(sleep),
               "project": str(clust_conf['project'])
               }

            
    # write job script
    sub_file = open(str(jobname)+'.sub.sh','w')
    sub_file.write(templ.format(**jobdict))
    sub_file.close()
    
    # finalize or remove optional lines
    if njobs <= 1:
        subprocess.check_call(['sed','-i','/^::PICO_ARRAY_ONLY::/d',str(sub_file.name)])
    else:
        subprocess.check_call(['sed','-i','s/^::PICO_ARRAY_ONLY:://',str(sub_file.name)])
    
    if threads <= 1:
        subprocess.check_call(['sed','-i','/^::PICO_THREAD_ONLY::/d',str(sub_file.name)])
    else:
        subprocess.check_call(['sed','-i','s/^::PICO_THREAD_ONLY:://',str(sub_file.name)])
    
    if njobs <= 1 and threads <= 1:
        subprocess.check_call(['sed','-i','/^::PICO_THREADARRAY_ONLY::/d',str(sub_file.name)])
    else:
        subprocess.check_call(['sed','-i','s/^::PICO_THREADARRAY_ONLY:://',str(sub_file.name)])        


    # command to run
    if hold_str != "":    
        launch_str = clust_conf['sub_cmd']+' '+hold_str+' '+str(sub_file.name)
    else:
        launch_str = clust_conf['sub_cmd']+' '+str(sub_file.name)
    
    # record            
    print launch_str

    # run
    if not testonly:
        p = subprocess.Popen(launch_str.split(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        out, err = p.communicate()
        print out
        return(p.returncode)

    else:          
        return 0


####################################
#
# Get cluster configuration file
# 
####################################

def read_clust_conf():
    
    import os
    
    conf_file = os.environ['HOME']+"/picopili.conf"
    configs = read_conf(conf_file)
    cluster = configs['cluster']
    
    pico_bin = os.path.dirname(os.path.realpath(__file__))
    clust_dir = os.path.dirname(pico_bin) + '/cluster_templates'
    
    assert os.path.isdir(clust_dir), "Unable to find cluster job submission template directory %s" % str(clust_dir)

    # load queue configuration info
    # - submission syntax, queue names, job holds
    clust_conf = read_conf(str(clust_dir)+'/'+str(cluster)+'.conf')   
    
    return clust_conf


####################################
#
# Save / load job configurations
# 
####################################

def init_sendjob_dict():
    
    sendjob_dict = {
        "jobname": None,
#        "arrayfile": None,
#        "cmd": None,
        "logname": None,
        "logloc": None,
        "mem": None,
        "walltime": None,
        "njobs": None,
        "maxpar": None,
        "threads": None,
        "wait_file": None,
        "wait_name": None,
#        "cluster": None,
        "sleep": None,
#        "testonly": None
    }
    
    return sendjob_dict



def save_job(jfile, cmd_templ, job_dict, sendjob_dict):
    
    import cPickle as pickle
    
    with open(jfile, 'wb') as pickle_out:
        pickle.dump(cmd_templ, pickle_out, -1)
        pickle.dump(job_dict, pickle_out, -1)
        pickle.dump(sendjob_dict, pickle_out, -1)

    return 0



def load_job(jfile):

    import cPickle as pickle
    
    with open(jfile, 'rb') as pickle_in:
        cmd_templ = pickle.load(pickle_in)
        job_dict = pickle.load(pickle_in)
        sendjob_dict = pickle.load(pickle_in)
        
    return cmd_templ, job_dict, sendjob_dict
        
    

####################################
#
# Parse arguments from ricopili interface if invoked directly
# 
####################################
if __name__ == "__main__":

    # conditional imports    
    import argparse

    # setup arguments matching usage in imp_prep.pl
    parser = argparse.ArgumentParser(prog='blueprint.py',
                                     formatter_class=lambda prog:
                                     argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))

    arg_core = parser.add_argument_group('Job Description')

    arg_core.add_argument('--name','--na',
#                        aliases=['--na'],
                        type=str,
                        help='job name',
                        required=True)    
    arg_core.add_argument('--array',
                        type=str,
                        help='file containing command lines to be run',
                        default=None,
                        required=False)
    arg_core.add_argument('--b','-b','--blueprint',
#                        aliases=['--blueprint','--cmd'],
                        type=str,
                        help='command line to be run',
                        default=None,
                        required=False)

    
    arg_old = parser.add_argument_group('Ricopili Backwards Compatibility')
    
    arg_old.add_argument('--job','-j','--j',
#                        aliases=['--j'],
                        action='store_true',
                        help='indicates ricopili call (repurposed)')    
    parser.add_argument('--noerr',
                        action='store_true',
                        help='no output to ./errandout (for ricopili comparitibility)')
    parser.add_argument('--direct','--di',
#                        aliases=['--di'],
                        action='store_true',
                        help='start job without reading prefixes')

    arg_req = parser.add_argument_group('Resource Requirements')

    arg_req.add_argument('--mem',
                        type=int,
                        help='memeory requirement for each job, in Mb',
                        default=2000,
                        required=False)    
    arg_req.add_argument('--walltime','--wa',
 #                       aliases=['--wa'],
                        type=int,
                        help='walltime for each job, in hours',
                        default=1,
                        required=False)                            
#    arg_req.add_argument('--week',
#                        type=int,
#                        help='use week/long queues',
#                        default=None,
#                        required=False)    
    arg_req.add_argument('--njob',
                        type=int,
                        help='max number of jobs to be submitted',
                        default=1000,
                        required=False)                        
    arg_req.add_argument('--maxpar',
                        type=int,
                        help='maximum number of jobs to run in parallel',
                        default=10000,
                        required=False)
#    arg_req.add_argument('--multi',
#                        type=str,
#                        help='number of jobs to parallelize, and the number of threads to use for each parallel job (comma separated)',
#                        default=None,
#                        required=False)
    arg_req.add_argument('--fwt',
                        type=str,
                        help='file listing job dependencies to wait for before launching job',
                        default=None,
                        required=False)
    arg_req.add_argument('--wait-name',
                        type=str,
                        help='name of job dependency',
                        default=None,
                        required=False)                        

    arg_test = parser.add_argument_group('Dev Testing')
    
    parser.add_argument('--testonly',
                        action='store_true',
                        help='Skip job submission',
                        default=False)   
    
    
    args = parser.parse_args()

    # get queue
    conf_file = os.environ['HOME']+"/picopili.conf"
    configs = read_conf(conf_file)
    queue = configs['queue']
    
    # set logfile name
    if args.noerr:
        logloc = os.getcwd()
    else:
        logloc = os.getcwd()+'/errandout/'
    
    # ignore arguments for direct
    if args.direct:
        args.njob=None
        args.walltime=None
        args.mem=None
        
    
    send_job(jobname=args.name,
             arrayfile=args.array,
             cmd=args.b,
             logloc=logloc,
             mem=args.mem,
             walltime=args.walltime,
#             week=None,
             njobs=args.njob,
             maxpar=args.maxpar,
#             multi=None,
             wait_file=args.fwt,
             wait_name=args.wait_name,
             cluster=queue,
             testonly=args.testonly)

