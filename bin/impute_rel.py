#! /usr/bin/env python

####################################
# impute_rel.py
# written by Raymond Walters, January 2016
"""
Impute GWAS data with related individuals
"""
# 
# Is a wrapper for:
# - shape_rel.py
# - imp2_rel.py
# - bg_imp.py
# - agg_imp.py
# 
####################################



import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import os
from args_impute import *
from py_helpers import unbuffer_stdout #, read_conf, file_tail, link, warn_format
from blueprint import send_job
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############
parser = argparse.ArgumentParser(prog='impute_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parserphase, parserimpute, parserchunk, parserref, parserbg, parsercluster])
                    
args = parser.parse_args()


# TODO: full sanity check of the args here


# print args
print '\nBasic settings:'
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
if args.addout is not None:
    print '--addout '+str(args.addout)


print '\nReference Alignment:'
print '--popname '+str(args.popname)
print '--sfh '+str(args.sfh)
print '--fth '+str(args.fth)
print '--ref-info '+str(args.ref_info)


print '\nShapeit arguments:'
print '--window '+str(args.mem_req)
if args.no_duohmm:
    print '--no-duohmm '
print '--shape-seed '+str(args.shape_seed)


print '\nShapeit resources:'
print '--mem-req '+str(args.mem_req)
print '--threads '+str(args.threads)


print '\nIMPUTE2 arguments:'
print '--Ne '+str(args.Ne)
print '--buffer '+str(args.buffer)
if args.imp_seed is not None and str(args.imp_seed) != '' and int(args.imp_seed) > 0:
    print '--seed '+str(args.imp_seed)


print '\nGenomic chunks:'
print '--Mb-size '+str(args.Mb_size)
print '--snp_size '+str(args.snp_size)
print '--chr_info_file '+str(args.chr_info_file)


print '\nImputation Reference Files:'
print '--ref-maps '+str(args.ref_maps)
print '--ref-haps '+str(args.ref_haps)
print '--ref-legs '+str(args.ref_legs)
print '--ref-samps '+str(args.ref_samps)


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



if str(args.addout) != '' and args.addout is not None:
    outdot = str(args.out)+'.'+str(args.addout)
else:
    outdot = str(args.out)


#############
# print '\n...Checking dependencies...'
#############


# TODO: here



#############
print '\n...Submitting first task...'
#############

rp_bin = os.path.dirname(os.path.realpath(__file__))
next_call = str(rp_bin) + '/shape_rel.py '+' '.join(sys.argv[1:])+' --full-pipe'

shape_log = 'shape.'+str(outdot)+'.sub.log'

# TODO: consider queue/mem
send_job(jobname='shape.'+str(outdot),
         cmd=next_call,
         logname=shape_log,
         mem=int(args.mem_req * 1000),
         walltime=168, # week
         sleep=args.sleep)


# finish
print '\n############'
print '\n'
print 'SUCCESS!'
print 'All jobs submitted.\n'
exit(0)


# eof
