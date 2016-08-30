#! /usr/bin/env python

####################################
# agg_gwas.py
# written by Raymond Walters, December 2015
"""
Aggregates GWAS results across genomic chunks
"""
# Overview:
# 1) get chunk info
# 2) get meta information
# 3) aggregate chunks (include file of filtered top hits)
#
####################################

# TODO: logging




####################################
# Setup
# a) load python dependencies
# b) get variables/arguments
# c) read config file
# d) check dependencies
####################################

import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

import os
import subprocess
import argparse
import gzip
# from warnings import warn
# from glob import glob
from math import log10, sqrt
from args_gwas import *
from py_helpers import unbuffer_stdout, file_len
# , read_conf, link
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='agg_gwas.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parseragg])

arg_file = parser.add_argument_group('Input Files')
arg_other = parser.add_argument_group('Other Settings')

arg_file.add_argument('--chunk-file', 
                    type=str,
                    metavar='FILE',
                    help='file defining chunks used for parallelized GWAS',
                    required=True)
arg_file.add_argument('--freq-file', 
                    type=str,
                    metavar='FILE',
                    help='file with case/control allele frequencies for full data (from \'plink --freq case-control --nonfounders\')',
                    required=True)

arg_other.add_argument('--model', 
                    type=str.lower,
                    choices=['dfam','gee','gmmat','gmmat-fam'],
                    help='Which GWAS testing method was used. Current options are plink \'--dfam\' (generalized TDT-alike), GEE (generalized estimating equations), GMMAT (logistic mixed model with GRM for variance component), or GMMAT-fam (logistic mixed model with GRM and family clusters).',
                    required=False,
                    default='gee')

args = parser.parse_args()


# derived args
rp_bin = os.path.dirname(os.path.realpath(__file__))

if args.addout is not None and str(args.addout) != '':
    outdot = str(args.out)+'.'+str(args.addout)
else:
    outdot = str(args.out)

outname = outdot +'.gwas.'+str(args.model)+'.txt.gz'

bim_file = str(args.bfile) + '.bim'

logp_int = -1*int( log10(float(args.p_th2)) )
filtoutname = outdot +'.gwas.'+str(args.model)+'.p'+str(logp_int)+'_sort.txt.gz'





# TODO: report args
# TODO: check dependencies












# read chunk def file
chunks = {}
mis_chunks = {}

chunks_in = open(args.chunk_file, 'r')
dumphead = chunks_in.readline()
for line in chunks_in:
    (chrom, start, end, chname) = line.split()
    chunks[str(chname)] = [str(chrom), int(start), int(end)]

    # verify output file exists
    if args.model == 'gee':
        ch_out = 'gee.'+str(outdot)+'.'+str(chname)+'.auto.R'
    elif args.model == 'dfam':
        ch_out = 'dfam.'+str(outdot)+'.'+str(chname)+'.dfam'
    elif args.model == 'gmmat':
    	ch_out = 'gmmat_score.'+str(outdot)+'.'+str(chname)+'.R.txt'
    elif args.model == 'gmmat-fam':
   	ch_out = 'gmmatfam_score.'+str(outdot)+'.'+str(chname)+'.R.txt'
    
    # record chunks with no output
    if not os.path.isfile(ch_out):
        mis_chunks[str(chname)] = [str(chrom), int(start), int(end)]
    elif not file_len(ch_out) > 10:
        mis_chunks[str(chname)] = [str(chrom), int(start), int(end)]

chunks_in.close()

###############
# if there are missing chunks, restart their gwas and resub agg script
###############
if len(mis_chunks) > 0:
    nummiss = len(mis_chunks)
    print 'Missing results for %d GWAS jobs. Preparing to resubmit...' % nummiss
    
    # just missing chunks for task array
    # fail if already tried
    tmp_chunk_file_name = 'tmp_missing_'+str(nummiss)+'_chunks.'+str(outdot)+'.txt'

    if os.path.isfile(tmp_chunk_file_name):
        print '\n####################'
        print 'ERROR:'
        print 'Found previous attempt to resubmit %d failed chunks.' % int(nummiss)
        print 'GWAS is likely stuck.'
        print 'See %s for failing chunks.' % (tmp_chunk_file_name)
        print 'Exiting...\n'
        exit(1)

    # else setup resubmission
    tmp_chunk_file = open(tmp_chunk_file_name, 'w')
    tmp_chunk_file.write(' '.join(['CHR','START','END','NAME']) + '\n')

    for ch in mis_chunks.keys():
        tmp_chunk_file.write(' '.join([str(mis_chunks[ch][0]), str(mis_chunks[ch][1]), str(mis_chunks[ch][2]), str(ch)]) + '\n')    
    tmp_chunk_file.close()
    
    print 'List of missing chunks: %s' % tmp_chunk_file.name
    
    # copy original submit script
    # replace chunk list, name, number of tasks
    orig_uger_file = open(str(outdot)+'.gwas_chunks.sub.sh', 'r')
    new_uger_file = open(str(outdot)+'.gwas_chunks.resub_'+ str(nummiss)+'_chunks.sub.sh', 'w')
    
    for line in orig_uger_file:
        if '#$ -t ' in line:
            new_uger_file.write('#$ -t 1-'+str(nummiss)+'\n')
            continue
#        elif '#$ -tc ' in line:
#            if nummiss < 20:
#                new_uger_file.write('#$ -tc 5 \n')
#            elif nummiss < 50:
#                new_uger_file.write('#$ -tc 10 \n')
#            elif nummiss < 100:
#                new_uger_file.write('#$ -tc 25 \n')
#            else:
#                new_uger_file.write('#$ -tc 40 \n')
#            new_uger_file.write('#$ -tc 5 \n')
        elif '#$ -l m_mem_free' in line:
	    new_uger_file.write('#$ -l m_mem_free=8g \n')
        else:
            line=line.replace(args.chunk_file, tmp_chunk_file.name)
            line=line.replace('.$TASK_ID.','.tmp'+str(nummiss)+'.$TASK_ID.')
            line=line.replace('#$ -N gwas.chunks.'+str(outdot), '#$ -N gwas.chunks.'+str(outdot)+'.resub_'+str(nummiss))
            new_uger_file.write(line)
            
    orig_uger_file.close()
    new_uger_file.close()

    print ' '.join(['qsub',new_uger_file.name]) + '\n'
    subprocess.check_call(' '.join(['qsub',new_uger_file.name]), shell=True)
    print 'GWAS jobs resubmitted for %d chunks.\n' % nummiss
    
    
    print '\n...Replacing this agg job in the queue...'

    agg_log = 'agg.'+str(outdot)+'.resub_'+str(nummiss)+'.qsub.log'
    uger_agg = ' '.join(['qsub',
                            '-hold_jid','gwas.chunks.'+str(outdot)+'.resub_'+str(nummiss),
                            '-q', 'long',
                            '-l', 'm_mem_free=4g',
                            '-N', 'agg_'+str(outdot),
                            '-o', agg_log,
                            str(rp_bin)+'/uger.sub.sh',
                            str(10), # hardcoded since chunks shouldn't normally need a sleep argument
                            ' '.join(sys.argv[:])])
    
    print uger_agg + '\n'
    subprocess.check_call(uger_agg, shell=True)

    print '\n############'
    print '\n'
    print 'All jobs submitted.\n'
    exit(0)


###############
# if no missing chunks, proceed collecting info for aggregation
###############

# chnames = chunks.keys()
# sort chunk keys to aggregate in chr/bp order
chnames = [k for k, v in sorted(chunks.iteritems(), key=lambda (key,value): float(value[0]) * 1e12 + float(value[1]))]



### get meta info, index on SNP
# bim
# for gee: a2
# for dfam: bp
# gmmat: nothing
if args.model == 'gee':
    a2_info = {}
elif args.model == 'dfam':
    bp_info = {}

if args.model == 'gee' or args.model == 'dfam':
	bim = open(bim_file, 'r')
	for line in bim:
	    (chrom, snp, cm, bp, a1, a2) = line.split()
    
	    if args.model == 'gee':
	        a2_info[str(snp)] = str(a2)
	    elif args.model == 'dfam':
	        bp_info[str(snp)] = int(bp)

	bim.close()
	print 'bim loaded'

# frq.cc
# for both: maf_a, maf_u, n_a, n_u
maf_a_info = {}
maf_u_info = {}
n_a_info = {}
n_u_info = {}

frq = open(args.freq_file, 'r')
dumphead = frq.readline()
for line in frq:
    (chrom, snp, a1, a2, mafa, mafu, nchra, nchru) = line.split()
    maf_a_info[str(snp)] = mafa
    maf_u_info[str(snp)] = mafu
    n_a_info[str(snp)] = int(nchra) / 2
    n_u_info[str(snp)] = int(nchru) / 2
frq.close()
print 'frq loaded'

# info, ngt if available
info_info = {}
ngt_info = {}
if args.info_file is not None:
    info_in = open(str(args.info_file), 'r')
    dumphead = info_in.readline()
    for line in info_in:
        (snp, old_snpid, chrom, bp, a1, a2, info, exp_a1, gt) = line.split()
        info_info[str(snp)] = info
        ngt_info[str(snp)] = gt

    info_in.close()
    print 'info loaded'


### create output files
out_file = gzip.open(outname, 'wb')
filt_file = gzip.open(filtoutname+'.tmp.gz', 'wb')

if args.model == 'gee':
    out_head = ['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A', 'FRQ_U', 'INFO', 'BETA', 'SE', 'CHISQ', 'P', 'N_CAS', 'N_CON', 'ngt']
elif args.model == 'dfam':
    out_head = ['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A', 'FRQ_U', 'INFO', 'OBSERVED', 'EXPECTED', 'CHISQ', 'P', 'N_CAS', 'N_CON', 'ngt']
elif args.model == 'gmmat' or args.model == 'gmmat-fam':
    out_head = ['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A', 'FRQ_U', 'INFO', 'SCORE', 'VAR', 'Z', 'CHISQ', 'P', 'N_CAS', 'N_CON', 'ngt']

filt_head = out_head

# header
out_file.write('\t'.join(out_head) + '\n')
filt_file.write('\t'.join(filt_head) + '\n')

print 'starting chunk loop'
# loop chunks to aggregate
for ch in chnames:
    # open output file
    if args.model == 'gee':
        chunk_res = open('gee.'+str(outdot)+'.'+str(ch)+'.auto.R', 'r')
    elif args.model == 'dfam':
        chunk_res = open('dfam.'+str(outdot)+'.'+str(ch)+'.dfam', 'r')
        dumphead = chunk_res.readline()
    elif args.model == 'gmmat':
        chunk_res = open('gmmat_score.'+str(outdot)+'.'+str(ch)+'.R.txt', 'r')
        dumphead = chunk_res.readline()          
    elif args.model == 'gmmat-fam':
        chunk_res = open('gmmatfam_score.'+str(outdot)+'.'+str(ch)+'.R.txt', 'r')
        dumphead = chunk_res.readline()       
    
    for line in chunk_res:
        # read results
        if args.model == 'gee':
            (chrom, snp, bp, a1, beta, se, chisq, p, n, m) = line.split()
            a2 = a2_info.pop(str(snp))
            
        elif args.model == 'dfam':
            (chrom, snp, a1, a2, obs, exp, chisq, p) = line.split()
            bp = bp_info.pop(str(snp))
        
        elif args.model == 'gmmat' or args.model == 'gmmat-fam':
	    (chrom, snp, cm, bp, a1, a2, n, af2, scoretest, scorevar, p) = line.split()

        # get meta info
        frqa = maf_a_info.pop(str(snp))
        frqu = maf_u_info.pop(str(snp))
        na = n_a_info.pop(str(snp))
        nu = n_u_info.pop(str(snp))
        
        # info_info will be empty if no file specified
        if str(snp) in info_info:
            info = info_info.pop(str(snp))
            ngt = ngt_info.pop(str(snp))
        else:
            info = 'NA'
            ngt = 'NA'
 
 
        # filter on meta info
        # note: do here so SNPs popped off all relevant info dicts
        if float(frqa) < float(args.maf_a_th):
            continue
        elif float(frqa) > 1.0 - float(args.maf_a_th):
            continue
        elif float(frqu) < float(args.maf_u_th):
            continue
        elif float(frqu) > 1.0 - float(args.maf_u_th):
            continue
        elif str(info) != 'NA' and float(info) < float(args.info_th):
            continue
            
 
        # construct output
        if args.model == 'gee':
            # ditch gee results with implausible SEs (likely errors / numerical instability)
            if str(se) == 'NA' or float(se) > float(args.max_se):
                continue
            else:
                outline = [chrom, snp, bp, a1, a2, frqa, frqu, info, beta, se, chisq, p, na, nu, ngt]
            
        elif args.model == 'dfam':
            outline = [chrom, snp, bp, a1, a2, frqa, frqu, info, obs, exp, chisq, p, na, nu, ngt]
        
	elif args.model == 'gmmat' or args.model == 'gmmat-fam':
	    if str(p) == "NA":
	        continue
	    else:
	        z = -1.0*float(scoretest)/sqrt(float(scorevar))
		chisq = float(z)*float(z)
                outline = [chrom, snp, bp, a1, a2, frqa, frqu, info, -1.0*float(scoretest), scorevar, z, chisq, p, na, nu, ngt]

        outline = [str(val) for val in outline]        
        
        # output
        out_file.write('\t'.join(outline)+'\n')
            
        if p != 'NA' and float(p) < args.p_th2:
            filt_file.write('\t'.join(outline)+'\n')

    chunk_res.close()
    print 'chunk %s complete' % str(ch) 

out_file.close()
filt_file.close()
# final file data
# gee: chr, snp, bp, a1, a2, frq_a, frq_u, info, beta, se, chi, p, nca, nco, ngt
# dfam: chr, snp, bp, a1, a2, frq_a, frq_u, info, obs, exp, chi, p, nca, nco, ngt
# gmmat: chr, snp, bp, a1, a2, frq_a, frq_u, info, score, var, z, chi, p, nca, nco, ngt


# sort filtered file
if args.model == 'dfam' or args.model == 'gee':
    pcol = '12,12'
elif args.model == 'gmmat' or args.model == 'gmmat-fam':
    pcol = '13,13'
subprocess.check_call(' '.join([
                            'zless',
                            filtoutname+'.tmp.gz',
                            '|',
                            'sort','-g','-k',pcol,
                            '|',
                            'gzip',
                            '>', filtoutname]), shell=True)


# TODO: check completion/remove tmp?
      

print '\n############'
print '\n'
print 'SUCCESS!'
exit(0)
