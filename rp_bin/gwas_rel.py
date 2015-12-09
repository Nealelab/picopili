#! /usr/bin/env python

####################################
# gwas_rel.py
# written by Raymond Walters, December 2015
"""
Runs GWAS of plink family data
"""
# Overview:
# 1) split files for parallelization
# 2) gwas chunks
# 3) aggregate and format (add a2, maf, info, etc)
# 4) plots/output
#
####################################


import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import argparse
import subprocess
import os
from warnings import warn
from args_gwas import *
from py_helpers import link, unbuffer_stdout, read_conf
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='gwas_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parsergwas, parserchunk, parseragg, parsersoft])
                    
args = parser.parse_args()


# get directory containing current script
# (to get absolute path for scripts)
rp_bin = os.path.dirname(os.path.realpath(__file__))
chunker_ex = rp_bin+'/chunk_snps.py'

# get desired gwas task script
if args.model == 'gee':
    gwas_ex = rp_bin+'/gwas_gee.py'
elif args.model == 'dfam':
    gwas_ex = rp_bin+'/gwas_dfam.py'
else:
    raise ValueError('Invalid \'--model\'. Must be one of \'gee\' or \'dfam\'.')


# get useful modified args
if args.addout is not None and str(args.addout) != '':
    outdir = 'gwas_'+str(args.out)+'_'+str(args.addout)
    addout_txt = ['--addout',str(args.addout)]
    outdot = str(args.out)+'.'+str(args.addout)
else:
    outdir = 'gwas_'+str(args.out)
    addout_txt = ['','']
    outdot = str(args.out)


# report settings in use
print '\nBasic settings:'
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
if args.addout is not None:
    print '--addout '+str(args.addout)

print '\nAssociation Testing:'
print '--model '+str(args.model)
print '--covar '+str(args.covar)
if args.covar_number is not None:
    print '--covar-number '+str(args.covar_number)

print '\nAnalysis Subset:'
if args.keep is not None:
    print '--keep '+str(args.keep)
else:
    print '--remove '+str(args.remove)

print '\nParallel Jobs:'
print '--snp-chunk '+str(args.snp_chunk)

print '\nCluster Settings:'
print '--sleep '+str(args.sleep)


# TODO: these
print '\nWARNING: THESE ARGUMENTS ARE NOT CURRENTLY FUNCTIONAL FROM gwas_rel.py:'
print '--no-cleanup'
print '--extract'
print '--exclude'
print '--info-th'
print '--rserve-active'
print '--r-ex'
print '--rplink-ex'
print '\n'


#############
print '\n...Reading ricopili config file...'
#############

### read plink loc from config
conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

plinkx = configs['p2loc']+"plink"


#############
print '\n...Checking dependencies...'
#############



# TODO: here





print '\n'
print '############'
print 'Begin!'
print '############'

######################
print '\n...Setting up working directory ./%s...' % str(outdir)
######################

wd = os.getcwd()
os.mkdir(outdir)
os.chdir(outdir)
link(wd+'/'+str(args.bfile)+'.bed', str(args.bfile)+'.bed', 'input plink bed file')
link(wd+'/'+str(args.bfile)+'.bim', str(args.bfile)+'.bim', 'input plink bim file')
link(wd+'/'+str(args.bfile)+'.fam', str(args.bfile)+'.fam', 'input plink fam file')




# TODO: need to also propogate keep/exclude files, either here or in later args
# TODO: allow SNP extract/exclude exclusion before chunking




######################
print '\n...Creating genomic chunks to parallelize GWAS...'
######################

# create chunks
chunk_call = [chunker_ex,
              '--bfile',str(args.bfile),
              '--out',str(args.out),
              addout_txt[0],addout_txt[1],
              '--Mb-size',str(1),
              '--snp-size',str(args.snp_chunk),
              '--ignore-centromeres',
              '--allow-small-chunks']
chunk_call = filter(None,chunk_call)

chunk_log = open('chunk.'+str(outdot)+'.log', 'w')
print ' '.join(chunk_call) + '\n'
subprocess.check_call(chunk_call, stderr=subprocess.STDOUT, stdout=chunk_log)
chunk_log.close()


# read in chunk info from output
chunks = {}

chunk_file = open(outdot+'.chunks.txt', 'r')
dumphead = chunk_file.readline()
for line in chunk_file:
    (chrom, start, end, chname) = line.split()
    chunks[str(chname)] = [str(chrom), int(start), int(end)]

chunk_file.close()
print 'Chunk definitions written to %s' % './'+str(outdir)+'/'+str(chunk_file.name)

######################
print '\n...Dividing SNPs into genomic chunks...'
######################

# fast chunk-finding
# note: the shortcut of trying the last_chunk first means overlappings chunks are unlikely to be caught
def find_chunk(snpchrom, snpbp, last_chunk):
    if str(snpchrom) == str(chunks[last_chunk][0]) and int(snpbp) >= int(chunks[last_chunk][1]) and int(snpbp) <= int(chunks[last_chunk][2]):
        return [last_chunk]
    else:
        return [key for key,value in chunks.iteritems() if str(snpchrom)==str(value[0]) and int(snpbp) >= int(value[1]) and int(snpbp) <= int(value[2])]

# track current chunk name/file, number of SNPs without a chunk assignment
cur_chunk = open(str(outdot)+'.snps.'+str(chunks.keys()[0])+'.txt', 'w')
last_chunk = str(chunks.keys()[0])
omitsnp = 0

# find chunk for each SNP, write to file
bim = open(str(args.bfile)+'.bim', 'r')
for line in bim:
    (chrom, snp, cm, bp, a1, a2) = line.split()

    # find matching chunk    
    snp_chunk = find_chunk(chrom, bp, last_chunk)
    last_chunk = str(snp_chunk[0])

    # sanity checks
    if len(snp_chunk) > 1:
        raise ValueError('SNP %s matches multiple chunks (%r). Chunks should be non-overlapping?' % (str(snp), snp_chunk))
    elif len(snp_chunk) == 0:
        if omitsnp == 0:
            omitfile = open(str(outdot)+'.omitted_snps.txt', 'w')
        omitfile.write(str(snp)+'\n')
    # write chunk file
    # use name check to avoid opening/closing file for each SNP 
    # - (assumes SNPs are likely sorted, so few file switches here)
    else:
        if cur_chunk.name != str(outdot)+'.snps.'+str(snp_chunk[0])+'.txt':
            cur_chunk.close()
            cur_chunk = open(str(outdot)+'.snps.'+str(snp_chunk[0])+'.txt', 'a')
            
        cur_chunk.write(str(snp)+'\n')

bim.close()
print 'SNP lists for each chunk written to %s' % './'+str(outdir)+'/'+str(outdot)+'.snps.[CHUNK].txt'
if omitsnp > 0:
    omitfile.close()
    warn('%d SNPs not assigned to any chunks. Unexpected chromosome locations? List written to %s' % (int(omitsnp), './'+str(outdir)+'/'+str(outdot)+'.omitted_snps.txt'))



######################
print '\n...Submitting GWAS for all chunks...'
######################

# gwas each chunk
# need to write submit script to include chunk name parsing
# TODO: consider making queue/resources flexible
uger_gwas_template = """#!/usr/bin/env sh
#$ -j y
#$ -cwd
#$ -V
#$ -N {jname}
#$ -q short
#$ -l m_mem_free=2g
#$ -t 1-{nchunk}
#$ -o {outlog}

source /broad/software/scripts/useuse
reuse -q Anaconda
sleep {sleep}

cname=`awk -v a=${{SGE_TASK_ID}} 'NR==a+1{{print $4}}' {cfile}`

{gwas_ex} --bfile {bfile} --out {argout} --extract {outdot}.snps.${{cname}}.txt {optargs}

# eof
"""
gwasargs = ''
if args.addout is not None:
    gwasargs = gwasargs + ' --addout '+str(args.addout)+'.${cname}'
else:
    gwasargs = gwasargs + ' --addout ${cname}'
if args.covar is not None:
    gwasargs = gwasargs + ' --covar '+str(args.covar)
if args.covar_number is not None:
    gwasargs = gwasargs + ' --covar-number '+' '.join(args.covar_number)
if args.keep is not None:
    gwasargs = gwasargs + ' --keep '+str(args.keep)
if args.remove is not None:
    gwasargs = gwasargs + ' --remove '+str(args.remove)
# TODO: pass through cleanup/Rserve/executables
    
nchunk = len(chunks.keys())
jobdict = {"jname": 'gwas.chunks.'+str(outdot),
           "nchunk": str(nchunk),
           "outlog": str('gwas.chunks.'+str(outdot)+'.$TASK_ID.qsub.log'),
           "sleep": str(args.sleep),
           "cfile": chunk_file.name,
           "gwas_ex": str(gwas_ex),
           "bfile": str(args.bfile),
           "argout": str(args.out),
           "outdot": str(outdot),
           "optargs": str(gwasargs)
           }

uger_gwas = open(str(outdot)+'.gwas_chunks.sub.sh', 'w')
uger_gwas.write(uger_gwas_template.format(**jobdict))
uger_gwas.close()

print ' '.join(['qsub',uger_gwas.name]) + '\n'
subprocess.check_call(' '.join(['qsub',uger_gwas.name]), shell=True)
print 'GWAS jobs submitted for %d chunks.\n' % nchunk



######################
print '\n...Preparing meta-data for aggregation...'
# - create .frq file for aggregation script
# TODO: any prep for info score file?
######################

frq_call = [plinkx,
            '--bfile',str(args.bfile),
            '--freq','case-control','--nonfounders',
            '--out','freqinfo.'+str(outdot)]
if args.keep is not None:
    frq_call.extend(['--keep',str(args.keep)])
elif args.remove is not None:
    frq_call.extend(['--remove',str(args.remove)])

freqname = 'freqinfo.'+str(outdot)+'.frq.cc'

frq_log = open('freqinfo.'+str(outdot)+'.plink.log', 'w')
subprocess.check_call(frq_call, stderr=subprocess.STDOUT, stdout=frq_log)
frq_log.close()


######################
print '\n...Queuing GWAS results aggregation script...'
# call to agg_gwas.py
# TODO: info score file
######################

agg_log = 'agg.'+str(outdot)+'.qsub.log'
agg_call = [str(rp_bin)+'/agg_gwas.py',
            '--bfile',str(args.bfile),
            '--out',str(args.out),
            addout_txt[0],addout_txt[1],
            '--maf-a-th',str(args.maf_a_th),
            '--maf-u-th',str(args.maf_u_th),
            '--p-th2',str(args.p_th2),
            '--chunk-file',str(chunk_file.name),
            '--freq-file',str(freqname),
            '--model',str(args.model)]
agg_call = filter(None,agg_call)

uger_agg = ' '.join(['qsub',
                        '-hold_jid','gwas.chunks.'+str(outdot),
                        '-q', 'long',
                        '-l', 'm_mem_free=4g',
                        '-N', 'agg_'+str(outdot),
                        '-o', agg_log,
                        str(rp_bin)+'/uger.sub.sh',
                        str(args.sleep),
                        ' '.join(agg_call)])

subprocess.check_call(uger_agg, shell=True)


# TODO:
# queue summarization script (plots, etc)




# finish
print '\n############'
print '\n'
print 'SUCCESS!'
print 'All jobs submitted.\n'
exit(0)

#uger_chunk = ' '.join(['qsub',
#                        '-hold_jid',str(name),
#                        '-q', 'short',
#                        '-N', str('chunk_'+out),
#                        '-o', chunk_log,
#                        str(rp_bin)+'/uger.sub.sh',
#                        str(sleep),
#                        str(chunk_call)])