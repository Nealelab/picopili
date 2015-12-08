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




import subprocess
import os
from warnings import warn
from py_helpers import link

bfile = 'plink_bfile'
rp_bin = '/home/unix/rwalters/github/picopili/rp_bin'

gwas_ex = rp_bin + '/gwas_gee.py'
chunker_ex = rp_bin+'/chunk_snps.py'
snp_size = 20000
out = 'test'
addout = 'test1'
sleep = 5
covar = None
covar_number = None
keep = None
remove = None


if addout is not None and str(addout) != '':
    outdir = 'gwas_'+str(out)+'_'+str(addout)
    addout_txt = ['--addout',str(addout)]
    outdot = str(out)+'.'+str(addout)
else:
    outdir = 'gwas_'+str(out)
    addout_txt = ['','']
    outdot = str(out)


# setup working directory
wd = os.getcwd()
os.mkdir(outdir)
os.chdir(outdir)
link(wd+'/'+str(bfile)+'.bed', str(bfile)+'.bed', 'input plink bed file')
link(wd+'/'+str(bfile)+'.bim', str(bfile)+'.bim', 'input plink bim file')
link(wd+'/'+str(bfile)+'.fam', str(bfile)+'.fam', 'input plink fam file')



# TODO: allow SNP extract/exclude exclusion before chunking


# create chunks
chunk_call = [chunker_ex,
              '--bfile',str(bfile),
              '--out',str(out),
              addout_txt[0],addout_txt[1],
              '--Mb-size',str(1),
              '--snp-size',str(snp_size),
              '--ignore-centromeres',
              '--allow-small-chunks']
chunk_call = filter(None,chunk_call)

chunk_log = open('chunk.'+str(outdot)+'.log', 'w')
print ' '.join(chunk_call) + '\n'
subprocess.check_call(chunk_call, stderr=subprocess.STDOUT, stdout=chunk_log)
chunk_log.close()


# create keep lists for each chunk
chunk_file = open(outdot+'.chunks.txt', 'r')


print "Loading chunks\n"
# read in all chunk info
chunks = {}
dumphead = chunk_file.readline()
for line in chunk_file:
    (chrom, start, end, chname) = line.split()
    chunks[str(chname)] = [str(chrom), int(start), int(end)]
    
chunk_file.close()

print "Filtering SNPs\n"
# filter SNPs into chunks
bim = open(str(bfile)+'.bim', 'r')
omitsnp = 0
cur_chunk = open(str(outdot)+'.snps.'+str(chunks.keys()[0])+'.txt', 'w')
last_chunk = str(chunks.keys()[0])

# fast chunk-finding
# note: the shortcut of trying the last_chunk first means overlappings chunks are unlikely to be caught
def find_chunk(snpchrom, snpbp, last_chunk):
    if str(snpchrom) == str(chunks[last_chunk][0]) and int(snpbp) >= int(chunks[last_chunk][1]) and int(snpbp) <= int(chunks[last_chunk][2]):
        return [last_chunk]
    else:
        return [key for key,value in chunks.iteritems() if str(snpchrom)==str(value[0]) and int(snpbp) >= int(value[1]) and int(snpbp) <= int(value[2])]

# find for each SNP, write to file
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
    # use name check to avoid open/close file for each SNP 
    # - (assume SNPs likely sorted, so few switches here)
    else:
        if cur_chunk.name != str(outdot)+'.snps.'+str(snp_chunk[0])+'.txt':
            cur_chunk.close()
            cur_chunk = open(str(outdot)+'.snps.'+str(snp_chunk[0])+'.txt', 'a')
            
        cur_chunk.write(str(snp)+'\n')

bim.close()
if omitsnp > 0:
    omitfile.close()
    warn('%d SNPs not assigned to any chunks. Unexpected chromosome locations? List written to %s' % (int(omitsnp), str(outdot)+'.omitted_snps.txt'))




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

cname=`awk -v a=${{SGE_TASK_ID}} 'NR==a{{print $4}}' {cfile}`

{gwas_ex} --bfile {bfile} --out {argout} --extract {outdot}.snps.${{cname}}.txt {optargs}

# eof
"""
gwasargs = ''
if addout is not None:
    gwasargs = gwasargs + ' --addout '+str(addout)+'.${cname}'
else:
    gwasargs = gwasargs + ' --addout ${cname}'
if covar is not None:
    gwasargs = gwasargs + ' --covar '+str(covar)
if covar_number is not None:
    gwasargs = gwasargs + ' --covar-number '+' '.join(covar_number)
if keep is not None:
    gwasargs = gwasargs + ' --keep '+str(keep)
if remove is not None:
    gwasargs = gwasargs + ' --remove '+str(remove)
# TODO: pass through cleanup/Rserve/executables
    
nchunk = len(chunks.keys())
jobdict = {"jname": 'gwas.chunks.'+str(outdot),
           "nchunk": str(nchunk),
           "outlog": str('gwas.chunks.'+str(outdot)+'.$TASK_ID.log'),
           "sleep": str(sleep),
           "cfile": chunk_file.name,
           "gwas_ex": str(gwas_ex),
           "bfile": str(bfile),
           "argout": str(out),
           "outdot": str(outdot),
           "optargs": str(gwasargs)
           }

uger_gwas = open(str(outdot)+'.gwas_chunks.sub.sh', 'w')
uger_gwas.write(uger_gwas_template.format(**jobdict))
uger_gwas.close()

print ' '.join(['qsub',uger_gwas.name]) + '\n'
subprocess.check_call(' '.join(['qsub',uger_gwas.name]), shell=True)
print 'GWAS jobs submitted for %d chunks.\n' % nchunk





# TODO:
# queue aggregation script
# queue summarization script (plots, etc)




# finish
print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)

#uger_chunk = ' '.join(['qsub',
#                        '-hold_jid',str(name),
#                        '-q', 'short',
#                        '-N', str('chunk_'+out),
#                        '-o', chunk_log,
#                        str(rp_bin)+'/uger.sub.sh',
#                        str(sleep),
#                        str(chunk_call)])