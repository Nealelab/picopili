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
from py_helpers import link, unbuffer_stdout, read_conf, find_from_path
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
elif args.model == 'gmmat':
    gwas_ex = rp_bin+'/gmmat_logit_covar_grmonly.R'
elif args.model == 'gmmat-fam':
    gwas_ex = rp_bin+'/gmmat_logit_covar.R'
else:
    raise ValueError('Invalid \'--model\'. Must be one of \'gee\', \'dfam\', \'gmmat\', or \'gmmat-fam\'.')


# get useful modified args
if args.addout is not None and str(args.addout) != '':
    outdir = 'gwas_'+str(args.out)+'_'+str(args.addout)
    addout_txt = ['--addout',str(args.addout)]
    outdot = str(args.out)+'.'+str(args.addout)
else:
    outdir = 'gwas_'+str(args.out)
    addout_txt = ['','']
    outdot = str(args.out)

if args.strict_bfile is None or str(args.strict_bfile) == "None":
    args.strict_bfile = str(args.bfile)

# report settings in use
print '\nBasic settings:'
print '--bfile '+str(args.bfile)
if args.model == 'gmmat' or args.model == 'gmmat-fam':
    print '--strict-bfile '+str(args.strict_bfile)
print '--out '+str(args.out)
if args.addout is not None:
    print '--addout '+str(args.addout)

print '\nAssociation Testing:'
print '--model '+str(args.model)
if args.pheno is not None:
    print '--pheno '+str(args.pheno)
print '--covar '+str(args.covar)
if args.covar_number is not None:
    print '--covar-number '+str(args.covar_number)

print '\nAnalysis Subset:'
if args.keep is not None:
    print '--keep '+str(args.keep)
else:
    print '--remove '+str(args.remove)
print '--maf-a-th '
print '--maf-u-th '
if args.info_file is not None:
    print '--info-th '+str(args.info_th)
    print '--info-file '+str(args.info_file)

print '\nParallel Jobs:'
print '--snp-chunk '+str(args.snp_chunk)

print '\nCluster Settings:'
print '--sleep '+str(args.sleep)
print '--r-ex '+str(args.r_ex)
print '--rscript-ex '+str(args.rscript_ex)
print '--rplink-ex '+str(args.rplink_ex)

# TODO: these
# note extract/exclude work in gwas scripts, but need to be handled here pre-chunking
# info score, maf thresholds left until agg since we read that data there anyway
print '\nWARNING: THESE ARGUMENTS ARE NOT CURRENTLY FUNCTIONAL FROM gwas_rel.py:'
print '--no-cleanup'
print '--extract'
print '--exclude'
# print '--info-th'
# print '--rserve-active'

print '\n'


if int(args.snp_chunk) % 100 == 0:
    warn("Setting --snp-chunk to a multiple of 100 is NOT recommended currently, due to potential plink --R errors.")


#############
print '\n...Reading ricopili config file...'
#############

### read plink loc from config
conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

plinkx = configs['p2loc']+"plink"

if args.model == 'gmmat' or args.model == 'gmmat-fam':
    if args.rscript_ex == None or args.rscript_ex == "None":
        args.rscript_ex = find_from_path('Rscript', 'Rscript')


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

if str(args.strict_bfile) != str(args.bfile):
    link(wd+'/'+str(args.strict_bfile)+'.bed', str(args.strict_bfile)+'.bed', 'plink bed file for GRM')
    link(wd+'/'+str(args.strict_bfile)+'.bim', str(args.strict_bfile)+'.bim', 'plink bim file for GRM')
    link(wd+'/'+str(args.strict_bfile)+'.fam', str(args.strict_bfile)+'.fam', 'plink fam file for GRM')

if args.covar is not None and str(args.covar) != "None":
    link(wd+'/'+str(args.covar), str(args.covar), 'covariate file')

if args.pheno is not None and str(args.pheno) != "None":
    link(wd+'/'+str(args.pheno), str(args.pheno), 'phenotype file')
    
if args.info_file is not None and str(args.info_file) != "None":
    link(wd+'/'+str(args.info_file), str(args.info_file), 'info score file')

# TODO: need to also propogate keep/exclude files, either here or in later args
# TODO: allow SNP extract/exclude exclusion before chunking




######################
print '\n...Creating genomic chunks to parallelize GWAS...'
######################

# create chunks
# note: <=16000 chunks is ideal for targeted range of Rserve port numbers
#       and is enough for chunks >500 SNPs with 8 million imputed markers
#       but can in theory handle up to 64510 (ie. 1025-65534) if needed before 
#       port collisions become an issue
chunk_call = [chunker_ex,
              '--bfile',str(args.bfile),
              '--out',str(args.out),
              addout_txt[0],addout_txt[1],
              '--Mb-size',str(1),
              '--snp-size',str(args.snp_chunk),
              '--ignore-centromeres',
              '--allow-small-chunks',
              '--max-chunks',str(64000)]
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
# track chrs (for GRM building in gmmat)
chrs = []

# find chunk for each SNP, write to file
bim = open(str(args.bfile)+'.bim', 'r')
for line in bim:
    (chrom, snp, cm, bp, a1, a2) = line.split()

    # reocrd novel chrs
    if int(chrom) not in chrs:
        chrs.append(int(chrom))

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
if args.model == 'gmmat' or args.model == 'gmmat-fam':
    print '\n...Building GRMs for mixed model (leave-one-chromosome-out)...'
######################  
    
    for ch in xrange(1,23):
        # skip for chromosomes with no data
        if ch not in chrs:
            continue

        # verify file doesn't already exist        
        grm_file = 'grm.'+str(outdot)+'.loco_chr'+str(ch)+'.rel.gz'
        if not os.path.isfile(str(grm_file)):
            grm_call = [plinkx,
                        '--bfile',str(args.strict_bfile),
                        '--maf',str(.01),
                        '--geno', str(.01),                        
                        '--make-rel','square','gz',
                        '--make-founders',
                        '--autosome',
                        '--not-chr',str(ch),
                        '--out','grm.'+str(outdot)+'.loco_chr'+str(ch)]
    
            if args.keep is not None:
                grm_call.extend(['--keep',str(args.keep)])
            elif args.remove is not None:
                grm_call.extend(['--remove',str(args.remove)])     

            # run plink, keep log
            grm_log = open('grm.'+str(outdot)+'.loco_chr'+str(ch)+'.plink.log', 'w')
            subprocess.check_call(grm_call, stderr=subprocess.STDOUT, stdout=grm_log)
            grm_log.close()
            print 'Created GRM omitting chr %d: %s' % (int(ch), str(grm_file))

            # verify successfully created
            if not os.path.isfile(str(grm_file)):
                raise IOError('Failed to create GRM omitting chr %d: %s' % (int(ch), str(grm_file)))
        
        # log in case existing file found
        else:
            print 'Found existing GRM omitting chr %d: %s' % (int(ch), str(grm_file))



######################
    if args.covar is not None:
        print '\n...Preparing covariates...' 
        # only for gmmat, plink handles covariates
######################

        cov_in = open(str(args.covar), 'r')
        covhead = cov_in.readline().split()
        
        if 'FID' not in covhead or 'IID' not in covhead:
            ValueError('Covariate file is missing header with FID and IID (required for GMMAT).')
    
        # extract covariates to use, if necessary        
        if args.covar_number is not None:
            cov_out = open(str(args.covar)+'.sub.txt', 'w')
            
            # split eg. 1-10,12
            cov1 = [x.split(',') for x in args.covar_number]
    
            # flatten list of lists
            cov2 = [item for sublist in cov1 for item in sublist]
            
            # remove empty strings
            cov3 = [x.replace(',','') for x in cov2 if x is not None and str(x) != '']
    
            # interpret ranges
            # credit: http://stackoverflow.com/questions/6405208/how-to-convert-numeric-string-ranges-to-a-list-in-python
            covnum = []
            for part in cov3:
                    if '-' in part:
                            a, b = part.split('-')
                            a, b = int(a), int(b)
                            covnum.extend(range(a, b+1))
                    else:
                            a = int(part)
                            covnum.append(a)
    
            abscol = [0, 1]
            # +1 here is +2 for FID/IID cols, -1 for python's base-0 index
            abscol.extend([x+1 for x in covnum])
    
            # write header line
            head_out = ' '.join([str(covhead[int(x)]) for x in abscol])
            cov_out.write(head_out + '\n')
    
            # write data lines
            for line in cov_in:
                line_out = ' '.join([str(line.split()[int(x)]) for x in abscol])
                cov_out.write(line_out + '\n')
            
            cov_out.close()
            cov_in.close()
            
        else:
            cov_in.close()
        
    # no covariates
    else:
        raise ValueError('GMMAT without covariates not currently implemented.\n')





######################
print '\n...Submitting GWAS for all chunks...'
######################

# gwas each chunk
# need to write submit script to include chunk name parsing
# TODO: consider making queue/resources flexible

if args.model == 'gee' or args.model == 'dfam':
    uger_gwas_template = """#!/usr/bin/env sh
#$ -j y
#$ -cwd
#$ -V
#$ -N {jname}
#$ -q short
#$ -l m_mem_free=4g
#$ -t 1-{nchunk}
#$ -tc 200
#$ -o {outlog}

source /broad/software/scripts/useuse
reuse -q Anaconda
sleep {sleep}

cname=`awk -v a=${{SGE_TASK_ID}} 'NR==a+1{{print $4}}' {cfile}`

{misc}

{gwas_ex} --bfile {bfile} --out {argout} --extract {outdot}.snps.${{cname}}.txt {optargs}

# eof
"""

# alternative template for GMMAT
# Rscript --no-save --no-restore 
# ~/github/picopili/bin/gmmat_logit_covar.R 
# fgwa.chr15_073_077_head               (bfile stem)
# fgw_grm_loco15.rel.gz                 (square GRM from plink matching bfile)
# ../fgwa_eur_1KGp3_postimp.pca.txt     (covariate file)
# test1                                 (output name)
# > test_gmm.log
elif args.model == 'gmmat' or args.model == 'gmmat-fam':
    uger_gwas_template = """#!/usr/bin/env sh
#$ -j y
#$ -cwd
#$ -V
#$ -N {jname}
#$ -q short
#$ -l m_mem_free=4g
#$ -t 1-{nchunk}
#$ -tc 200
#$ -o {outlog}

source /broad/software/scripts/useuse
sleep {sleep}

cname=`awk -v a=${{SGE_TASK_ID}} 'NR==a+1{{print $4}}' {cfile}`
chrnum=`awk -v a=${{SGE_TASK_ID}} 'NR==a+1{{print $1}}' {cfile}`

{plinkx} --bfile {bfile} --extract {outdot}.snps.${{cname}}.txt {optargs} --make-bed --out {outdot}.${{cname}}

{rsc} --no-save --no-restore {gwas_ex} {outdot}.${{cname}} grm.{outdot}.loco_chr${{chrnum}}.rel.gz {covarsub} {outdot}.${{cname}} > {outdot}.${{cname}}.gmmat.R.log

# eof
"""


gwasargs = ''

if args.pheno is not None:
    gwasargs = gwasargs + ' --pheno '+str(args.pheno)
if args.keep is not None:
    gwasargs = gwasargs + ' --keep '+str(args.keep)
if args.remove is not None:
    gwasargs = gwasargs + ' --remove '+str(args.remove)

# these args not passed for gmmat
if args.model == 'gee' or args.model == 'dfam':
    if args.addout is not None:
        gwasargs = gwasargs + ' --addout '+str(args.addout)+'.${cname}'
    else:
        gwasargs = gwasargs + ' --addout ${cname}'
    if args.covar is not None:
        gwasargs = gwasargs + ' --covar '+str(args.covar)
    if args.covar_number is not None:
        gwasargs = gwasargs + ' --covar-number '+' '.join(args.covar_number)

    gwasargs = gwasargs + ' --r-ex '+str(args.r_ex)+' --rplink-ex '+str(args.rplink_ex)

# TODO: pass through cleanup

    
nchunk = len(chunks.keys())
jobdict = {"jname": 'gwas.chunks.'+str(outdot),
           "nchunk": str(nchunk),
           "outlog": str('gwas.chunks.'+str(outdot)+'.$TASK_ID.qsub.log'),
           "sleep": str(args.sleep),
           "cfile": chunk_file.name,
           "misc": '',
           "gwas_ex": str(gwas_ex),
           "bfile": str(args.bfile),
           "argout": str(args.out),
           "outdot": str(outdot),
           "optargs": str(gwasargs),
           "plinkx": str(plinkx),
           "covarsub": str(args.covar)+'.sub.txt',
           "rsc": str(args.rscript_ex)
           }

# for gee, need to specify Rserve port for each job
# targeting IANA range 49152-65535 
# (assuming here will be < 16k jobs; gwas_gee.py handles overflow check)           
if args.model == 'gee':
    jobdict['misc'] = 'rport=$((49151+SGE_TASK_ID))'
    jobdict['optargs'] = str(gwasargs) +' --port $rport'


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
            '--silent',
            '--memory', str(2000),
            '--out','freqinfo.'+str(outdot)]
if args.keep is not None:
    frq_call.extend(['--keep',str(args.keep)])
elif args.remove is not None:
    frq_call.extend(['--remove',str(args.remove)])
if args.pheno is not None and str(args.pheno) != "None":
    frq_call.extend(['--pheno',str(args.pheno)])

freqname = 'freqinfo.'+str(outdot)+'.frq.cc'

frq_log = open('freqinfo.'+str(outdot)+'.plink.log', 'w')
subprocess.check_call(frq_call, stderr=subprocess.STDOUT, stdout=frq_log)
frq_log.close()


######################
print '\n...Queuing GWAS results aggregation script...'
# call to agg_gwas.py
######################

if args.info_file is not None and str(args.info_file) != "None":
    info_file_txt = ['--info-file',str(args.info_file),'--info-th',str(args.info_th)]
else:
    info_file_txt = ['','','','']

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
            info_file_txt[0],info_file_txt[1],info_file_txt[2],info_file_txt[3],
            '--max-se',str(args.max_se),
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

print uger_agg + '\n'
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
