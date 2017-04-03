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
from textwrap import dedent
from args_gwas import parserbase, parsergwas, parserchunk, parseragg, parsersoft
from py_helpers import link, unbuffer_stdout, find_exec, read_conf
from blueprint import send_job, read_clust_conf, init_sendjob_dict, save_job
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
elif args.model == 'logistic':
    gwas_ex = rp_bin+'/gwas_logis.py'
else:
    raise ValueError('Invalid \'--model\'. Must be one of \'gee\', \'dfam\', \'gmmat\', \'gmmat-fam\', or \'logistic\'.')


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
print '--maf-a-th '+str(args.maf_a_th)
print '--maf-u-th '+str(args.maf_u_th)
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
print '\n...Checking dependencies...'
#############

plinkx = find_exec('plink',key='p2loc')

if args.model == 'gmmat' or args.model == 'gmmat-fam':
    if args.rscript_ex == None or args.rscript_ex == "None":
        args.rscript_ex = find_exec('Rscript', key='rscloc')

# get cluster configuration
# needed for specifying logfile names with clust_conf['log_task_id']
conf_file = os.environ['HOME']+"/picopili.conf"
configs = read_conf(conf_file)
cluster = configs['cluster']
clust_conf = read_clust_conf()

# TODO: here
# TODO: move to before logging




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
        # only for gmmat, plink handles covariates for GEE/logistic
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
# TODO: consider making queue/resources flexible
######################

# basic template, depending on model
# cbopen/cbclose are placeholders for real curly braces, 
#     to survive .format() here and in send_job
if args.model == 'gee' or args.model == 'dfam' or args.model == 'logistic':
    gwas_templ = dedent("""\
    cname=`awk -v a={task} 'NR==a+1{cbopen}print $4{cbclose}' {cfile}`
    {misc}
    {gwas_ex} --bfile {bfile} --out {argout} --extract {outdot}.snps.${cbopen}cname{cbclose}.txt {optargs}    
    """)

elif args.model == 'gmmat' or args.model == 'gmmat-fam':
    gwas_templ = dedent("""\
    cname=`awk -v a={task} 'NR==a+1{cbopen}print $4{cbclose}' {cfile}`
    chrnum=`awk -v a={task} 'NR==a+1{cbopen}print $1{cbclose}' {cfile}`

    {plinkx} --bfile {bfile} --extract {outdot}.snps.${cbopen}cname{cbclose}.txt {optargs} --make-bed --out {outdot}.${cbopen}cname{cbclose}

    {rsc} --no-save --no-restore {gwas_ex} {outdot}.${cbopen}cname{cbclose} grm.{outdot}.loco_chr${cbopen}chrnum{cbclose}.rel.gz {covarsub} {outdot}.${cbopen}cname{cbclose} > {outdot}.${cbopen}cname{cbclose}.gmmat.R.log
    """)

# alternative template for GMMAT
# Rscript --no-save --no-restore 
# ~/github/picopili/bin/gmmat_logit_covar.R 
# fgwa.chr15_073_077_head               (bfile stem)
# fgw_grm_loco15.rel.gz                 (square GRM from plink matching bfile)
# ../fgwa_eur_1KGp3_postimp.pca.txt     (covariate file)
# test1                                 (output name)
# > test_gmm.log


# optional arguments
gwasargs = ''
if args.pheno is not None:
    gwasargs = gwasargs + ' --pheno '+str(args.pheno)
if args.keep is not None:
    gwasargs = gwasargs + ' --keep '+str(args.keep)
if args.remove is not None:
    gwasargs = gwasargs + ' --remove '+str(args.remove)

# model-specific arguments not passed for gmmat
if args.model == 'gee' or args.model == 'dfam' or args.model == 'logistic':
    if args.addout is not None:
        gwasargs = gwasargs + ' --addout '+str(args.addout)+'.${{cname}}'
    else:
        gwasargs = gwasargs + ' --addout ${{cname}}'
    if args.covar is not None:
        gwasargs = gwasargs + ' --covar '+str(args.covar)
    if args.covar_number is not None:
        gwasargs = gwasargs + ' --covar-number '+' '.join(args.covar_number)

# TODO: not needed for dfam, but is currently in it's args
if args.model == 'gee' or args.model == 'dfam':
    gwasargs = gwasargs + ' --r-ex '+str(args.r_ex)+' --rplink-ex '+str(args.rplink_ex)

# model specific arguments for gee to specify Rserve port for each job
# targeting IANA range 49152-65535 
# (assuming here will be < 16k jobs; gwas_gee.py handles overflow check)  
if args.model == 'gee':
    misc_txt = 'rport=$((49151+{task}))'
    gwasargs = str(gwasargs) +' --port $rport'
else:
    misc_txt = ''

# TODO: pass through cleanup


# fill in template
jobdict = {"task": "{task}",
           "cfile": chunk_file.name,
           "misc": str(misc_txt),
           "gwas_ex": str(gwas_ex),
           "bfile": str(args.bfile),
           "argout": str(args.out),
           "outdot": str(outdot),
           "optargs": str(gwasargs),
           "plinkx": str(plinkx),
           "covarsub": str(args.covar)+'.sub.txt',
           "rsc": str(args.rscript_ex),
	   "cbopen":'{{',
	   "cbclose":'}}',
           }

nchunk = len(chunks.keys())


# store job information for possible resubs
job_store_file = 'gwas.chunks.'+str(outdot)+'.pkl'

clust_dict = init_sendjob_dict()
clust_dict['jobname'] = 'gwas.chunks.'+str(outdot)
clust_dict['logname'] = str('gwas.chunks.'+str(outdot)+'.'+str(clust_conf['log_task_id'])+'.sub.log')
clust_dict['mem'] = 4000
clust_dict['walltime'] = 2
clust_dict['njobs'] = int(nchunk)
clust_dict['maxpar'] = 200
clust_dict['sleep'] = args.sleep

save_job(jfile=job_store_file, cmd_templ=gwas_templ, job_dict=jobdict, sendjob_dict=clust_dict)


# submit job
gwas_cmd = gwas_templ.format(**jobdict)

jobres = send_job(jobname='gwas.chunks.'+str(outdot),
         	  cmd=gwas_cmd,
	          logname=str('gwas.chunks.'+str(outdot)+'.'+str(clust_conf['log_task_id'])+'.sub.log'),
	          mem=4000,
	          walltime=2,
	          njobs=int(nchunk),
	          maxpar=200,
	          sleep=args.sleep)

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

agg_log = 'agg.'+str(outdot)+'.sub.log'
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


send_job(jobname='agg_'+str(outdot),
         cmd=' '.join(agg_call),
         logname=agg_log,
         mem=8000,
         walltime=30,
         wait_name='gwas.chunks.'+str(outdot),
         wait_num=str(jobres).strip(),
         sleep=args.sleep)

# TODO:
# queue summarization script (plots, etc)




# finish
print '\n############'
print '\n'
print 'SUCCESS!'
print 'All jobs submitted.\n'
exit(0)

# eof
