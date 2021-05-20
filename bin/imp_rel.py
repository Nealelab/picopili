#! /usr/bin/env python

####################################
# imp_rel.py
# written by Raymond Walters, January 2016
"""
Runs IMPUTE 2 or 4 for GWAS data with related individuals
"""
# Overview:
# 1) Parse arguments
#    - check dependencies, print args, etc
# 2) Define genomic chunks
# 3) Run impute
#
# TODO: better reference specification options
# TODO: chrx imputation
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
from textwrap import dedent
from args_impute import parserbase, parserimpute, parserref, parserchunk, parsercluster, parserjob, parserphase
from py_helpers import unbuffer_stdout, file_len, link, find_exec, test_exec, read_conf
from blueprint import send_job, read_clust_conf, init_sendjob_dict, save_job
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############
parser = argparse.ArgumentParser(prog='imp_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parserimpute, parserref, parserchunk, parsercluster, parserjob, parserphase])

parser.add_argument('--ref-dir',
			type=str,
			metavar='DIRECTORY',
			help='Directory containing imputation reference files (haps, legends, sample, and maps). ' +
				'Used as prefix for specifying full paths of --ref-maps, --ref-haps, --ref-legs, and --ref-samps',
			required=False,
			default=None)


args, extra_args = parser.parse_known_args()


# get useful modified args
if args.addout is not None and str(args.addout) != '':
    addout_txt = ['--addout',str(args.addout)]
    outdot = str(args.out)+'.'+str(args.addout)
else:
    addout_txt = ['','']
    outdot = str(args.out)
    
if args.imp_seed is not None and str(args.imp_seed) != '' and int(args.imp_seed) > 0:
    if args.imp_version==4:
        print "\n\nWARNING: --seed is ignored for IMPUTE4 \n"
        seedtxt = ''
    else:
        seedtxt = '-seed '+str(args.imp_seed)
else:
    seedtxt = ''


# report settings in use
print '\nBasic settings:'
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
if args.addout is not None:
    print '--addout '+str(args.addout)

print '\nIMPUTE arguments:'
print '--imp-version '+str(args.imp_version)
print '--Ne '+str(args.Ne)
print '--buffer '+str(args.buffer)
if seedtxt != '':
    print '--seed '+str(args.imp_seed)

print '\nImputation reference files:'
print '--ref-maps '+str(args.ref_maps)
print '--ref-haps '+str(args.ref_haps)
print '--ref-legs '+str(args.ref_legs)
print '--ref-samps '+str(args.ref_samps)

print '\nGenomic chunks:'
print '--Mb-size '+str(args.Mb_size)
print '--snp_size '+str(args.snp_size)
print '--chr_info_file '+str(args.chr_info_file)

print '\nCluster settings:'
print '--sleep '+str(args.sleep)
if args.full_pipe:
    print '--full-pipe'


#############
print '\n...Checking dependencies...'
#############

# get cluster configuration
# needed for specifying logfile names with clust_conf['log_task_id']
conf_file = os.environ['HOME']+"/picopili.conf"
configs = read_conf(conf_file)
cluster = configs['cluster']
clust_conf = read_clust_conf()

# from config
if args.imp_version==2:
    impute_ex = find_exec('impute2',key='i2loc')
elif args.imp_version==4:
    impute_ex = find_exec('impute4',key='i4loc')
    qctool_ex = find_exec('qctool', key='qctoolloc')
else:
    raise ValueError("Currently --imp-version can only be 2 or 4")
shapeit_ex = find_exec('shapeit',key='shloc')

# get directory containing current script
# (to get absolute path for scripts)
rp_bin = os.path.dirname(os.path.realpath(__file__))
chunker_ex = rp_bin+'/chunk_snps.py'
test_exec(chunker_ex,'picopili chunking script')

if args.ref_dir is not None:
	# verify exists
	assert os.path.isdir(args.ref_dir), "Failed to find imputation reference directory %s" % args.ref_dir

	# prepend to references accordingly
	args.ref_maps = str(args.ref_dir) +'/' + args.ref_maps
	args.ref_haps = str(args.ref_dir) +'/' + args.ref_haps
	args.ref_legs = str(args.ref_dir) +'/' + args.ref_legs
	args.ref_samps = str(args.ref_dir) +'/' + args.ref_samps


# TODO: here
# .hg19.ch.fl.bim for chunking
# imp. references
# executables


# directories
wd = os.getcwd()
shape_dir = wd + '/phase_chr'




print '\n'
print '############'
print 'Begin!'
print '############'



######################
print '\n...Verifying pre-phasing was successful...'
######################

bad_chr = []

for chrom in xrange(1,23):
    haps_out = str(shape_dir)+'/'+str(outdot)+'.chr'+str(chrom)+'.phased.haps'
    samp_out = str(shape_dir)+'/'+str(outdot)+'.chr'+str(chrom)+'.phased.sample'

    if not os.path.isfile(haps_out):
        bad_chr.append(str(chrom))
    elif not os.path.isfile(samp_out):
        bad_chr.append(str(chrom))
        

# if any shapeit jobs failed, 
# resubmit them and re-queue this job
if bad_chr:    
    num_chr = len(bad_chr)
    print 'Missing pre-phasing results for %d chromosomes.' % num_chr
    
    if not args.full_pipe:
        print 'Missing: %s' % ','.join(bad_chr)
        print 'Exiting\n'
        exit(1)
    # else continue to resub
    print 'Preparing to resubmit...'
    # note: assuming required shapeit args will be in args
    #   if running under --full-pipe
    #   TODO: add check on this
    #   (mem_req, threads, no_duohmm, window, shape_seed) 
    
    os.chdir(shape_dir)
        
    # verify haven't already tried this resub
    phase_sub_name = 'shapeit.'+str(outdot)+'.resub_'+str(num_chr)+'.sub.sh'
    if os.path.isfile(phase_sub_name):
        print '\n####################'
        print 'ERROR:'
        print 'Found previous attempt to resubmit %d failed chromosomes.' % int(num_chr)
        print 'Pre-phasing is likely stuck.'
        print 'Problem chromosomes: %s' % (','.join(bad_chr))
        print 'Exiting...\n'
        exit(1)

    # setup submit script
    # with "chr_list" to get have adaptive chromosome list
    cmd_templ = dedent("""\
    chrs=({chr_list})
    chrom=${cbopen}chrs[{task}-1]{cbclose}

    {shape_ex} {bed} {map} {ph_ref_txt} {window} {duo_txt} {thread_str} {seed_str} {outmax} {shapelog}    
    """)

#    shape_call = [shapeit_ex,
#                  '--input-bed', chrstem+'.bed', chrstem+'.bim', chrstem+'.fam',
#                  '--input-map', args.map_dir+'/genetic_map_chr\$tasknum_combined_b37.txt',
#                  '--input-ref', refstem+'_chr\$tasknum.hap.gz', refstem+'_chr\$tasknum.legend.gz', refstem+'.sample',
#                  '--window', str(args.window),
#                  '--duohmm',
#                  '--thread', str(args.threads),
#                  '--seed', str(args.shape_seed),
#                  '--output-max', outstem+'.phased.haps', outstem+'.phased.sample',
#                  '--output-log', outstem+'.shape.log']
    
    # manage additional arg pieces
    chrstem = str(args.bfile)+'.hg19.ch.fl.chr${chrom}'
    outstem = str(outdot)+'.chr${chrom}'
    if args.no_duohmm:
        duo_txt = ''
    else:
        duo_txt = '--duohmm'

if args.no_phaseref:
    ph_ref_txt =''
else:
    ph_ref_txt ='--input-ref '+ \
                    str(args.ref_haps).replace('###','${chrom}')+' '+ \
                    str(args.ref_legs).replace('###','${chrom}')+' '+ \
                    str(args.ref_samps).replace('###','${chrom}')
    
    # fill in shapeit template
    jobdict = {"task": "{task}",
               "chr_list": ' '.join(bad_chr),
               "shape_ex": str(shapeit_ex),
               "bed": '--input-bed '+str(chrstem)+'.bed '+str(chrstem)+'.bim '+str(chrstem)+'.fam',
               "map": '--input-map '+str(args.ref_maps).replace('###','${chrom}'),
               "ref": str(ph_ref_txt),
               "window": '--window '+str(args.window),
               "duo_txt": str(duo_txt),
               "thread_str": '--thread '+str(args.threads),
               "seed_str": '--seed '+str(args.shape_seed),
               "outmax": '--output-max '+str(outstem)+'.phased.haps '+str(outstem)+'.phased.sample',
               "shapelog": str(outstem)+'.shape.resub_'+str(num_chr)+'.log',
	       "cbopen":'{{',
	       "cbclose":'}}',
               }    
    shape_cmd = cmd_templ.format(**jobdict)

    # submit
    jobres = send_job(jobname='shapeit.'+str(outdot)+'.resub_'+str(num_chr),
	              cmd=shape_cmd,
	              logname='shapeit.'+str(outdot)+'.resub_'+str(num_chr)+'.sub.'+str(clust_conf['log_task_id'])+'.log',
	              mem=int(args.mem_req)*1000,
	              walltime=30,
	              njobs=int(num_chr),
	              threads=args.threads,
		      sleep=args.sleep)

    print 'Pre-phasing jobs re-submitted for %d chromosomes.\n' % num_chr



    # put this job back in the queue
    print '\n...Replacing this imputation job in the queue...'
    
    os.chdir(wd)
    imp_log = 'imp_chunks.'+str(outdot)+'.sub.log'

    send_job(jobname='imp.chunks.'+str(outdot),
             cmd=' '.join(sys.argv[:]),
             logname=imp_log,
             mem=8000,
             walltime=2, # week
             wait_name='shapeit.'+str(outdot)+'.resub_'+str(num_chr),
	     wait_num=str(jobres).strip(),
             sleep=args.sleep)

    print '\n############'
    print '\n'
    print 'All jobs submitted.\n'
    exit(0)




######################
print '\n...Setting up working directories...'
######################

# working directories
chunk_dir = wd + '/chunks_for_imp'
imp_dir = wd + '/imp_sub'
os.mkdir(chunk_dir)
os.mkdir(imp_dir)
print 'Genomic chunks: %s' % chunk_dir
print 'Imputation: %s' % imp_dir



######################
print '\n...Creating genomic chunks for imputation...'
######################

os.chdir(chunk_dir)
link(str(shape_dir)+'/'+str(args.bfile)+'.hg19.ch.fl.bim', str(args.bfile)+'.hg19.ch.fl.bim', 'reference-aligned bim file')

# create chunks
chunk_call = [chunker_ex,
              '--bfile',str(args.bfile)+'.hg19.ch.fl',
              '--out',str(args.out),
              addout_txt[0],addout_txt[1],
              '--Mb-size',str(args.Mb_size),
              '--snp-size',str(args.snp_size),
              '--chr-info-file',str(args.chr_info_file)]
chunk_call = filter(None,chunk_call)

chunk_log = open('chunk.'+str(outdot)+'.log', 'w')
print ' '.join(chunk_call) + '\n'
subprocess.check_call(chunk_call, stderr=subprocess.STDOUT, stdout=chunk_log)
chunk_log.close()



######################
print '\n...Submitting IMPUTE job array...'
# - require "allow_large_regions" to handle centromere boundaries
# - write submit script to include chunk name parsing
# TODO: consider making queue/resources flexible
#
# impute2 
# -use_prephased_g 
# -known_haps_g chr20.phased.haps 
# -h chr20.reference.hap.gz 
# -l chr20.reference.legend.gz 
# -m chr20.gmap.gz 
# -int 35e6 36e6 
# -Ne 20000
# -allow_large_regions
# -o_gz
# -o chr20.imputed
# -seed 12345
#
######################

os.chdir(imp_dir)
link(str(chunk_dir)+'/'+str(outdot)+'.chunks.txt', str(outdot)+'.chunks.txt', 'genomic chunk results')

# job script
imp_templ = dedent("""\
    cchr=`awk -v a={task} 'NR==a+1{cbopen}print $1{cbclose}' {cfile}`
    cstart=`awk -v a={task} 'NR==a+1{cbopen}print $2{cbclose}' {cfile}`
    cend=`awk -v a={task} 'NR==a+1{cbopen}print $3{cbclose}' {cfile}`
    cname=`awk -v a={task} 'NR==a+1{cbopen}print $4{cbclose}' {cfile}`

    {impute_ex} {version_args} {g_arg} {in_haps} -h {ref_haps} -l {ref_leg} -m {map} -int ${cbopen}cstart{cbclose} ${cbopen}cend{cbclose} -buffer {buffer} -Ne {Ne} -o_gz -o {out} {seedtxt}

    {info_cmd}
""")

if args.imp_version==2:
    version_args = '-allow_large_regions'
    g_arg = '-use_prephased_g -known_haps_g'
    info_cmd = ''
elif args.imp_version==4:
    version_args = '-no_maf_align'
    g_arg = '-g'
    info_cmd = qctool_ex + ' -g {out}.gen.gz -snp-stats -osnp {out}.qctool_info.txt -log {out}.qctool_info.log'
    info_cmd = info_cmd.format(**{"out": str(outdot)+'.imp.${{cname}}'})

# fill in template
jobdict = {"task": "{task}",
           "cfile": str(outdot)+'.chunks.txt',
           "impute_ex": str(impute_ex),
           "version_args": str(version_args),
           "g_arg": str(g_arg),
           "in_haps": str(shape_dir)+'/'+str(outdot)+'.chr${{cchr}}.phased.haps',
           "ref_haps": str(args.ref_haps).replace('###','${{cchr}}'),
           "ref_leg": str(args.ref_legs).replace('###','${{cchr}}'),
           "map": str(args.ref_maps).replace('###','${{cchr}}'),
           "Ne": str(args.Ne),
           "buffer": str(args.buffer),
           "out": str(outdot)+'.imp.${{cname}}',
           "seedtxt": str(seedtxt),
           "cbopen":'{{',
           "cbclose":'}}',
           "info_cmd": str(info_cmd)
           }


# get number of chunks (-1 is for header)
nchunks = file_len(outdot+'.chunks.txt')-1


# store job information for possible resubs
job_store_file = 'imp.chunks.'+str(outdot)+'.pkl'

clust_dict = init_sendjob_dict()
clust_dict['jobname'] = 'imp.chunks.'+str(outdot)
clust_dict['logname'] = str('imp.chunks.'+str(outdot)+'.'+str(clust_conf['log_task_id'])+'.sub.log')
clust_dict['mem'] = 8000
clust_dict['walltime'] = 2
clust_dict['njobs'] = int(nchunks)
clust_dict['sleep'] = args.sleep

save_job(jfile=job_store_file, cmd_templ=imp_templ, job_dict=jobdict, sendjob_dict=clust_dict)


# submit
cmd_imp = imp_templ.format(**jobdict)

jobres2 = send_job(jobname='imp.chunks.'+str(outdot),
         	   cmd=cmd_imp,
         	   logname=str('imp.chunks.'+str(outdot)+'.'+str(clust_conf['log_task_id'])+'.sub.log'),
         	   mem=8000,
         	   walltime=2,
         	   njobs=int(nchunks),
         	   sleep=args.sleep)
print 'Imputation jobs submitted for %d chunks.\n' % nchunks



###
# submit next imputation task
###
if args.full_pipe:
    ######################
    print '\n...Queuing best-guess script...'
    ######################
    
    os.chdir(wd)
    next_call = str(rp_bin) + '/bg_imp.py '+' '.join(sys.argv[1:])

    bg_log = 'bg_imp.'+str(outdot)+'.sub.log'

    # TODO: consider queue/mem for agg
    send_job(jobname='bg.chunks.'+str(outdot),
             cmd=next_call,
             logname=bg_log,
             mem=8000,
             walltime=2, # week
             wait_name='imp.chunks.'+str(outdot),
	     wait_num=str(jobres2).strip(),
             sleep=args.sleep)



print '\n############'
print '\n'
print 'SUCCESS!\n'
print 'All jobs submitted.\n'
exit(0)
# eof
