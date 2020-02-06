#! /usr/bin/env python

####################################
# impute_rel_chr23.py
# written by Raymond Walters, January 2016
# edited by Nikolas Baya for chrX/chr23 imputation, February 2020
"""
Impute GWAS data with related individuals
"""
# 
# Is a wrapper for:
# - shape_rel_chr23.py
# - imp2_rel_chr23.py
# - bg_imp_chr23.py
# - agg_imp_chr23.py
# 
####################################



import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
import os
import argparse

from args_impute import parserbase, parserphase, parserimpute, parserchunk, parserref, arg_ref, parserbg, parsercluster
from py_helpers import unbuffer_stdout
from blueprint import send_job
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

arg_ref.add_argument('--ref-dir',
                         type=str,
                         metavar='DIRECTORY',
                         help='Directory containing imputation reference files (haps, legends, sample, and maps). ' + 
                              'Used as prefix for specifying full paths of --ref-maps, --ref-haps, --ref-legs, and --ref-samps',
                         required=False,
                         default=None)

parser = argparse.ArgumentParser(prog='impute_rel_chr23.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parserphase, parserimpute, parserchunk, parserref, parserbg, parsercluster])
                    
args = parser.parse_args()

if args.ref_dir is not None:
    # verify exists
    assert os.path.isdir(args.ref_dir), "Failed to find imputation reference directory %s" % args.ref_dir
    
    # prepend to references accordingly    
    args.ref_maps = str(args.ref_dir) +'/' + args.ref_maps
    args.ref_haps = str(args.ref_dir) +'/' + args.ref_haps
    args.ref_legs = str(args.ref_dir) +'/' + args.ref_legs
    args.ref_samps = str(args.ref_dir) +'/' + args.ref_samps

# reference recommendation
def print_ref_rec():
    print '\nIf you do not have an imputation reference available, the 1000 Genomes'
    print 'Phase 3 reference panel provided by IMPUTE is directly compatible with'
    print 'picopili and broadly covers most major continental populations.'
    print '\nDirect download:'
    print 'wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz'
    print '\nWARNING: download filesize is > 12 GB\n'    

# check these references exist
if not os.path.isfile(args.ref_maps.replace('###','23')):
    print "Failed to verify genetic maps exist."
    print_ref_rec()
    raise IOError("No chr 23 genetic map: %s" % args.ref_maps.replace('###','23'))
    
if not os.path.isfile(args.ref_haps.replace('###','23')):
    print "Failed to verify reference haplotypes exist."
    # print rec, since is possible have genetic map but not imputation panel
    print_ref_rec()
    raise IOError("No chr 23 reference haplotypes: %s" % args.ref_haps.replace('###','23'))
    
if not os.path.isfile(args.ref_legs.replace('###','23')):
    # not printing ref_rec here since at this point have verified haplotypes exist
    raise IOError("Failed to verify reference legend files exist (tested for chr 23 at %s)" % args.ref_legs.replace('###','23'))
    
if not os.path.isfile(args.ref_samps.replace('###','23')):
    # not printing ref_rec here since at this point have verified haplotypes exist
    raise IOError("Failed to verify reference sample file exists (tested for chr 23 at %s)" % args.ref_samps.replace('###','23'))


# more flexible handling for info file for shapeit, since could be external
if not os.path.isfile(args.ref_info.replace('###','23')):
        
        if args.ref_dir is not None and os.path.isfile(str(args.ref_dir) +'/' + args.ref_info.replace('###','23')):
            args.ref_info = str(args.ref_dir) +'/' + args.ref_info
            
        else:
            print "Reference information file for phasing not found (tested for chr 23: %s)." % args.ref_info.replace('###','23')
            if args.ref_dir is not None:
                print "Tried both relative path and in --ref-dir %s" % str(args.ref_dir)
            
            if args.ref_dir == "1000GP_Phase3_chr###.legend.gz":
                print "For 1000 Genomes Phase 3 reference from IMPUTE the required file is "
                print "the same as the reference legend.\n"
                print "Maybe you wanted to add this\n?"
                # verified above that the legend file exists
                print "--ref-info %s\n" % args.ref_legs
            
            raise IOError("Failed to verify phasing info file exists (tested for chr 23 at %s)" % args.ref_info.replace('###','23'))


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
print '--ref-dir '+str(args.ref_dir)
print '--ref-maps '+str(args.ref_maps)
print '--ref-haps '+str(args.ref_haps)
print '--ref-legs '+str(args.ref_legs)
print '--ref-samps '+str(args.ref_samps)


print '\nBest-guess genotype calling:'
if args.hard_call_th is None:
    print '--bg-th '+str(args.bg_th)
else:
    print '--hard-call-th '+str(args.hard_call_th)
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
next_call = str(rp_bin) + '/shape_rel_chr23.py '+' '.join(sys.argv[1:])+' --full-pipe'

shape_log = 'shape.'+str(outdot)+'.sub.log'

# TODO: consider queue/mem
send_job(jobname='shape.'+str(outdot),
         cmd=next_call,
         logname=shape_log,
         mem=int(args.mem_req * 1000),
         walltime=30,
         sleep=args.sleep)


# finish
print '\n############'
print '\n'
print 'SUCCESS!'
print 'All jobs submitted.\n'
exit(0)


# eof
