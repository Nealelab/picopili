#! /usr/bin/env python

####################################
# admix_rel.py
# written by Raymond Walters, July 2015
"""
Estimates relatedness for admixed sample using ADMIXTURE, REAP
"""
# Overview:
# 1) Input source and target plink files
#    - assume matching QCed LD pruned SNPs, source is unrelated subset of target IDs
# 2) Run ADMIXTURE on unrelated set
# 3) Select population exemplars based on admixture proportions
# 4) Run ADMIXTURE on full data in supervised mode using pop exemplars
# 5) Estimate relatedness with REAP
# 6) Generate diagnostic plots
#    - exemplars on PCA (if PCA available)
#    - final admixture on PCA (if PCA available)
#    - IBD0/IBD1 for REAP
#
####################################


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

### load requirements
import os
import subprocess
from distutils import spawn
import argparse
from string import ascii_uppercase
from py_helpers import read_conf, file_len, unbuffer_stdout
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='admix_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))

arg_base = parser.add_argument_group('Basic Arguments')
arg_admix = parser.add_argument_group('Admixture Settings')
arg_reap = parser.add_argument_group('Relatedness Settings')
arg_plot = parser.add_argument_group('Plot Settings')
arg_exloc = parser.add_argument_group('Software Executable Locations')

arg_base.add_argument('--unrel-bfile', 
                    type=str,
                    metavar='FILESTEM',
                    help='File stem for plink bed/bim/fam files ' + \
                         'with unrelated individuals to estimate admixture.',
                    required=True)
arg_base.add_argument('--target-bfile', 
                    type=str,
                    metavar='FILESTEM',
                    help='file stem for plink bed/bim/fam files. ' + \
                         'Relatedness will be estimated for these samples. ' + \
                         'All individuals from --unrel-bfile should also be present in this data.',
                    required=True)
arg_base.add_argument('--out',
                    type=str,
                    metavar='OUTNAME',
                    help='base name for output files; recommend 4 character stem to match ricopili',
                    required=True)
arg_base.add_argument('--outdir',
                    type=str,
                    metavar='DIRNAME',
                    help='Directory for output files. Will create if needed. ' + \
                         'Uses ./OUTNAME_admix_rel if unspecified',
                    required=False)
arg_base.add_argument('--no-cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')

arg_admix.add_argument('--npops',
                    type=int,
                    metavar='INT',
                    help='Number of ancestral populations for admixture',
                    required=False,
                    default=4)
arg_admix.add_argument('--prop-th',
                    type=float,
                    metavar='FLOAT',
                    help='Minimum admixture proportion to select an individual ' + \
                         'from the unrelated set as an exemplar for the ancestral population',
                    required=False,
                    default=.9)
arg_admix.add_argument('--min-exemplar',
                    type=int,
                    metavar='INT',
                    help='Minimum number of individuals from unrelated set ' + \
                         'assigned as exemplars for each ancestral population',
                    required=False,
                    default=20)
arg_admix.add_argument('--multithread-cores',
                    type=int,
                    metavar='INT',
                    help='Number of cores to use for multi-threading in admixture analysis',
                    required=False,
                    default=1)

arg_reap.add_argument('--min-rel',
                    type=float,
                    metavar='FLOAT',
                    help='Minimum pi-hat relatedness level to include in output. ' + \
                         'Default is halfway between 3rd and 4th degree relatives.',
                    required=False,
                    default=.09375)
                    
arg_plot.add_argument('--plot-admix-pca',
                    type=str,
                    metavar='FILE',
                    help='PCA file for the target data; used for plotting admixture results. ' + \
                         'If no file given, will skip plotting.', 
                    required=False,
                    default=None)

arg_exloc.add_argument('--rscript-ex',
                    type=str,
                    metavar='PATH',
                    help='path to Rscript executable, tries reading from PATH if unspecified',
                    required=False,
                    default=None)
arg_exloc.add_argument('--admixture-ex',
                    type=str,
                    metavar='PATH',
                    help='path to ADMIXTURE executable',
                    required=False,
                    default="/humgen/atgu1/fs03/shared_resources/shared_software/bin/admixture")
arg_exloc.add_argument('--reap-ex',
                    type=str,
                    metavar='PATH',
                    help='path to REAP executable',
                    required=False,
                    default="/humgen/atgu1/fs03/shared_resources/shared_software/bin/REAP")

#parser.add_argument('--test-sub',
#                    action='store_true',
#                    help='Test run without submitting tasks',
#                    required=False)

args = parser.parse_args()

# set dependent defaults
if args.outdir == None or args.outdir == "None":
    args.outdir = str(args.out)+'_admix_rel'


# print settings
print 'Using settings:'
print '--unrel-bfile '+args.unrel_bfile
print '--target-bfile '+args.target_bfile
print '--out '+args.out
print '--outdir '+args.outdir
print '--npops '+str(args.npops)
print '--prop-th '+str(args.prop_th)
print '--min-exemplar '+str(args.min_exemplar)
print '--min-rel '+str(args.min_rel)
print '--plot-admix-pca '+str(args.plot_admix_pca)



#############
print '\n...Reading ricopili config file...'
#############

### read plink loc from config
# not getting R here since ricopili.conf currently relies on platform info
conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

plinkx = configs['p2loc']+"plink"


#############
print '\n...Checking dependencies...'
# check exists, executable
#############

# get R from path
if args.rscript_ex == None or args.rscript_ex == "None":
    args.rscript_ex = spawn.find_executable("Rscript")
# if still not found
if args.rscript_ex == None or args.rscript_ex == "None":
    raise AssertionError('Unable to find Rscript in search path')
assert os.path.isfile(args.rscript_ex), "Rscript not found at %r" % args.rscript_ex
assert os.access(args.rscript_ex, os.X_OK), "Rscript not executable (%r)" % args.rscript_ex
print "Rscript found: %s" % args.rscript_ex

# PCA plot script, if needed
if not (args.plot_admix_pca == None or args.plot_admix_pca == "None"):
    Rplotpcax = spawn.find_executable("plot_pca.Rscript")
    if Rplotpcax == None:
        raise AssertionError('Unable to find plot_pca.Rscript in search path')
    print "PCA plotting script found: %s" % Rplotpcax

# plink
assert os.path.isfile(plinkx), "Plink not found at %r" % plinkx
assert os.access(plinkx, os.X_OK), "Plink not executable (%r)" % plinkx
print "Plink found: %s" % plinkx

# admixture
assert os.path.isfile(args.admixture_ex), "ADMIXTURE not found at %r" % args.admixture_ex
assert os.access(args.admixture_ex, os.X_OK), "ADMIXTURE not executable (%r)" % args.admixture_ex
print "ADMIXTURE found: %s" % args.admixture_ex

# reap
assert os.path.isfile(args.reap_ex), "REAP not found at %r" % args.reap_ex
assert os.access(args.reap_ex, os.X_OK), "REAP not executable (%r)" % args.reap_ex
print "REAP found: %s" % args.admixture_ex

# pca file
if not (args.plot_admix_pca==None or args.plot_admix_pca=="None"):
    assert os.path.isfile(args.plot_admix_pca), "PCA file does not exist (%r)" % args.plot_admix_pca
    assert '/' not in args.target_bfile, "--plot-admix-pca must specify only a file, not a path"

# verify bfiles are files, not paths
assert '/' not in args.unrel_bfile, "--unrel-bfile must specify only a file stem, not a path"
assert '/' not in args.target_bfile, "--target-bfile must specify only a file stem, not a path"



print '\n'
print '############'
print 'Begin!'
print '############'

#############
print '\n...Setting up working directory (%s)...' % str(args.outdir)
#############

wd = os.getcwd()

if not os.path.isdir(str(args.outdir)):
    os.makedirs(str(args.outdir))

os.chdir(args.outdir)

# link ref plink files
os.symlink(str(wd+'/'+args.unrel_bfile+'.bed'), str(args.unrel_bfile+'.bed'))
os.symlink(str(wd+'/'+args.unrel_bfile+'.bim'), str(args.unrel_bfile+'.bim'))
os.symlink(str(wd+'/'+args.unrel_bfile+'.fam'), str(args.unrel_bfile+'.fam'))

# link target plink files
os.symlink(str(wd+'/'+args.target_bfile+'.bed'), str(args.target_bfile+'.bed'))
os.symlink(str(wd+'/'+args.target_bfile+'.bim'), str(args.target_bfile+'.bim'))
os.symlink(str(wd+'/'+args.target_bfile+'.fam'), str(args.target_bfile+'.fam'))

# verify links
if not os.path.isfile(str(args.unrel_bfile+'.bed')):
    raise IOError("Failed to link bed file with unrelated individuals (%r)" % str(args.unrel_bfile+'.bed'))
if not os.path.isfile(str(args.unrel_bfile+'.bim')):
    raise IOError("Failed to link bim file with unrelated individuals (%r)" % str(args.unrel_bfile+'.bim'))
if not os.path.isfile(str(args.unrel_bfile+'.bed')):
    raise IOError("Failed to link fam file with unrelated individuals (%r)" % str(args.unrel_bfile+'.fam'))
if not os.path.isfile(str(args.target_bfile+'.bed')):
    raise IOError("Failed to link bed file with target individuals (%r)" % str(args.target_bfile+'.bed'))
if not os.path.isfile(str(args.target_bfile+'.bim')):
    raise IOError("Failed to link bim file with target individuals (%r)" % str(args.target_bfile+'.bim'))
if not os.path.isfile(str(args.target_bfile+'.bed')):
    raise IOError("Failed to link fam file with target individuals (%r)" % str(args.target_bfile+'.fam'))

# link pca file, if provided
if not (args.plot_admix_pca==None or args.plot_admix_pca=="None"):

    os.symlink(str(wd+'/'+args.plot_admix_pca), str(args.plot_admix_pca))
    
    if not os.path.isfile(str(args.plot_admix_pca)):
        raise IOError("Failed to link PCA file (%r)" % str(args.plot_admix_pca))



#############
print '\n...Running Admixture on unrelated dataset...'
#############

admix_call = [args.admixture_ex,
              str(args.unrel_bfile+'.bed'),
              str(args.npops),
              '-j'+str(args.multithread_cores)]
admix_unrel_log = open(str('admix_'+args.out+'_unrel.log'), 'w')

print str(' '.join(admix_call))
subprocess.check_call(admix_call, stdout=admix_unrel_log)

admix_unrel_log.close()



#############
print '\n...Selecting exemplars for each ancestral population...'
#############
# - identify population assignment (including "-") for each input individual
# - confirm whether there are enough IDs assigned to each populations
# - match population assignments to 

# label for populations are popA, popB, popC, ...
popnames = [str('pop'+ascii_uppercase[i]) for i in range(args.npops)]

# define function returning popname or '-' based on largest proportion
# Note: ties broken in favor of first pop listed in names (possible if th <= 0.5)
def maxpop(props, names, th):
    whichmax = props.index(max(props))
    if props[whichmax] > th:
        outpop = names[whichmax]
    else:
        outpop = '-'
    return outpop

# get list of selected pop for each individual in admixture results
ind_pops = []
admix_pops_file = str(args.unrel_bfile+'.'+str(args.npops)+'.Q')
with open(admix_pops_file, 'r') as f:
    # map() required to read probs as float instead of string
    ind_pops = [maxpop(props=map(float,line.split()), names=popnames, th=args.prop_th) for line in f]

# sanity check parsing
nfam = file_len(str(args.unrel_bfile+'.fam'))
if len(ind_pops) != nfam:
    raise ValueError('Number of individuals parsed from admixture results (%d in %s) ' + \
                     'and fam file of unrelateds (%d in %s) do not match.' % (len(ind_pops), admix_pops_file, int(nfam), str(args.unrel_bfile+'.fam')))

# check have sufficient exemplars
popcounts = [ind_pops.count(popnames[i]) for i in range(args.npops)]
lackingpops = [popcounts[i] < args.min_exemplar for i in range(args.npops)]

print 'Exemplars per population:'
for i in range(args.npops):
    print str(popnames[i] + ': ' + str(popcounts[i]))
print 'Unassigned: '+str(ind_pops.count('-'))

if any(lackingpops):
    print '\n###########\n'
    print 'ERROR: One or more populations with insufficient number of exemplars (<'+str(args.min_exemplar)+').'
    print '\nConsider rerunning with fewer ancestral populations (here: '+str(args.npops)+'), \n' + \
          'a looser threshold for selecting population exemplars (here: '+str(args.prop_th)+'), \n' + \
          'or fewer required exemplars per ancestral population in the unrelated set ' + \
          '(here :'+str(args.min_exemplar)+').\n'
    exit(1)

# match exemplar pop status with FID/IIDs




# run supervised admix

# prep files for reap
# - tped
# ./plink --file mydata --recode12 --output-missing-genotype 0 --transpose --out newfile

# run reap using admix results
# use -m to reduce output
# see sec. 6 of reap documentation


# Generate diagnostic plots
# - exemplars on PCA (if PCA available)
# - final admixture on PCA (if PCA available)
# - IBD0/IBD1 for REAP

# final cleanup

# finish
