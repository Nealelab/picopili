#! /usr/bin/env python

print '...Importing packages...'
# load requirements
import os
import subprocess
from distutils import spawn
import argparse



# init vars that may be set as functions of others
primusx = ""
flashpcax = ""
rscriptx = ""
pcadir = ""

print '...Parsing arguments...' 
# parse arguments
parser = argparse.ArgumentParser(prog='pca_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))
parser.add_argument('--bfile', 
                    type=str,
                    metavar='FILESTEM',
                    help='file stem for input plink bed/bim/fam',
                    required=True)
parser.add_argument('--out',
                    type=str,
                    metavar='OUTNAME',
                    help='base name for output; recommend 4 character stem to match ricopili',
                    required=True)
parser.add_argument('--rel-deg',
                    type=int,
                    metavar='INT',
                    help='relatedness degree threshold for defining \"unrelated\" set',
                    required=False,
                    default=3)
parser.add_argument('--npcs',
                    type=int,
                    metavar='INT',
                    help='number of principal components to compute',
                    required=False,
                    default=10)
parser.add_argument('--pcadir',
                    type=str,
                    metavar='DIRNAME',
                    help='name for PCA output directory, defaults to pca_imus_OUTNAME',
                    required=False)
# parser.add_argument('--plink-ex',
#                    type=str,
#                    metavar='PATH',
#                    help='path to plink executable, read from ~/ricopili.conf if unspecified',
#                    required=False)
parser.add_argument('--rscript-ex',
                    type=str,
                    metavar='PATH',
                    help='path to Rscript executable, tries reading from PATH if unspecified',
                    required=False,
                    default=spawn.find_executable("Rscript"))
parser.add_argument('--primus-ex',
                    type=str,
                    metavar='PATH',
                    help='path to PRIMUS executable',
                    required=False,
                    default=os.environ['HOME']+"/PRIMUS_v1.8.0/bin/run_PRIMUS.pl")
parser.add_argument('--flashpca-ex',
                    type=str,
                    metavar='PATH',
                    help='path to flashpca executable',
                    required=False,
                    default="/humgen/atgu1/fs03/shared_resources/shared_software/bin/flashpca")

args, pass_through_args = parser.parse_known_args()

#set remaining defaults#
if pcadir == "":
    pcadir = 'pca_imus_' + args.out
# mem?

 

print '...reading ricopili config file...'
### read plink loc from config

conf_file = os.environ['HOME']+"/ricopili.conf"

configs = {}
with open(conf_file, 'r') as f:
    for line in f:
        (key, val) = line.split()
        configs[str(key)] = val

plinkx = configs['p2loc']+"plink"



print '...Checking dependencies...'
# R from path
if rscriptx == None:
    raise AssertionError('Unable to find Rscript in search path')
    
# file exists
assert os.path.isfile(args.primus_ex), "PRIMUS not found at %r" % args.primus_ex
assert os.path.isfile(args.flashpca_ex), "FlashPCA not found at %r" % args.flashpca_ex
assert os.path.isfile(plinkx), "Plink not found at %r" % plinkx
assert os.path.isfile(args.rscript_ex), "Rscript not found at %r" % args.rscript_ex

# file executable
assert os.access(args.primus_ex, os.X_OK), 'FlashPCA not executable (' + args.primus_ex + ')'
assert os.access(args.flashpca_ex, os.X_OK), 'FlashPCA not executable (' + args.flashpca_ex + ')'
assert os.access(plinkx, os.X_OK), 'Plink not executable (' + plinkx + ')'
assert os.access(args.rscript_ex, os.X_OK), 'Rscript not executable (' + args.rscript_ex + ')'



print '############'
print 'Begin!'
print '############\n'

print '...Computing IMUS (unrelated) set...'



print '...Setting up PCA directory...'




print '...Extracting IMUS set from data...'





print '...Computing PCA with IMUS individuals...'





print '...Projecting PCs for remaining individuals...'




print '...Plotting PCs...'




print '############'
print '\n'
print 'SUCCESS!\n'
exit(0)
