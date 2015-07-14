#! /usr/bin/env python

print '...Importing packages...'
# load requirements
import os
import subprocess
from distutils import spawn
import argparse



# init vars that may be set as functions of others
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
parser.add_argument('--smartpca-ex',
                    type=str,
                    metavar='PATH',
                    help='path to smartpca executable',
                    required=False,
                    default="/humgen/atgu1/fs03/shared_resources/shared_software/EIG6.0beta_noreq/bin/smartpca")

args, pass_through_args = parser.parse_known_args()

#set remaining defaults
if args.pcadir == None:
    pcadir = 'pca_imus_' + args.out
else:
    pcadir = args.pcadir
    
# mem?

 

print '...reading ricopili config file...'
### read plink loc from config
# not getting R here since ricopili.conf currently relies on platform info

conf_file = os.environ['HOME']+"/ricopili.conf"

configs = {}
with open(conf_file, 'r') as f:
    for line in f:
        (key, val) = line.split()
        configs[str(key)] = val

plinkx = configs['p2loc']+"plink"



print '...Checking dependencies...'
# R from path
if args.rscript_ex == None:
    raise AssertionError('Unable to find Rscript in search path')
    
# file exists
assert os.path.isfile(args.primus_ex), "PRIMUS not found at %r" % args.primus_ex
assert os.path.isfile(args.flashpca_ex), "FlashPCA not found at %r" % args.flashpca_ex
assert os.path.isfile(args.smartpca_ex), "SmartPCA not found at %r" % args.smartpca_ex
assert os.path.isfile(plinkx), "Plink not found at %r" % plinkx
assert os.path.isfile(args.rscript_ex), "Rscript not found at %r" % args.rscript_ex

# file executable
assert os.access(args.primus_ex, os.X_OK), "FlashPCA not executable (%r)" % args.primus_ex
assert os.access(args.flashpca_ex, os.X_OK), "FlashPCA not executable (%r)" % args.flashpca_ex
assert os.access(args.smartpca_ex, os.X_OK), "FlashPCA not executable (%r)" % args.smartpca_ex
assert os.access(plinkx, os.X_OK), "Plink not executable (%r)" % plinkx
assert os.access(args.rscript_ex, os.X_OK), "Rscript not executable (%r)" % args.rscript_ex



print '############'
print 'Begin!'
print '############\n'

print '...Computing IMUS (unrelated) set...'
primelog = 'primus_' + args.out + '_imus.log'
subprocess.check_call([args.primus_ex,
                       "--file", args.bfile,
                       "--genome",
                       "--degree_rel_cutoff", str(args.rel_deg),
                       "--no_PR",
                       "--plink_ex", plinkx,
                       "--smartpca_ex", args.smartpca_ex,
                       "&>", primelog])

# verify successful output
primedir = os.getcwd() + '/' + args.bfile + '_PRIMUS'
imus_file = primedir + '/' + args.bfile + '_cleaned.genome_maximum_independent_set'

if not os.path.isdir(primedir):
    raise IOError("Expected PRIMUS output directory %r not found" % primedir)
elif not os.path.isfile(imus_file):
    raise IOError("Failed to create IMUS set (missing %r)" % imus_file)



print '...Setting up PCA directory...'




print '...Extracting IMUS set from data...'





print '...Computing PCA with IMUS individuals...'





print '...Projecting PCs for remaining individuals...'




print '...Plotting PCs...'




print '############'
print '\n'
print 'SUCCESS!\n'
exit(0)
