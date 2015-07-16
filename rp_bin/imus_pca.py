#! /usr/bin/env python

####################################
# imus_pca.py
# written by Raymond Walters, July 2015
"""
Runs PCA for GWAS data with related individuals
"""
# Overview:
# 1) Input QCed plink bed/bim/fam
# 2) Define set of unrelated individuals using PRIMUS
# 3) Compute PCA on the unrelated set and projects to remainder
# 4) Plot projected PCs
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
from py_helpers import file_len, read_conf


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

### init vars that may be set as functions of others
pcadir = ""

### parse arguments
parser = argparse.ArgumentParser(prog='imus_pca.py',
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
                    help='name for PCA output directory, defaults to OUTNAME_imus_pca',
                    required=False)
parser.add_argument('--no_cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')
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
#parser.add_argument('--smartpca-ex',
#                    type=str,
#                    metavar='PATH',
#                    help='path to smartpca executable',
#                    required=False,
#                    default="/humgen/atgu1/fs03/shared_resources/shared_software/EIG6.0beta_noreq/bin/smartpca")

args = parser.parse_args()

# set remaining defaults
if args.pcadir == None:
    pcadir = args.out + '_imus_pca'
else:
    pcadir = args.pcadir
    
wd = os.getcwd()

# print settings
print 'Using settings:'
print '--bfile '+args.bfile
print '--out '+args.out
print '--rel-deg '+str(args.rel_deg)
print '--npcs '+str(args.npcs)


 
#############
print '\n...Reading ricopili config file...'
#############

### read plink loc from config
# not getting R here since ricopili.conf currently relies on platform info
conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

plinkx = configs['p2loc']+"plink"
smartpcax = configs['eloc']+"/smartpca"


#############
print '\n...Checking dependencies...'
# check exists, executable
#############

# R from path
if args.rscript_ex == None:
    raise AssertionError('Unable to find Rscript in search path')
assert os.path.isfile(args.rscript_ex), "Rscript not found at %r" % args.rscript_ex
assert os.access(args.rscript_ex, os.X_OK), "Rscript not executable (%r)" % args.rscript_ex
print "Rscript found: %s" % args.rscript_ex

# primus
assert os.path.isfile(args.primus_ex), "PRIMUS not found at %r" % args.primus_ex
assert os.access(args.primus_ex, os.X_OK), "PRIMUS not executable (%r)" % args.primus_ex
print "PRIMUS found: %s" % args.primus_ex

# plink
assert os.path.isfile(plinkx), "Plink not found at %r" % plinkx
assert os.access(plinkx, os.X_OK), "Plink not executable (%r)" % plinkx
print "Plink found: %s" % plinkx
    
# smartpca
assert os.path.isfile(smartpcax), "Eigensoft smartpca not found at %r" % smartpcax
assert os.access(smartpcax, os.X_OK), "Eigensoft smartpca not executable (%r)" % smartpcax
print "Smartpca found: %s" % smartpcax



print '\n'
print '############'
print 'Begin!'
print '############'

####################################
# Compute maximum unrelated set
# a) run PRIMUS
# b) verify ran successfully
####################################

#############
print '\n...Computing IMUS (unrelated) set...'
#############

primelog = 'primus_' + args.out + '_imus.log'
subprocess.check_call([args.primus_ex,
                       "--file", args.bfile,
                       "--genome",
                       "--degree_rel_cutoff", str(args.rel_deg),
                       "--no_PR",
                       "--plink_ex", plinkx,
                       "--smartpca_ex", smartpcax,
                       "&>", primelog])

# verify successful output
primedir = os.getcwd() + '/' + args.bfile + '_PRIMUS'
imus_file = args.bfile + '_cleaned.genome_maximum_independent_set'
imus_dirfile = primedir + '/' + imus_file

if not os.path.isdir(primedir):
    raise IOError("Expected PRIMUS output directory %r not found" % primedir)
elif not os.path.isfile(imus_dirfile):
    raise IOError("Failed to create IMUS set (missing %r)" % imus_dirfile)



####################################
# Start output directory
# a) create PCA directory with links to input files, PRIMUS results
# b) verify links
# c) extr
####################################

#############
print '\n...Setting up PCA directory...'
#############

if not os.path.exists(pcadir):
    os.makedirs(pcadir)

os.chdir(pcadir)

# setup file links
os.symlink(imus_dirfile,imus_file)
os.symlink(wd+'/'+args.bfile+'.bed',args.bfile+'.bed')
os.symlink(wd+'/'+args.bfile+'.bim',args.bfile+'.bim')
os.symlink(wd+'/'+args.bfile+'.fam',args.bfile+'.fam')

# verify links
if not os.path.isfile(imus_file):
    raise IOError("Failed to link IMUS file (%r)" % imus_file)
elif not os.path.isfile(args.bfile+'.bed'):
    raise IOError("Failed to link bed file (%r)" % str(args.bfile+'.bed') )
elif not os.path.isfile(args.bfile+'.bim'):
    raise IOError("Failed to link bim file (%r)" % str(args.bfile+'.bim') )
elif not os.path.isfile(imus_file):
    raise IOError("Failed to link fam file (%r)" % str(args.bfile+'.fam') )


#############
print '\n...Extracting IMUS set from data...'
#############

bfile_imus = args.bfile + '.imus'
subprocess.check_call([plinkx,
                       "--bfile", args.bfile,
                       "--keep", imus_file,
                       "--freq",
                       "--silent",
                       "--make-bed",
                       "--allow-no-sex",
                       "--out", bfile_imus])



####################################
# Compute PCA using unrelated set and project to full sample
# a) create pedind file labelling IMUS vs. RELATEDS for projection 
# b) setup smartpca par file
# c) run smartpca to compute PCs on IMUS, project to remainder
# d) process output
####################################

# load IDs from imus bim file, in format FID:IID
imus_ids = []
with open(str(bfile_imus+'.fam'), 'r') as f:
    imus_ids = [':'.join(line.split()[0:2]) for line in f]


### process fam file from full data
# - read fam file
# - create short id, print conversion to file
# - compare FID/IID to IMUS set
# - print pedind file for smartpca (shortfid, shortiid, 0, 0, U, [IMUS or RELATEDS])

# init conversion file
id_conv = open(str(args.bfile)+'.pedind.ids.txt', 'w')
id_conv.write('FID IID pca_id' + '\n')

# init pedind file
pedind = open(str(args.bfile)+'.pedind', 'w')

# process fam file by line
bfile_fam = open(str(args.bfile+'.fam'), 'r')
i=0
for line in bfile_fam:
    # iterate line counter, used for short id
    i += 1
    
    # read
    (longfid, longiid, pat, mat, sex, phen) = line.split()

    # assign and record short identifier
    shortfid = str(i)
    shortiid = str(i)
    pcaid = shortfid +':'+ shortiid
    id_conv.write(longfid +' '+ longiid +' '+ pcaid + '\n')

    # get FID:IID identifier to compare to imus
    bfile_id = longfid +':'+ longiid
    
    # write pedind with "IMUS" if match, "RELATED" if not
    if any(bfile_id == refid for refid in imus_ids):
        pedind.write(shortfid +' '+ shortiid +' 0 0 U IMUS')
    else:
        pedind.write(shortfid +' '+ shortiid +' 0 0 U RELATEDS')
    

pedind.close()
id_conv.close()


### create par file
par = open(str(bfile_imus + '.pca.par'), 'w')

par.write('genotypename:     '+str(args.bfile+'.bed')+'\n')
par.write('snpname:          '+str(args.bfile+'.bim')+'\n')
par.write('indivname:        '+str(args.bfile+'.pedind')+'\n')
par.write('poplistname:      '+str(args.bfile+'.refpoplist.txt')+'\n')
par.write('evecoutname:      '+str(args.bfile+'.pca.txt')+'\n')
par.write('snpweightoutname  '+str(args.bfile+'.pca_snpw.txt')+'\n')
par.write('fastmode:         '+'YES'+'\n')
par.write('altnormstyle:     '+'NO'+'\n')
par.write('numoutevec:       '+str(args.npcs)+'\n')
par.write('numoutlieriter:   '+str(0)+'\n')

par.close()


### create poplist file
poplist = open(str(args.bfile+'.refpoplist.txt'), 'w')
poplist.write("IMUS")
poplist.close()


### run smartpca
subprocess.check_call([smartpcax, 
                       "-p", str(bfile_imus + '.pca.par') ])

exit(0)

### TODO: process output
# - convert fid:iids back
# - fix header

####################################
# Plot projected PCs
# a) TODO
####################################

#############
print '\n...Plotting PCs...'
#############

#########
######### TODO: add r plotting
#########


######### TODO: no cm



####################################
# Clean up files
# TODO
####################################

os.chdir(wd)

if not args.no_cleanup:
    
    #############
    print '\n...Clean-up interim files...'
    #############

# zip 1:
# [str(args.out + '.projpca.pc' + str(i) + '.log') for i in xrange(1,args.npcs+1)]

# zip 2:
# str(pcadir + '/' + args.out+'_imus_pca.evec.txt')
# str(pcadir + '/' + args.out+'_imus_pca.eval.txt')
# str(pcadir + '/' + args.out+'_imus_pca.pve.txt')
    
    #############
    print 'Compress:'
    #############
    subprocess.check_call(["gzip", "-fv", str(pcadir + '/' + pc_out_nam)])
    
    #############
    print 'Remove interim:'
    #############
    subprocess.check_call(["rm", "-v",
                           str(pcadir + '/' + pc_files_nam),
                           str(pcadir + '/' + snpw)])    
    
    #############
    print 'Remove if exist:'
    #############
    # allowing failure, since files may or may not exists
    subprocess.call(["rm", "-v",
                     str(pcadir + '/' + bfile_imus +'.nosex'),
                     str(pcadir + '/' + bfile_imus +'.hh'),
                     [str(pcadir + '/' + args.out + '.projpca.pc' + str(i) + '.nosex') for i in xrange(1,args.npcs+1)],
                     [str(pcadir + '/' + args.out + '.projpca.pc' + str(i) + '.hh') for i in xrange(1,args.npcs+1)]])


    

print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)
