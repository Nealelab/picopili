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


print '...Parsing arguments...' 
# parse arguments
parser = argparse.ArgumentParser(prog='pca_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))
#...#
                                 
# plink_ex
# primus_ex
# rscript_ex


# plink file stem
# name
# pcadir
# npca




#set remaining defaults#
if flashpcax == "":
    flashpcax = "/humgen/atgu1/fs03/shared_resources/shared_software/bin/flashpca"
    
if primusx == "":
    primusx = "~/PRIMUS_v1.8.0/bin/run_PRIMUS.pl"
    
if rscriptx =="":
    rscriptx = spawn.find_executable("Rscript") # check success later
 
 

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
assert os.path.isfile(primusx), 'PRIMUS not found at ' + primusx
assert os.path.isfile(flashpcax), 'FlashPCA not found at ' + flashpcax
assert os.path.isfile(plinkx), 'Plink not found at ' + plinkx
assert os.path.isfile(rscriptx), 'Rscript not found at ' + rscriptx

# file executable
assert os.access(primusx, os.X_OK), 'FlashPCA not executable (' + primusx + ')'
assert os.access(flashpcax, os.X_OK), 'FlashPCA not executable (' + flashpcax + ')'
assert os.access(plinkx, os.X_OK), 'Plink not executable (' + plinkx + ')'
assert os.access(rscriptx, os.X_OK), 'Rscript not executable (' + rscriptx + ')'


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
