#! /usr/bin/env python

print '...Importing packages...'
# load requirements
import subprocess
import argparse



# init vars that may be set as functions of others



print '...Parsing arguments...' 
# parse arguments
parser = argparse.ArgumentParser(prog='pca_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))
#...#




print '...Checking dependencies...'



print '...Setting up directory...'



print '...Computing IMUS (unrelated) set...'




print '...Extracting IMUS set from data...'





print '...Computing PCA with IMUS individuals...'





print '...Projecting PCs for remaining individuals...'




print '...Plotting PCs...'




print '############'
print '\n'
print 'SUCCESS!\n'
exit(0)
