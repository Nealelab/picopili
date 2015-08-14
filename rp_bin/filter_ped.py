#! /usr/bin/env python

####################################
# filter_ped.py
# written by Raymond Walters, August 2015
"""
Filter family-based GWAS data based on cryptic relatedness
"""
# Overview:
# 1) Input relatedness estimates
#    - allow REAP formats
# TODO: add plink, king formats, allow running new estimates
# 2) Filter cross-FID relatedness
# 3) Filter within-FID unrelateds
# 4) Return exclusion lists with reason, flags for followup
#
####################################


####################################
# Setup
# a) load python dependencies
# b) get variables/arguments
####################################

import sys
#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Importing packages...'
#############

### load requirements
# import os
# import subprocess
import argparse
# from string import ascii_uppercase
# from glob import glob
# from numpy import digitize
from py_helpers import unbuffer_stdout
# file_len, test_exec, read_conf, find_from_path, link, gz_confirm
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='filter_ped.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))

arg_base = parser.add_argument_group('Basic Arguments')
arg_ibd = parser.add_argument_group('Relatedness File Settings')
arg_prefwt = parser.add_argument_group('Filtering Preferences')

arg_base.add_argument('--input-ibd', 
                      type=str,
                      metavar='FILE',
                      help='file containing IBD relatedness estimates',
                      required=True)           
arg_base.add_argument('--bfile', 
                      type=str,
                      metavar='FILESTEM',
                      help='file stem for plink bed/bim/fam files',
                      required=True)
arg_base.add_argument('--out',
                      type=str,
                      metavar='OUTNAME',
                      help='base name for output files; recommend 4 character stem to match ricopili',
                      required=True)
#arg_base.add_argument('--no-cleanup',
#                      action='store_true',
#                      help='skip cleanup of interim files')

arg_ibd.add_argument('--format',
                      type=str.lower,
                      choices=['reap'],
                      help='format of the input IBD file',
                      required=False,
                      default='reap') 
arg_ibd.add_argument('--min-rel',
                    type=float,
                    metavar='FLOAT',
                    help='Minimum pi-hat relatedness level to treat as related. ' + \
                         'Default is halfway between 3rd and 4th degree relatives.',
                    required=False,
                    default=.09375)

# arg_prefwt.add_argument()

args = parser.parse_args()




# TODO: add dependency checks, args printing, reading config files, etc





print '\n'
print '############'
print 'Begin!'
print '############'

#############
print '\n...Loading fam file...'
# Assume plink fam format
# Load values for:
#  - dict of IIDs in each FID
#  - dict of fam file entry for each FID:IID key
#############

fid_members = {}
fam_info = {}

with open(str(args.bfile) + '.fam', 'r') as fam_file:
    for line in fam_file:
        (fid, iid, pat, mat, sex, phen) = line.split()
        
        joint_id = str(fid) + ':' + str(iid)
        
        # record fid/iid existence
        if str(fid) not in fid_members:
            fid_members[str(fid)] = [str(iid)] 
            
        elif str(iid) not in fid_members[str(fid)]:
            fid_members[str(fid)].append(str(iid))
        # else: already listed
               
        # record full fam entry for id
        if joint_id in fam_info:
            raise ValueError('Duplicated FID:IID value in fam file (%s): %s' % (str(args.bfile)+'.fam', joint_id))
        else:
            fam_info[joint_id] = [fid, iid, pat, mat, sex, phen]



# TODO: add option to compute IBD on the fly here?

#############
print '\n...Parsing relatedness estimates...'
# handle reap file format
# load in values for:
# - dict of relationships for each individual
# - dict of IBD info for each relationship pair
#############

iid_relatives = {}
rel_info = {}

### parse reap files
if str(args.format) == 'reap':
    
    # option of taking gzipped reap files
    if str(args.input_ibd).endswith('.gz'):
        import gzip
        ibd_file = gzip.open(str(args.input_ibd), 'rb')
        
    else:
        ibd_file = open(str(args.input_ibd), 'r')
    
    # read header and verify it matches expected REAP format
    head = ibd_file.readline().split()
    if not head == ['FID1','IID1','FID2','IID2','N_SNP','IBD0_PROB','IBD1_PROB','IBD2_PROB','KINCOEF']:
        raise IOError('First line of %s does not match expected header format for REAP' % str(args.input_ibd))
    
    # process lines
    for line in ibd_file:
        # parse
        (fid1, iid1, fid2, iid2, nsnp, ibd0, ibd1, ibd2, kin) = line.split()
        
        # build needed variables
        pi = 2.0 * float(kin)
        jointid1 = str(fid1) + ':' + str(iid1)
        jointid2 = str(fid2) + ':' + str(iid2)
        pair = jointid1 + ':::' + jointid2
        crossfid = bool (fid1 != fid2)

        # skip relationships that don't hit minimum relatedness threshold
        if pi < float(args.min_rel):
            next

        # record info on the related pair
        # mark cross-fid relatedness
        rel_info[pair] = {'fam1': fid1, 'id1': iid1, 'joint1': jointid1, 
                          'fam2': fid2, 'id2': iid2, 'joint2': jointid2,
                          'ibds': (ibd0,ibd1,ibd2), 'pihat': pi, 'cross_fid': crossfid}
        
        # record relationship for the individuals
        if jointid1 in iid_relatives:
            iid_relatives[jointid1].append(pair)
        else:
            iid_relatives[jointid1] = [pair]
        
        if jointid2 in iid_relatives:
            iid_relatives[jointid2].append(pair)
        else:
            iid_relatives[jointid2] = [pair]
        
    ibd_file.close()
            

# should be impossible if --format choices correctly parsed
else:
    raise ValueError('Unsupported format for IBD file: %s' % str(args.format))


# print cross-fid relatives
print [(info['joint1'],info['joint2']) for rel,info in rel_info.iteritems() if info['cross_fid'] == True ]


exit(0)
