#! /usr/bin/env python

####################################
# ped_confirm.py
# written by Raymond Walters, November 2015
"""
Confirm reported pedigrees consistent with relatedness-based reconstruction 
"""
# Overview:
# 1) Run PRIMUS to estimate pedigrees from relatedness
# 2) Compare estimated to reported pedigrees
#       - split fam file by FID
#       - add dummy parents if needed
#       - run find_expected_pedigree.pl from PRIMUS
# 3) Parse logs to find possible mismatches
# 4) Ouput results
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
import os
import subprocess
import argparse
import re
import warnings
from args_ped import parserbase, parseribd, parserexloc
from py_helpers import unbuffer_stdout, test_exec
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='ped_confirm.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parseribd,parserexloc])

args = parser.parse_args()

wd = os.getcwd()

# print settings
print 'Using settings:'
print '--input-ibd '+str(args.input_ibd)
print '--bfile '+str(args.bfile)
print '--out '+str(args.out)
print '--format '+str(args.format)
print '--min-rel '+str(args.min_rel)
print '--max-gens '+str(args.max_gens)
print ' '

# verify input files exist
assert os.path.isfile(args.input_ibd), "IBD/relatedness file does not exist (%r)" % args.input_ibd
assert os.path.isfile(str(args.bfile)+'.fam'), "Plink fam file does not exist (%s)" % str(args.bfile)+'.fam'

# test executables
test_exec(args.primus_ex, 'PRIMUS')
test_exec(args.findped_ex, 'PRIMUS pedigree matching script')
print ' '

# unzip relatedness file if needed
if args.input_ibd.endswith('.gz'):
    ibd_txtfile = str(args.input_ibd) + '.txt'
    print 'Unzipping IBD relatedness file to %s' % ibd_txtfile
    ibd_out = open(ibd_txtfile, 'w')
    subprocess.check_call(['gunzip','-c',str(args.input_ibd)],stdout=ibd_out)
    ibd_out.close()
else:
    ibd_txtfile = str(args.input_ibd)

assert os.path.isfile(ibd_txtfile), "Failed to extract IBD/relatedness file (%r)" % args.input_ibd



print '\n'
print '############'
print 'Begin!'
print '############'


#############
print '\n...Constructing pedigrees from IBD/relatedness...'
# - Run using primus
# - args depend on ibd file format (default here: reap)
#############

primus_peds_log = open(str(args.out)+'.primus_peds.log', 'w')

if args.format=='reap':
    pr_input_text = '--input FILE='+str(ibd_txtfile)+' IBD0=6 IBD1=7 IBD2=8 RELATEDNESS=9'
elif args.format=='plink' or args.format=='plink_full':
    pr_input_text = '--input FILE='+str(ibd_txtfile)+' IBD0=7 IBD1=8 IBD2=9 PI_HAT=10'
else:
    raise ValueError('Unsupported IBD file format (%s)' % str(args.format))

famfile_text = 'FILE='+str(args.bfile)+'.fam'

subprocess.check_call([str(args.primus_ex),
                       str(pr_input_text),
                       '--rel_threshold',str(args.min_rel),
                       '--no_IMUS',
                       '--max_gens',str(args.max_gens),
                       '--sexes',famfile_text,'SEX=5',
                       '--affections',famfile_text,'AFFECTION=6',
                       '--output_dir',str(args.out)+'_primus_peds'],
                       stderr=subprocess.STDOUT,
                       stdout=primus_peds_log)

primus_peds_log.close()

print 'IBD-based pedigrees written to %s' % ('./'+str(args.out)+'_primus_peds/')

#############
print '\n...Preparing fam files...'
# - load plink fam file
# -- dict of IIDs in each FID
# -- dict of fam file entry for each FID:IID key
# - Format fam files for PRIMUS pedigree matching
# -- Requires separate fam file per FID
# -- Sibs must have parents
# TODO: find more elegant solution for adding dummy parents?
# - currently: add father/mother if >1 IID in FID and no IIDs have listed parents, assume are siblings
#############

fid_members = {}
fam_info = {}

# load in fam file
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


# make folder for all the fam-based pedigrees
fids_to_check = []
fam_dir = str(args.out)+'_fam_peds'
if not os.path.exists(fam_dir):
    os.makedirs(fam_dir)

# loop FIDs to write individual fam files
for fid in fid_members:
    num_indivs = len(fid_members[fid])

    # no check needed for single individuals
    if num_indivs < 2:
        continue
    
    
    else:
        fids_to_check.append(str(fid))
        
        # setup output file
        fam_out = open('./'+str(fam_dir)+'/'+str(args.out)+'.fid_'+str(fid)+'.fam','w')
        
        # check if dummy parents needed. Criteria:
        # - 2+ IIDs in fam
        # - no IIDs with listed parent (e.g. likely siblings)
        hasParents = False
        for indiv in fid_members[fid]:
            ind_id = str(fid) + ':' + str(indiv)
            paternal_id = str(fam_info[ind_id][2])
            maternal_id = str(fam_info[ind_id][3])
            if paternal_id != "0" or maternal_id != "0":
                hasParents = True
                break

        # loop individuals to write, add dummy parents if needed
        for indiv in fid_members[fid]:
            ind_id = str(fid) + ':' + str(indiv)
            if hasParents:
                fam_out.write(' '.join(fam_info[ind_id]) + '\n')
            else:
                dum_pat = str(fid)+'_father'
                dum_mat = str(fid)+'_mother'
                fam_out.write(' '.join([
                                        str(fid),
                                        str(indiv),
                                        dum_pat,
                                        dum_mat,
                                        fam_info[ind_id][4],
                                        fam_info[ind_id][5]
                                ]) + '\n')
            
        fam_out.close()

print 'Fam files for %d families written to %s' % (len(fids_to_check), fam_dir)


#############
print '\n...Comparing IBD pedigrees to fam file pedigrees...'
# using find_expected_pedigree.pl from PRIMUS
# loop FIDs, parse output to find matches vs. non-matches
#############

comp_dir = str(args.out)+'_ped_comparisons'
if not os.path.exists(comp_dir):
    os.makedirs(comp_dir)

pr_file = './'+str(args.out)+'_primus_peds/Summary_'+str(ibd_txtfile)+'.txt'

results_file = open(str(args.out)+'.pedigree_matches.txt', 'w')
results_file.write(' '.join(["FID","hasMatch","network"]) + '\n')

solo_fids = []
match_fids = []
nonmatch_fids = []
nonparse_fids = []

# loop all fids for completness of results file
for fid in fid_members:
    
    if fid not in fids_to_check:
        results_file.write(' '.join([str(fid), 'Solo', 'NA']) + '\n')
        solo_fids.append(str(fid))
        continue

    # run pedigree comparison    
    fid_fam = './'+str(fam_dir)+'/'+str(args.out)+'.fid_'+str(fid)+'.fam'
    fid_log = './'+comp_dir+'/fid_'+str(fid)+'_findped.log'
    fid_log_file = open(fid_log, 'w')
    subprocess.check_call([args.findped_ex,
                           pr_file,
                           fid_fam],
                           stderr=subprocess.STDOUT,
                           stdout=fid_log_file)
    fid_log_file.close()

    
    # parse results
    # key is "MATCH!!!" or "NO MATCHING PEDIGREE FOUND" in output log
    # if has match, record network number
    fidParsed = False
    fid_results = open(fid_log, 'r')
    for line in fid_results:
        if re.match("MATCH!!!(.*)", line):
            linetxt = line.split()
            network = linetxt[1][7:] # entry is e.g. "network27", just want 27
            results_file.write(' '.join([str(fid), 'True', str(network)]) + '\n')
            match_fids.append(str(fid))
            fidParsed = True
            break
        elif re.match("NO MATCHING PEDIGREE FOUND(.*)", line):
            results_file.write(' '.join([str(fid), 'False', 'NA']) + '\n')
            nonmatch_fids.append(str(fid))
            fidParsed = True
            break

    fid_results.close()
    
    if not fidParsed:
        results_file.write(' '.join([str(fid), 'Error', 'Error']) + '\n')
        nonparse_fids.append(str(fid))
        warnings.warn("Failed to parse results for %s (%s). Check for errors." % (str(fid), str(fid_log)))
        

results_file.close()

# summarize results for logging
print 'Found %d reported pedigrees matching relatedness structure' % len(match_fids)
print 'Found %d reported pedigrees NOT matching relatedness structure' % len(nonmatch_fids)
print 'Skipped %d reported pedigrees with only one individual' % len(solo_fids)
print 'Failed to parse results for %d reported pedigrees' % len(nonparse_fids)
print ' '
print 'Results written to: %s' % str(results_file.name)

print '\n\n'
print '############'
print 'Finished!'
print '############'
print '\n'
exit(0)
