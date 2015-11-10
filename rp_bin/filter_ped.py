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
import random
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
arg_prefWt = parser.add_argument_group('Filtering Preferences')

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

arg_prefWt.add_argument('--case-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of case status when selecting individuals to ' + \
                             'keep among cryptically related pairs',
                        required=False,
                        default=5.0)
arg_prefWt.add_argument('--con-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of control status when selecting individuals to ' + \
                             'keep among cryptically related pairs',
                        required=False,
                        default=2.0)
arg_prefWt.add_argument('--miss-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of missing phenotype when selecting individuals to ' + \
                             'keep among cryptically related pairs',
                        required=False,
                        default=1.0)
arg_prefWt.add_argument('--fam-case-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of a related case in pedigree when selecting ' + \
                             'individuals to keep among cryptically related pairs',
                        required=False,
                        default=5.0)
arg_prefWt.add_argument('--fam-con-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of a related control in pedigree when selecting ' + \
                             'individuals to keep among cryptically related pairs',
                        required=False,
                        default=2.0)
arg_prefWt.add_argument('--fam-miss-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of a related individual in pedigree with a ' + \
                             'missing phenotype when selecting individuals to keep ' + \
                             'among cryptically related pairs',
                        required=False,
                        default=1.0)
arg_prefWt.add_argument('--cross-fid-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of relationship to individuals in a different ' + \
                             'pedigree when selecting individuals to keep among ' + \
                             'cryptically related pairs',
                        required=False,
                        default=-10.0)
arg_prefWt.add_argument('--geno-weight',
                        type=float,
                        metavar='FLOAT',
                        help='Value of genotyping rate when selecting individuals ' + \
                             'to keep among cryptically related pairs. NOT CURRENTLY USED.',
                        required=False,
                        default=0.1)
arg_prefWt.add_argument('--rand-weight',
                        type=float,
                        metavar='FLOAT',
                        help='When selecting individuals to keep among cryptically ' + \
                             'related pairs, the range of the random value used to ' + \
                             'break ties. NOTE: a large value here can override the ' + \
                             'preference from the case/control/pedigree weights.',
                        required=False,
                        default=1e-5)
arg_prefWt.add_argument('--seed',
                        type=int,
                        metavar='INT',
                        help='Random seed. Applies to random values for breaking ties ' + \
                             'among cryptically related individuals',
                        required=False,
                        default=123)

args = parser.parse_args()

# removing tie-breaker would likely cause problems 
assert args.rand_weight != 0.0, 'To prevent unexpected behavior with tied preferences, --rand-weight must be non-zero.'

# make dict of weights for easier use
pref_weights = {
        'case' : float(args.case_weight),
        'control' : float(args.con_weight),
        'missing': float(args.miss_weight),
        'fam_case' : float(args.fam_case_weight),
        'fam_control' : float(args.fam_con_weight),
        'fam_missing' : float(args.fam_miss_weight),
        'rel_cross' : float(args.cross_fid_weight),
        'geno_rate' : float(args.geno_weight),
        'rand_range' : float(args.rand_weight),
}

# set random seed
random.seed(int(args.seed))


# magic numbers for parent/offspring tagging
PO_PI_MAX = 0.6
PO_PI_MIN = 0.4
PO_IBD0_MAX = 0.25
PO_IBD2_MAX = 0.25
PO_IBD1_MIN = 0.75


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
            iid_relatives[jointid1].append(jointid2)
        else:
            iid_relatives[jointid1] = [jointid2]
        
        if jointid2 in iid_relatives:
            iid_relatives[jointid2].append(jointid1)
        else:
            iid_relatives[jointid2] = [jointid1]
        
    ibd_file.close()
            

# should be impossible if --format choices correctly parsed
else:
    raise ValueError('Unsupported format for IBD file: %s' % str(args.format))


#############
print '\n...Flagging possible unmarked parent/offspring pairs...'
#############

def isPossiblePO(pair_info, 
                 PO_PI_MIN=PO_PI_MIN, 
                 PO_PI_MAX=PO_PI_MAX, 
                 PO_IBD0_MAX=PO_IBD0_MAX, 
                 PO_IBD1_MIN=PO_IBD1_MIN, 
                 PO_IBD2_MAX=PO_IBD2_MAX):

    assert 'pihat' in pair_info, 'Invalid argument for pair_info (missing entry \'pihat\')'
    assert 'ibds' in pair_info, 'Invalid argument for pair_info (missing entry \'ibds\')'

    if float(pair_info['ibds'][0]) <= PO_IBD0_MAX and \
       float(pair_info['ibds'][1]) >= PO_IBD1_MIN and \
       float(pair_info['ibds'][2]) <= PO_IBD2_MAX:
           return True
    else:
           return False

# print 'info --- %s' % str(rel_info['QC1025302436:QC1025302436:::QC1025302407:QC1025302407'])
# print 'trying QC1025302436:QC1025302436:::QC1025302407:QC1025302407 --- %s' % str(isPossiblePO(rel_info['QC1025302436:QC1025302436:::QC1025302407:QC1025302407']))

possibleParents = [rel for rel,info in rel_info.iteritems() if isPossiblePO(info)]

def isFamPO(pair_info, fam_info):
    if pair_info['fam1'] == pair_info['fam2']:
        if fam_info[pair_info['joint1']][2] == pair_info['id2']:
            return True # iid2 is father of iid1
        elif fam_info[pair_info['joint1']][3] == pair_info['id2']:
            return True # iid2 is mother of iid1
        elif fam_info[pair_info['joint2']][2] == pair_info['id1']:
            return True # iid1 is father of iid2
        elif fam_info[pair_info['joint2']][3] == pair_info['id1']:
            return True # iid1 is mother of iid2
        else:
            return False
    else:
        return False

# if any...
if possibleParents: 
    
    parents_file = open(str(args.out) + '.unmarkedparents.txt', 'w')
    parents_file.write(' '.join(["FID1","IID1","SEX1","FID2","IID2","SEX2","pihat","ibd0","ibd1","ibd2"]) + '\n')

    n_po = 0
    n_new_po = 0

    for po_pair in possibleParents:
        n_po += 1
        if isFamPO(rel_info[po_pair], fam_info):
            next
        else:
            n_new_po += 1
            po_joint1 = rel_info[po_pair]['joint1']
            po_joint2 = rel_info[po_pair]['joint2']
            parents_file.write(' '.join([
                                        str(fam_info[po_joint1][0]),
                                        str(fam_info[po_joint1][1]),
                                        str(fam_info[po_joint1][4]),
                                        str(fam_info[po_joint2][0]),
                                        str(fam_info[po_joint2][1]),
                                        str(fam_info[po_joint2][4]),
                                        str(rel_info[po_pair]['pihat']),
                                        str(rel_info[po_pair]['ibds'][0]),
                                        str(rel_info[po_pair]['ibds'][1]),
                                        str(rel_info[po_pair]['ibds'][2])
                                    ]) + '\n')
    
    parents_file.close()
    print 'Found %d putative parent-offspring pairs (%d not indicated in fam file).' % (n_po, n_new_po)



#############
print '\n...Flagging IDs to remove for cross-FID cryptic relatedness...'
#############

# get IDs with cross-FID relationships
cross_ids = [(info['joint1'],info['joint2']) for rel,info in rel_info.iteritems() if info['cross_fid'] == True ]

# if any...
if cross_ids:
    
    # get list of unique IDs in cross-FID list
    cross_id_list = list(set(list(sum(cross_ids, ()))))
    
    # log
    print 'Found '+str(len(cross_ids))+' cross-FID relationships involving '+str(len(cross_id_list))+' IDs.'

    # define function to score preference for keeping each individual
    # lowest score will get deleted
    def pref_score(ind_id, fam_dict, rel_dict, weight_dict):
        # init
        pref = 0.0
        ind_id = str(ind_id)
        
        # score self phenotype
        phen = fam_dict[ind_id][5]
        if int(phen) == 2:
            pref += weight_dict['case']
        elif int(phen) == 1:
            pref += weight_dict['control']
        elif int(phen) == 0 or int(phen) == -9:
            pref += weight_dict['missing']
        else:
            raise ValueError('Unexpected phenotype value %r for %s. Not case/control?' % (phen, ind_id))
    
        # score relationships
        for rel_id in rel_dict[ind_id]:
            rel_phen = int(fam_dict[rel_id][5])
            if fam_dict[rel_id][0] != fam_dict[ind_id][0]:
                pref += weight_dict['rel_cross']
            elif rel_phen == 2:
                pref += weight_dict['fam_case']
            elif rel_phen == 1:
                pref += weight_dict['fam_control']
            elif rel_phen == 0 or rel_phen == -9:
                pref += weight_dict['fam_missing']                

        # TODO: score geno rate

        return pref
    
    # loop removal until no cross-fid relationship left
    cryptex_file = open(str(args.out) + '.remove.crossFID.txt', 'w')
    while len(cross_id_list) > 0:
    
        # score each cross-FID related IID's prority for keep/remove
        prefs = [pref_score(indiv, fam_info, iid_relatives, pref_weights) for indiv in cross_id_list]
        
        # breaks ties randomly
        if len(prefs) != len(set(prefs)):
            prefs = [x + random.uniform(-1.0*pref_weights['rand_range'], 0) for x in prefs]
            
        # select individual with lowest preference to remove
        fail_id = cross_id_list[prefs.index(min(prefs))]
        fail_fid = fam_info[fail_id][0]
        fail_iid = fam_info[fail_id][1]
        
        # write selected individual to file
        cryptex_file.write(' '.join([str(fail_fid), str(fail_iid), 'pref='+str(min(prefs)), 'related='+';'.join(iid_relatives[fail_id])]) + '\n')
        
        ### update relationship info for next iteration
        # remove dict entries
        iid_relatives.pop(fail_id)
        fam_info.pop(fail_id)
        fid_members[fail_fid] = fid_members[fail_fid].remove(fail_iid)
        iid_relatives = {key: value.remove(fail_id) if fail_id in value else value for key, value in iid_relatives.iteritems()}    
        rel_info = {key: value for key, value in rel_info.iteritems() if (value['joint1'] != fail_id and value['joint2'] != fail_id)}
    
        # update cross-FID relatedness list
        cross_ids = [(info['joint1'],info['joint2']) for rel,info in rel_info.iteritems() if info['cross_fid'] == True ]
        cross_id_list = list(set(list(sum(cross_ids, ()))))
        
        print 'Flagged %s. Now %d cross-FID related pairs remaining.' % (fail_id, len(cross_ids))

    cryptex_file.close()



#############
print '\n...Flagging IIDs unrelated to reported FIDs...'
#############

nonrelex_file = open(str(args.out) + '.remove.nonFID.txt', 'w')
nonrelex_file.write(' '.join(["FID","IID"]) + '\n')

n_nonrelex = 0

# loop individuals in fam file
for indiv in fam_info:
    
    # print 'indiv = %s' % str(indiv)    
    
    # ok if IID is only member of FID
    if len(fid_members[fam_info[indiv][0]]) < 2:
        # print '    Found no others in FID'
        continue
    
    # if is related to at least 1 person...
    elif indiv in iid_relatives:
        # print '    Found some relatives'

        # check if at least one of them is in same FID
        num_in_fam = 0 
        for id2 in iid_relatives[indiv]:
            if fam_info[id2][0] == fam_info[indiv][0]:
                num_in_fam += 1
                # print '    Found a related family member'
                break
        
        if num_in_fam > 0:
            continue
    
    # if here, are other IIDs in FID but not related to any of them
    nonrelex_file.write(' '.join([str(fam_info[indiv][0]), str(fam_info[indiv][1])]) + '\n') 
    n_nonrelex += 1

nonrelex_file.close()

if n_nonrelex > 0:
    print 'Found %d IIDs not related to other individuals with same FID.' % n_nonrelex
    print 'List of flagged IIDs written to %s.' % str(nonrelex_file.name)




exit(0)
