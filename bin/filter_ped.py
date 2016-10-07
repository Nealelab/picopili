#! /usr/bin/env python

####################################
# filter_ped.py
# written by Raymond Walters, August 2015
"""
Filter family-based GWAS data based on cryptic relatedness
"""
# Overview:
# 1) Input relatedness estimates
#    - allow REAP, plink (regular or full) formats
# TODO: add king format, allow running new estimates
# 2) Flag possible parent/offspring pairs not already in fam file
# 3) Flag parent/offsping pairs from fam not supported by IBD results
# 4) Flag cross-FID relatedness
# 5) Flag within-FID unrelateds
# 6) Return exclusion lists with reason, flags for followup
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
import argparse
import random
import warnings
from args_ped import parserbase, parsergeno, parseribd, parserweights
from py_helpers import unbuffer_stdout
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='filter_ped.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase,parsergeno,parseribd,parserweights])

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
PO_IBD0_MAX = 0.2
PO_IBD1_MIN = 0.8
PO_IBD2_MAX = 0.2

# print settings
print 'Using settings:'
print '--input-ibd '+str(args.input_ibd)
print '--bfile '+str(args.bfile)
print '--geno '+str(args.geno)
print '--out '+str(args.out)
print '--format '+str(args.format)
print '--min-rel '+str(args.min_rel)
print '--case-weight '+str(args.case_weight)
print '--con-weight '+str(args.con_weight)
print '--miss-weight '+str(args.miss_weight)
print '--fam-case-weight '+str(args.fam_case_weight)
print '--fam-con-weight '+str(args.fam_con_weight)
print '--fam-miss-weight '+str(args.fam_miss_weight)
print '--cross-fid-weight '+str(args.cross_fid_weight)
print '--geno-weight '+str(args.geno_weight)
print '--rand-weight '+str(args.rand_weight)
print '--seed '+str(args.seed)


#############
print '\n...Checking dependencies...'
#############

# verify input files exist
assert os.path.isfile(args.input_ibd), "IBD/relatedness file does not exist (%r)" % args.input_ibd
assert os.path.isfile(str(args.bfile)+'.fam'), "Plink fam file does not exist (%s)" % str(args.bfile)+'.fam'

if str(args.geno) != 'NONE':
    assert os.path.isfile(args.geno), "Missingness rate file does not exist (%r)" % args.geno


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
fam_parents = []

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

        # record parent/offspring pairs
        if str(pat) != '0':
            joint_pat = str(fid) + ':' + str(pat)
            pat_pair = joint_pat + ':::' + joint_id
            fam_parents.append(pat_pair)
            
        if str(mat) != '0':
            joint_mat = str(fid) + ':' + str(mat)
            mat_pair = joint_mat + ':::' + joint_id
            fam_parents.append(mat_pair)        


# TODO: add option to compute IBD on the fly here?



#############
print '\n...Loading genotyping rate information...'
# Assume plink imiss format
# if no file, load dummy values with no missingness
#############

genorate = {}

if str(args.geno) == 'NONE':
    print 'Skipping (no file provided).'    
    for indiv in fam_info:
        genorate[indiv] = 1.0

else:
    genofile = open(str(args.geno), 'r')

    # verify expected header
    head = genofile.readline().split()
    if not head == ['FID','IID','MISS_PHENO','N_MISS','N_GENO','F_MISS']:
        raise IOError('First line of %s does not match expected header format for plink imiss file' % str(args.geno))
    
    # read per individual, indexed by FID:IID
    for line in genofile:
        (fid, iid, miss_phen, nmiss, ngeno, fmiss) = line.split()
        
        # id key
        ind_id = str(fid) + ':' + str(iid)
        
        # record
        genorate[ind_id] = 1.0 - float(fmiss)
    
    genofile.close()
    
    # check values present for all IDs
    for indiv in fam_info:
        if indiv in genorate:
            continue
        else:
            warnings.warn('Genotyping rate not loaded for %s. Setting to zero.' % str(indiv))
            genofile[indiv] = 1.0



#############
print '\n...Parsing relatedness estimates...'
# handle reap file format
# load in values for:
# - dict of relationships for each individual
# - dict of IBD info for each relationship pair
#############

iid_relatives = {}
rel_info = {}
    
# option of taking gzipped files
if str(args.input_ibd).endswith('.gz'):
    import gzip
    ibd_file = gzip.open(str(args.input_ibd), 'rb')
    
else:
    ibd_file = open(str(args.input_ibd), 'r')


    
# read header and verify it matches expected format
head = ibd_file.readline().split()
if str(args.format) == 'reap':
    if not head == ['FID1','IID1','FID2','IID2','N_SNP','IBD0_PROB','IBD1_PROB','IBD2_PROB','KINCOEF']:
        ibd_file.close()
        raise IOError('First line of %s does not match expected header format for REAP' % str(args.input_ibd))

elif str(args.format) == 'plink':
    if not head == ['FID1','IID1','FID2','IID2','RT','EZ','Z0','Z1','Z2','PI_HAT','PHE','DST','PPC','RATIO']:
        ibd_file.close()
        raise IOError('First line of %s does not match expected header format for Plink (require output from \'--genome full\')' % str(args.input_ibd))

elif str(args.format) == 'plink_full':
    if not head == ['FID1','IID1','FID2','IID2','RT','EZ','Z0','Z1','Z2','PI_HAT','PHE','DST','PPC','RATIO','IBS0','IBS1','IBS2','HOMHOM','HETHET']:
        ibd_file.close()
        raise IOError('First line of %s does not match expected header format for Plink \'--genome full\'' % str(args.input_ibd))

else: # should be prevented by argparse
    ibd_file.close()
    raise ValueError('Failed to handle IBD file format specification. Expected \'reap\' or \'plink\', was given \'%s\'.' % str(args.format) )

# process lines
for line in ibd_file:
    # parse, get needed values
    if args.format == 'reap':
        (fid1, iid1, fid2, iid2, nsnp, ibd0, ibd1, ibd2, kin) = line.split()
        pi = 2.0 * float(kin)

    elif args.format == 'plink':
        (fid1, iid1, fid2, iid2, rt, ez, ibd0, ibd1, ibd2, pi, phe, dst, ppc, ratio) = line.split()
    
    elif args.format == 'plink_full':
        (fid1, iid1, fid2, iid2, rt, ez, ibd0, ibd1, ibd2, pi, phe, dst, ppc, ratio, ibs0, ibs1, ibs2, homhom, hethet) = line.split()
        
    # build index variables
    jointid1 = str(fid1) + ':' + str(iid1)
    jointid2 = str(fid2) + ':' + str(iid2)
    pair = jointid1 + ':::' + jointid2
    crossfid = bool (fid1 != fid2)

    # skip relationships that don't hit minimum relatedness threshold
    if float(pi) < float(args.min_rel):
        continue

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



#############
print '\n...Flagging possible unmarked parent/offspring pairs...'
#############

def isPossiblePO(pair_info,
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


# get list of related pairs that are possibly parent/offspring
possibleParents = [rel for rel,info in rel_info.iteritems() if isPossiblePO(info)]

# get starting number of parent pairs
# new = genetic PO not in fam file
# bad = fam file PO not supported by genetic rel.
# mis = parent missing from fam file
n_fam_po = len(fam_parents)
n_gen_po = len(possibleParents)
n_new_po = 0
n_bad_po = 0
n_mis_po = 0

# if any genetic parent/offspring pairs
if possibleParents: 
       
    for po_pair in possibleParents:

        # get reversed version for checking in fam_parents[]   
        pair_joint_ids = po_pair.split(':::')
        po_pair_reverse = str(pair_joint_ids[1] + ':::' + pair_joint_ids[0])

        # if in fam as PO, mark and move on
        if isFamPO(rel_info[po_pair], fam_info):
            if po_pair in fam_parents:
                fam_parents.remove(po_pair)
            else:
                fam_parents.remove(po_pair_reverse)
            continue                
        
        # otherwise is new; flag in file
        else:
            n_new_po += 1
            
            # if first found, initialize output file
            if n_new_po == 1:
                parents_file = open(str(args.out) + '.unmarkedparents.txt', 'w')
                parents_file.write(' '.join(["FID1","IID1","SEX1","FID2","IID2","SEX2","pihat","ibd0","ibd1","ibd2"]) + '\n')
                        
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
    
    # finished flagging genetic PO pairs
    if n_new_po > 0:
        parents_file.close()
    
# log 
print 'Found %d putative parent-offspring pairs from relatedness info.' % int(n_gen_po)
if n_new_po > 0:
    print '%d pairs not in fam file written to %s.' % (int(n_new_po), str(parents_file.name))
elif n_gen_po > 0:
    print 'All putative pairs indicated in fam file.'


#############
print '\n...Flagging possibly incorrect parent/offspring pairs from fam file...'
#############

# any remaining in fam_parents are missing genetic relationship
if fam_parents:
    
    for po_pair in fam_parents:
        # get joint_ids of parent, offspring 
        # (are stored in that order in fam_parents)  
        pair_joint_ids = po_pair.split(':::')
        
        # parent not present in plink fam file
        if pair_joint_ids[0] not in fam_info.keys():
            n_mis_po += 1
            continue
        
        # otherwise parent exists, but relationship is wrong
        else:
            n_bad_po += 1
            
            # if needed, init output file
            if n_bad_po == 1:
                bad_par_file = open(str(args.out) + '.wrongparents.txt', 'w')
                bad_par_file.write(' '.join(["FID1","IID1","SEX1","FID2","IID2","SEX2","pihat","ibd0","ibd1","ibd2"]) + '\n')
           
            # to record
            po_joint1 = pair_joint_ids[0]
            po_joint2 = pair_joint_ids[1]
            
            # get relatedness info if present
            if po_pair in rel_info.keys():
                bad_par_file.write(' '.join([
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
                
            # possible relatedness file has pair in reverse order
            elif str(pair_joint_ids[1] + ':::' + pair_joint_ids[0]) in rel_info.keys():
                po_pair = str(pair_joint_ids[1] + ':::' + pair_joint_ids[0])
                bad_par_file.write(' '.join([
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
        
            # otherwise, pi/ibd=NA since below threshold for relationship file
            else:
                bad_par_file.write(' '.join([
                            str(fam_info[po_joint1][0]),
                            str(fam_info[po_joint1][1]),
                            str(fam_info[po_joint1][4]),
                            str(fam_info[po_joint2][0]),
                            str(fam_info[po_joint2][1]),
                            str(fam_info[po_joint2][4]),
                            str('NA'),
                            str('NA'),
                            str('NA'),
                            str('NA')
                        ]) + '\n')

    if n_bad_po > 0:
        bad_par_file.close()


# log
n_good_po = int(n_fam_po)-len(fam_parents)
print 'Found %d parent-offspring pairs indicated by fam file.' % int(n_fam_po)
if n_mis_po > 0:
    print '%d pairs with parent not present in fam file.' % int(n_mis_po)
print '%d are supported by genetic relatedness.' % (n_good_po)
if n_bad_po > 0:
    print '%d pairs with genetic relatedness inconsistent with parent/offspring relationship.' % int(n_bad_po)
    print 'Bad pairs written to %s.' % str(bad_par_file.name)



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
    def pref_score(ind_id, fam_dict, rel_dict, geno_dict, weight_dict):
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

        # score geno rate
        pref += weight_dict['geno_rate'] * float(geno_dict[ind_id])

        return pref
    
    # loop removal until no cross-fid relationship left
    cryptex_file = open(str(args.out) + '.remove.crossFID.txt', 'w')
    while len(cross_id_list) > 0:
    
        # score each cross-FID related IID's prority for keep/remove
        prefs = [pref_score(indiv, fam_info, iid_relatives, genorate, pref_weights) for indiv in cross_id_list]
        
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
        fid_members[fail_fid].remove(fail_iid)
#        iid_relatives = {key: value.remove(fail_id) if fail_id in value else value for key, value in iid_relatives.iteritems()}
        for value in iid_relatives.values():
            try:
                value.remove(fail_id)
            except ValueError:
                pass
        rel_info = {key: value for key, value in rel_info.iteritems() if (value['joint1'] != fail_id and value['joint2'] != fail_id)}
    
        # update cross-FID relatedness list
        cross_ids = [(info['joint1'],info['joint2']) for rel,info in rel_info.iteritems() if info['cross_fid'] == True ]
        cross_id_list = list(set(list(sum(cross_ids, ()))))
        
        print 'Flagged %s. Now %d cross-FID related pairs remaining.' % (fail_id, len(cross_ids))

    cryptex_file.close()
    print 'Task finished. Results written to %s.' % str(cryptex_file.name) 

else:
    print 'None found.'


#############
print '\n...Flagging IIDs unrelated to reported FIDs...'
#############

n_nonrelex = 0

# loop individuals in fam file
for indiv in fam_info:
    
    # print 'indiv = %s' % str(indiv)
    
    # ok if IID is only member of FID
    if len(fid_members[fam_info[indiv][0]]) < 2:
        # print '    Found no others in FID'
        continue
    
    # if is related to at least 1 person...
    # (latter check possible if all relatives removed in cross-FID filters)
    elif indiv in iid_relatives and iid_relatives[indiv] is not None:
        # print '    Found some relatives'
        # print iid_relatives[indiv]


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
    n_nonrelex += 1

    # only create output file, write header if needed
    if n_nonrelex == 1:
        nonrelex_file = open(str(args.out) + '.remove.nonFID.txt', 'w')
        nonrelex_file.write(' '.join(["FID","IID"]) + '\n')
    
    nonrelex_file.write(' '.join([str(fam_info[indiv][0]), str(fam_info[indiv][1])]) + '\n') 
#     nonrelex_file.write(' '.join([str(fam_info[indiv])]) + '\n') 
    
# close at end and report
if n_nonrelex > 0:
    nonrelex_file.close()

    print 'Found %d IIDs not related to other individuals with same FID.' % n_nonrelex
    print 'List of flagged IIDs written to %s.' % str(nonrelex_file.name)
else:
    print 'None found.'


print '\n\n'
print '############'
print 'Finished!'
print '############'
print '\n'
exit(0)
