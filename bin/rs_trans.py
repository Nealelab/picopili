#! /usr/bin/env python

####################################
# rs_trans.py
# written by Raymond Walters, January 2016
"""
Update rsIDs and info scores for best-guess calls
"""
# Overview:
# 1) Parse arguments
# 5) Convert chunks to plink files
# 6) Fix bim files (extract rsids, etc)
#    - farm 5 & 6 out to cluster
# 
####################################


###
# Arg parsing
###

import os
import argparse

parser = argparse.ArgumentParser(prog='rs_trans.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))

parser.add_argument('--chunk',
                    type=str,
                    help='Name of genomic chunk',
                    required=True)
parser.add_argument('--name',
                    type=str,
                    help='file name stem',
                    required=True)
parser.add_argument('--info-file',
                    type=str,
                    help='Imputation directory with _info file for chunk',
                    required=True)
parser.add_argument('--fam-trans',
                    type=str,
                    help='Fam ID translation file',
                    required=True)
parser.add_argument('--gt-bim',
                    type=str,
                    help='bim file for genotyped data (to identify which SNPs typed vs imputed)',
                    required=False)
parser.add_argument('--imp-version',
                    type=int,
                    choices=[2,4],
                    help="Version of IMPUTE (for expected info file format). Currently supports version 2 or 4. Version 4 expects qctool snpstats format.",
                    required=False,
                    default=2)

                    
args = parser.parse_args()

# rename args
chname = str(args.chunk)
outdot = str(args.name)


###
# process bim to update names, save translation
###

bim_trans = {}
bim_names = []
bim_raw = open(str(outdot) + '.bg.filtered.' + str(chname) + '.bim', 'r')
bim_out = open(str(outdot) + '.bg.filtered.' + str(chname) + '.bim_rsids', 'w')
for line in bim_raw:
    (chrom, snpid, cm, bp, a1, a2) = line.split()

    ####
    # logic for parsing SNP names
    # (1) Assume snpid is one of
    # - rsid
    # - rsid:bp:a2:a1
    # - chr:bp:a2:a1
    # - chr:bp:a1:<struct>:bp2
    # (strip leading "chr" on first field if present)
    #
    # (2) If first field is rsid, take it
    # else take first 2 (chr:bp)
    # 
    # (3) Check if extracted name is unique.
    # If not, revert to full name in bim
    #
    # (4) If name still isn't unique (a rare/
    # unexpected edge case), append m's until it is.
    # - Hack inspired by ricopili
        
    snp_split = snpid.split(':')

    if str(snp_split[0][:3]).lower() == "chr":
        snp_split[0] = snp_split[0][3:]
        
    if str(snp_split[0]) == str(chrom):
        rs = str(snp_split[0]) + ':' + str(snp_split[1])
    else:
        rs = str(snp_split[0])
    
    if rs in bim_names:
        if snpid in bim_names:
            # unexpected/rare edge case
            rs = str(snpid) + 'm'
            while rs in bim_names:
                rs = str(snpid) + 'm'
        else:
            rs = str(snpid)
    
    # write bim file with updated name
    bim_names.append(rs)
    bim_out.write(' '.join([str(chrom),str(rs),str(cm),str(bp),str(a1),str(a2)]) + '\n')

    # save translation for creating info file
    bim_trans[str(snpid)] = [str(rs), str(chrom), str(cm), str(bp), str(a1), str(a2)]

bim_raw.close()
bim_out.close()


###
# Get genotype status
###

gt_bim = open(os.path.realpath(str(args.gt_bim)), 'r')

gt_snps = []
for line in gt_bim:
    (chrom, snpid, cm, bp, a1, a2) = line.split()
    gt_snps.append(snpid)

gt_bim.close()


###
# Get info scores
###

# info_in = open(imp_dir+'/'+outdot+'.imp.'+str(chname)+'_info', 'r')
info_in = open(os.path.realpath(str(args.info_file)), 'r')
info_out = open(str(outdot)+'.bg.filtered.'+str(chname)+'.info', 'w')

if args.imp_version==2:
    dumphead = info_in.readline()
elif args.imp_version==4:
    dumphead = info_in.readline()
    while '#' in dumphead:
        dumphead = info_in.readline()

info_out.write(' '.join(['SNP','oldID','CHR','BP','A1','A2','INFO','EXP_A1','TYPED'])+'\n')
for line in info_in:
    
    if args.imp_version==2:
        (foo, snpid, bpinf, exp_frq, info, certain, imptype, info0, concord0, r2_0) = line.split()
    elif args.imp_version==4:
        
        # catch completion message at end of file
        if '#' in line:
            continue

        # assumes QC tool snpstats will provide:
        # alternate_ids rsid chromosome position alleleA alleleB comment HW_exact_p_value HW_lrt_p_value alleleA_count alleleB_count alleleA_frequency alleleB_frequency minor_allele_frequency minor_allele major_allele info impute_info missing_proportion A B AA AB BB NULL total
        (aid, snpid, chrinf, bpinf, aA, aB, com, hwe, hwe2, Acount, Bcount, Afreq, exp_frq, maf, minorA, majorA, qcinfo, info, miss, A, B, AA, AB, BB, null, total) = line.split()

    # check if SNP survived filtering        
    if str(snpid) not in bim_trans.keys():
        continue
    
    # if yes, record info
    else:
        # lookup used name, etc
        rs = str(bim_trans[str(snpid)][0])
        chrom = str(bim_trans[str(snpid)][1])
        bp = str(bim_trans[str(snpid)][3])
        a1 = str(bim_trans[str(snpid)][4])
        a2 = str(bim_trans[str(snpid)][5])
        
        # get genotyped status
        if snpid in gt_snps:
            gt = 1
        else:
            gt = 0

        # write to file
        info_out.write(' '.join([rs, snpid, chrom, bp, a1, a2, str(info), str(exp_frq), str(gt)]) + '\n')
        
        # shrink dict
        bim_trans.pop(str(snpid))
        
info_in.close()
info_out.close()



###
# translate fam file
###

# load translation file
fam_dict = {}
fam_trans_file = open(str(args.fam_trans), 'r')

for line in fam_trans_file:
    (tmpfid, tmpiid, realfid, realiid) = line.split()
    
    tmp_id = str(tmpfid) + ':' + str(tmpiid) 

    fam_dict[tmp_id] = [tmpfid, tmpiid, realfid, realiid]

fam_trans_file.close()



# process best-guess fam file
fam_in = open(str(outdot) + '.bg.filtered.' + str(chname) + '.fam', 'r')
fam_out = open(str(outdot) + '.bg.filtered.' + str(chname) + '.fam_trans', 'w')

for line in fam_in:
    (fid, iid, pat, mat, sex, phen) = line.split()
    
    # dict key
    joint = str(fid) + ':' + str(iid)

    # lookup full ID
    # note: assuming id will be in trans file (fails otherwise)
    orig_fid = str(fam_dict[joint][2])
    orig_iid = str(fam_dict[joint][3])
    
    # check for parents
    if str(mat) is not '0':
        joint_mat = str(fid) + ':' + str(mat)
        orig_mat = str(fam_dict[joint_mat][3])
    else:
        orig_mat = '0'        
    
    if str(pat) is not '0':
        joint_pat = str(fid) + ':' + str(pat)
        orig_pat = str(fam_dict[joint_pat][3])
    else:
        orig_pat = '0'
    
    # write translated file
    fam_out.write(' '.join([orig_fid, orig_iid, orig_pat, orig_mat, str(sex), str(phen)]) + '\n')

fam_in.close()
fam_out.close()




# finish
print '\n############'
print '\n'
print 'SUCCESS!'
exit(0)
# eof