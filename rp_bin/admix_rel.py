#! /usr/bin/env python

####################################
# admix_rel.py
# written by Raymond Walters, July 2015
"""
Estimates relatedness for admixed sample using ADMIXTURE, REAP
"""
# Overview:
# 1) Input source and target plink files
#    - assume matching QCed LD pruned SNPs, source is unrelated subset of target IDs
# 2) Run ADMIXTURE on unrelated set
# 3) Select population exemplars based on admixture proportions
# 4) Run ADMIXTURE on full data in supervised mode using pop exemplars
# 5) Estimate relatedness with REAP
# 6) Generate diagnostic plots
#    - exemplars on PCA (if PCA available)
#    - final admixture on PCA (if PCA available)
#    - IBD0/IBD1 for REAP
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
from string import ascii_uppercase
from py_helpers import read_conf, test_exec, file_len, unbuffer_stdout
unbuffer_stdout()


#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

parser = argparse.ArgumentParser(prog='admix_rel.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40))

arg_base = parser.add_argument_group('Basic Arguments')
arg_admix = parser.add_argument_group('Admixture Settings')
arg_reap = parser.add_argument_group('Relatedness Settings')
arg_plot = parser.add_argument_group('Plot Settings')
arg_exloc = parser.add_argument_group('Software Executable Locations')

arg_base.add_argument('--unrel-bfile', 
                    type=str,
                    metavar='FILESTEM',
                    help='File stem for plink bed/bim/fam files ' + \
                         'with unrelated individuals to estimate admixture.',
                    required=True)
arg_base.add_argument('--target-bfile', 
                    type=str,
                    metavar='FILESTEM',
                    help='file stem for plink bed/bim/fam files. ' + \
                         'Relatedness will be estimated for these samples. ' + \
                         'All individuals from --unrel-bfile should also be present in this data.',
                    required=True)
arg_base.add_argument('--out',
                    type=str,
                    metavar='OUTNAME',
                    help='base name for output files; recommend 4 character stem to match ricopili',
                    required=True)
arg_base.add_argument('--outdir',
                    type=str,
                    metavar='DIRNAME',
                    help='Directory for output files. Will create if needed. ' + \
                         'Uses ./OUTNAME_admix_rel if unspecified',
                    required=False,
                    default=None)
arg_base.add_argument('--no-cleanup',
                    action='store_true',
                    help='skip cleanup of interim files')

arg_admix.add_argument('--npops',
                    type=int,
                    metavar='INT',
                    help='Number of ancestral populations for admixture',
                    required=False,
                    default=4)
arg_admix.add_argument('--prop-th',
                    type=float,
                    metavar='FLOAT',
                    help='Minimum admixture proportion to select an individual ' + \
                         'from the unrelated set as an exemplar for the ancestral population',
                    required=False,
                    default=.9)
arg_admix.add_argument('--min-exemplar',
                    type=int,
                    metavar='INT',
                    help='Minimum number of individuals from unrelated set ' + \
                         'assigned as exemplars for each ancestral population',
                    required=False,
                    default=20)
arg_admix.add_argument('--multithread-cores',
                    type=int,
                    metavar='INT',
                    help='Number of cores to use for multi-threading in admixture analysis',
                    required=False,
                    default=1)

arg_reap.add_argument('--min-rel',
                    type=float,
                    metavar='FLOAT',
                    help='Minimum pi-hat relatedness level to include in output. ' + \
                         'Default is halfway between 3rd and 4th degree relatives.',
                    required=False,
                    default=.09375)
                    
arg_plot.add_argument('--plot-admix-pca',
                    type=str,
                    metavar='FILE',
                    help='PCA file for the target data; used for plotting admixture results. ' + \
                         'If no file given, will skip plotting.', 
                    required=False,
                    default=None)

arg_exloc.add_argument('--rscript-ex',
                    type=str,
                    metavar='PATH',
                    help='path to Rscript executable, tries reading from PATH if unspecified',
                    required=False,
                    default=None)
arg_exloc.add_argument('--admixture-ex',
                    type=str,
                    metavar='PATH',
                    help='path to ADMIXTURE executable',
                    required=False,
                    default="/humgen/atgu1/fs03/shared_resources/shared_software/bin/admixture")
arg_exloc.add_argument('--reap-ex',
                    type=str,
                    metavar='PATH',
                    help='path to REAP executable',
                    required=False,
                    default="/humgen/atgu1/fs03/shared_resources/shared_software/bin/REAP")

#parser.add_argument('--test-sub',
#                    action='store_true',
#                    help='Test run without submitting tasks',
#                    required=False)

args = parser.parse_args()

# set dependent defaults
if args.outdir == None or args.outdir == "None":
    args.outdir = str(args.out)+'_admix_rel'


# print settings
print 'Using settings:'
print '--unrel-bfile '+args.unrel_bfile
print '--target-bfile '+args.target_bfile
print '--out '+args.out
print '--outdir '+args.outdir
print '--npops '+str(args.npops)
print '--prop-th '+str(args.prop_th)
print '--min-exemplar '+str(args.min_exemplar)
print '--min-rel '+str(args.min_rel)
print '--plot-admix-pca '+str(args.plot_admix_pca)



#############
print '\n...Reading ricopili config file...'
#############

### read plink loc from config
# not getting R here since ricopili.conf currently relies on platform info
conf_file = os.environ['HOME']+"/ricopili.conf"
configs = read_conf(conf_file)

plinkx = configs['p2loc']+"plink"


#############
print '\n...Checking dependencies...'
# check exists, executable
#############

# get R from path
if args.rscript_ex == None or args.rscript_ex == "None":
    args.rscript_ex = spawn.find_executable("Rscript")
# if still not found
if args.rscript_ex == None or args.rscript_ex == "None":
    raise AssertionError('Unable to find Rscript in search path')
# verify exec
test_exec(args.rscript_ex, 'Rscript')

# IBD plot script
Rplotibdx = spawn.find_executable("plot_reap_ibd.Rscript")
if(Rplotibdx) == None:
    raise AssertionError('Unable to find plot_reap_ibd.Rscript in search path')
print "IBD plotting script found: %s" % Rplotibdx

# PCA plot script, if needed
if not (args.plot_admix_pca == None or args.plot_admix_pca == "None"):
    Rplotpcax = spawn.find_executable("plot_pca.Rscript")
    if Rplotpcax == None:
        raise AssertionError('Unable to find plot_pca.Rscript in search path')
    print "PCA plotting script found: %s" % Rplotpcax

# plink
test_exec(plinkx, 'Plink')

# admixture
test_exec(args.admixture_ex, 'ADMIXTURE')

# reap
test_exec(args.reap_ex, 'REAP')

# pca file
if not (args.plot_admix_pca==None or args.plot_admix_pca=="None"):
    assert os.path.isfile(args.plot_admix_pca), "PCA file does not exist (%r)" % args.plot_admix_pca
    assert '/' not in args.target_bfile, "--plot-admix-pca must specify only a file, not a path"

# verify bfiles are files, not paths
assert '/' not in args.unrel_bfile, "--unrel-bfile must specify only a file stem, not a path"
assert '/' not in args.target_bfile, "--target-bfile must specify only a file stem, not a path"



print '\n'
print '############'
print 'Begin!'
print '############'

#############
print '\n...Setting up working directory (%s)...' % str(args.outdir)
#############

wd = os.getcwd()

if not os.path.isdir(str(args.outdir)):
    os.makedirs(str(args.outdir))

os.chdir(args.outdir)

# link ref plink files
os.symlink(str(wd+'/'+args.unrel_bfile+'.bed'), str(args.unrel_bfile+'.bed'))
os.symlink(str(wd+'/'+args.unrel_bfile+'.bim'), str(args.unrel_bfile+'.bim'))
os.symlink(str(wd+'/'+args.unrel_bfile+'.fam'), str(args.unrel_bfile+'.fam'))

# link target plink files
os.symlink(str(wd+'/'+args.target_bfile+'.bed'), str(args.target_bfile+'.bed'))
os.symlink(str(wd+'/'+args.target_bfile+'.bim'), str(args.target_bfile+'.bim'))
os.symlink(str(wd+'/'+args.target_bfile+'.fam'), str(args.target_bfile+'.fam'))

# verify links
if not os.path.isfile(str(args.unrel_bfile+'.bed')):
    raise IOError("Failed to link bed file with unrelated individuals (%r)" % str(args.unrel_bfile+'.bed'))
if not os.path.isfile(str(args.unrel_bfile+'.bim')):
    raise IOError("Failed to link bim file with unrelated individuals (%r)" % str(args.unrel_bfile+'.bim'))
if not os.path.isfile(str(args.unrel_bfile+'.bed')):
    raise IOError("Failed to link fam file with unrelated individuals (%r)" % str(args.unrel_bfile+'.fam'))
if not os.path.isfile(str(args.target_bfile+'.bed')):
    raise IOError("Failed to link bed file with target individuals (%r)" % str(args.target_bfile+'.bed'))
if not os.path.isfile(str(args.target_bfile+'.bim')):
    raise IOError("Failed to link bim file with target individuals (%r)" % str(args.target_bfile+'.bim'))
if not os.path.isfile(str(args.target_bfile+'.bed')):
    raise IOError("Failed to link fam file with target individuals (%r)" % str(args.target_bfile+'.fam'))

# link pca file, if provided
if not (args.plot_admix_pca==None or args.plot_admix_pca=="None"):

    os.symlink(str(wd+'/'+args.plot_admix_pca), str(args.plot_admix_pca))
    
    if not os.path.isfile(str(args.plot_admix_pca)):
        raise IOError("Failed to link PCA file (%r)" % str(args.plot_admix_pca))



#############
print '\n...Running Admixture on unrelated dataset...'
#############

admix_call = [args.admixture_ex,
              str(args.unrel_bfile+'.bed'),
              str(args.npops),
              '-j'+str(args.multithread_cores)]
admix_unrel_log = open(str('admix_'+args.out+'_unrel.log'), 'w')

print str(' '.join(admix_call))
print 'Logging to ' + admix_unrel_log.name + '\n'
subprocess.check_call(admix_call, stdout=admix_unrel_log)

admix_unrel_log.close()



#############
print '\n...Selecting exemplars for each ancestral population...'
#############
# - identify population assignment (including "-") for each input individual
# - confirm whether there are enough IDs assigned to each populations
# - match population assignments to FID/IIDs
# - write .pops file for target bfile, .pops.info file 

# label for populations are popA, popB, popC, ...
popnames = [str('pop'+ascii_uppercase[i]) for i in range(args.npops)]

# define function returning popname or '-' based on largest proportion
# Note: ties broken in favor of first pop listed in names (possible if th <= 0.5)
def maxpop(props, names, th):
    whichmax = props.index(max(props))
    if props[whichmax] > th:
        outpop = names[whichmax]
    else:
        outpop = '-'
    return outpop

# get list of selected pop for each individual in admixture results
ind_pops = []
admix_pops_file = str(args.unrel_bfile+'.'+str(args.npops)+'.Q')
with open(admix_pops_file, 'r') as f:
    # map() required to read probs as float instead of string
    ind_pops = [maxpop(props=map(float,line.split()), names=popnames, th=args.prop_th) for line in f]

# sanity check parsing
nfam = file_len(str(args.unrel_bfile+'.fam'))
if len(ind_pops) != nfam:
    raise ValueError('Number of individuals parsed from admixture results (%d in %s) ' + \
                     'and fam file of unrelateds (%d in %s) do not match.' % (len(ind_pops), admix_pops_file, int(nfam), str(args.unrel_bfile+'.fam')))

# check have sufficient exemplars
popcounts = [ind_pops.count(popnames[i]) for i in range(args.npops)]
lackingpops = [popcounts[i] < args.min_exemplar for i in range(args.npops)]

print 'Exemplars per population:'
for i in range(args.npops):
    print str(popnames[i] + ': ' + str(popcounts[i]))
print 'Unassigned: '+str(ind_pops.count('-'))

if any(lackingpops):
    print '\n###########\n'
    print 'ERROR: One or more populations with insufficient number of exemplars (<'+str(args.min_exemplar)+').'
    print '\nConsider rerunning with fewer ancestral populations (here: '+str(args.npops)+'), \n' + \
          'a looser threshold for selecting population exemplars (here: '+str(args.prop_th)+'), \n' + \
          'or fewer required exemplars per ancestral population in the unrelated set ' + \
          '(here :'+str(args.min_exemplar)+').\n'
    exit(1)


### match exemplar pop status with FID/IIDs, record in dict
pop_dict = {}

# process fam file by line
ref_fam = open(str(args.unrel_bfile+'.fam'), 'r')
idnum=0
for line in ref_fam:
    # iterate line counter, used to get elements from ind_pops[]
    idnum += 1
    
    # read
    (fid, iid, pat, mat, sex, phen) = line.split()

    # use FID:IID identifier as key to record pop status
    bfile_id = fid +':'+ iid
    pop_dict[bfile_id] = ind_pops[idnum-1]

ref_fam.close()


### create pop file to match target fam file, pop info file
target_fam = open(str(args.target_bfile+'.fam'), 'r')
target_pop = open(str(args.target_bfile+'.pop'), 'w')
target_popinfo = open(str(args.target_bfile+'.pop.info'), 'w')

for line in target_fam:
    
    # read
    (targetfid, targetiid, pat, mat, sex, phen) = line.split()
    target_id = targetfid +':'+ targetiid
    
    # check dict
    if target_id in pop_dict:
        target_pop.write(pop_dict[target_id] + '\n')
        target_popinfo.write(targetfid + ' ' + targetiid + ' ' + target_id + ' unrel ' + pop_dict[target_id] + '\n')
    else:
        target_pop.write('-' + '\n')
        target_popinfo.write(targetfid + ' ' + targetiid + ' ' + target_id + ' target ' + '-' + '\n')


target_fam.close()
target_pop.close()
target_popinfo.close()



#############
print '\n...Running supervised admixture analysis in target data...'
#############

admix_super_call = [args.admixture_ex,
                    str(args.target_bfile+'.bed'),
                    str(args.npops),
                    '-j'+str(args.multithread_cores),
                    '--supervised']
admix_target_log = open(str('admix_'+args.out+'_target.log'), 'w')

print str(' '.join(admix_super_call))
print 'Logging to ' + admix_target_log.name + '\n'
subprocess.check_call(admix_super_call, stdout=admix_target_log)

admix_target_log.close()



#############
print '\n...Preparing admixture results for relatedness analysis...'
#############

### get tped genotypes of target dataset
reap_tped = str(args.target_bfile + '.tmp_recode')
subprocess.check_call([plinkx,
                       '--silent',
                       '--bfile', str(args.target_bfile),
                       '--recode12',
                       '--output-missing-genotype', '0',
                       '--transpose',
                       '--out', reap_tped])
                       
### attach FID/IIDs to mixture proportions

# verify files are same length
target_Qfile_nam = str(args.target_bfile + '.' + str(args.npops) + '.Q')                     
target_fam_nam = str(args.target_bfile + '.fam')

if not (file_len(target_Qfile_nam) == file_len(target_fam_nam)):
    raise ValueError('Length of admixture proportions ouput (%s) does not match fam file (%s). ' + \
                     'Error during output?' % (target_Qfile_nam, target_fam_nam))

# paste together columns, should be in same order (based on ADMIXTURE's ouptut format)
target_Q_file = open(target_Qfile_nam, 'r')
target_fam = open(target_fam_nam, 'r')
reap_mix_props = open(str(args.target_bfile + '.props.tmp.txt'), 'w')


for line in target_Q_file:
    qprops = line.split()
    (fid, iid, pat, mat, sex, phen) = target_fam.readline().split()
    
    reap_mix_props.write(str(fid) + ' ' + str(iid) + ' ' + ' '.join(qprops) + '\n')

reap_mix_props.close()

# verify hit end of fam file (readline should return false)
if target_fam.readline():
    raise IOError('Reached end of admixture proportions (%s) before end of fam file (%s). ' + \
                  'Problems with parsing?' % (target_Qfile_nam, target_fam_nam))

target_Q_file.close()
target_fam.close()


#############
print '\n...Estimating relatedness...'
#############

# note: REAP docs claim '-r 1' needed when using ADMIXTURE files
# but testing gives gives correct estimates with '-r 2', not '-r 1'
reap_call = [str(args.reap_ex),
             '-g', str(reap_tped + '.tped'),
             '-p', str(reap_tped + '.tfam'),
             '-a', str(args.target_bfile + '.props.tmp.txt'),
             '-f', str(args.target_bfile + '.' + str(args.npops) + '.P'),
             '-r', str(2),
             '-k', str(args.npops),
             '-m',
             '-t', str(args.min_rel)]
reap_log = open(str('reap_' + args.out + '.log'), 'w')

print str(' '.join(reap_call))
print 'Logging to ' + reap_log.name + '\n'
subprocess.check_call(reap_call, stdout=reap_log)

reap_log.close()



#############
print '\n...Generating diagnostic plots...'
#############

### IBD0/IBD1 points and density
# plot_reap_ibd.Rscript has args <input_file> <outname> <minimum relatedness>
r_ibd_log = open(str(args.out) + '.plot_ibd.log', 'w')
subprocess.check_call([Rplotibdx,
                       str('REAP_pairs_relatedness.txt'),
                       str(args.out),
                       str(args.min_rel)],
                       stderr=subprocess.STDOUT,
                       stdout=r_ibd_log)

r_ibd_log.close()
print 'IBD plots: %s.IBD.png, %s.IBD_density.png' % (args.out, args.out)

### exemplars on PCA

### admixture on PCA

# Generate diagnostic plots
# - exemplars on PCA (if PCA available)
# - final admixture on PCA (if PCA available)
# - IBD0/IBD1 for REAP

# final cleanup

# finish