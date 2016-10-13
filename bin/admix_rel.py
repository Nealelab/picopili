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
import argparse
from string import ascii_uppercase
from glob import glob
from numpy import digitize
from py_helpers import unbuffer_stdout, file_len, test_exec, find_exec, link, gz_confirm
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
                         'with unrelated individuals to estimate admixture.' + \
                         'Must specify either this or --admix-p.',
                    required=False)
arg_base.add_argument('--admix-p',
                    type=str,
                    metavar='FILE',
                    help='Admixture results .P file from sample of ' + \
                         'unrelated individuals. Can alternatively specify ' + \
                         '--unrel-bfile to run this initial admixture.',
                    required=False)
arg_base.add_argument('--admix-q',
                    type=str,
                    metavar='FILE',
                    help='Admixture results .Q file from sample of ' + \
                         'unrelated individuals. Required only if using ' + \
                         '--admix-p and --use-exemplars.',
                    required=False)
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
arg_admix.add_argument('--use-exemplars',
                    action='store_true',
                    help='Determine admixture in target sample based on ' + \
                          'supervised fit with a selection of population exemplars ' + \
                          'rather than a project of admixture solution in unrelateds. ' + \
                          '(Required for ADMIXTURE version < 1.3). Requires --unrel-bfile, ' + \
                          'and if using --admix-p also requires specifying --admix-q.')                    
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
                    default=None)
arg_exloc.add_argument('--reap-ex',
                    type=str,
                    metavar='PATH',
                    help='path to REAP executable',
                    required=False,
                    default=None)

args = parser.parse_args()

# set dependent defaults
if args.outdir == None or args.outdir == "None":
    args.outdir = str(args.out)+'_admix_rel'

if not (args.plot_admix_pca == None or args.plot_admix_pca == "None"):
    plot_pca = True
else:
    plot_pca = False


# plot settings
exemplar_pch = "20"
ref_pch = "1"
other_pch = "3"

exemplar_color = "orange"
ref_color = "black"
other_color = "gray50"

# Adapted from RColorBrewer YlGnBu to trim lighter colors, increase resolution
# t(
#   apply(
#       colorRamp(brewer.pal(9,"YlGnBu")[3:9],space="rgb")(1:10/10),
#       1, function(a) as.character(as.hexmode(round(a))
#   )
# )
col_gradient = ["#9CD8B8","#73C8BD","#4DBBC2","#33A7C2","#1D91C0","#2072B2","#2356A4","#243C98","#192B7C","#081D58"]


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
print '\n...Checking dependencies...'
# find, check exists, executable
#############

plinkx = find_exec('plink',key='p2loc')

if args.admixture_ex is None or args.admixture_ex == "None":
    args.admixture_ex = find_exec('admixture', key='admloc')

test_exec(args.admixture_ex, 'ADMIXTURE')

if args.rscript_ex is None or args.rscript_ex == "None":
    args.rscript_ex = find_exec('Rscript', key='rscloc')

if args.reap_ex is None or args.reap_ex == "None":
    args.reap_ex = find_exec('REAP', key='reaploc')

rp_bin = os.path.dirname(os.path.realpath(__file__))
Rplotibdx = rp_bin+'/plot_reap_ibd.Rscript'

if plot_pca:
    Rplotpcax = rp_bin+'/plot_pca.Rscript'

# check if running admixture for unrelateds
run_admix = True
if args.admix_p is not None and args.admix_p != "":
    run_admix = False

else:
    assert os.path.isfile(args.admix_p), "Admixture .P file %s does not exist." % str(args.admix_p)
    
    if args.use_exemplars:
        assert os.path.isfile(args.admix_q), "Admixture .Q file %s does not exist." % str(args.admix_q)


# check if have unrel-bfile if needed:
if args.unrel_bfile is None or args.unrel_bfile == "":
    
    if run_admix:
        raise parser.error('Must specify either --unrel-bfile or --admix-p.')

    if args.use_exemplars:
        raise parser.error('Must specify --unrel-bfile to define exemplars for --use-exemplars.')

else:
    assert '/' not in args.unrel_bfile, "--unrel-bfile must specify only a file stem, not a path"
    assert os.path.isfile(str(args.unrel_bfile)+'.bed'), "bed file for unrelated individuals %s does not exist." % str(args.unrel_bfile)+'.bed'
    assert os.path.isfile(str(args.unrel_bfile)+'.bim'), "bim file for unrelated individuals %s does not exist." % str(args.unrel_bfile)+'.bim'
    assert os.path.isfile(str(args.unrel_bfile)+'.fam'), "fam file for unrelated individuals %s does not exist." % str(args.unrel_bfile)+'.fam'


# verify executables
test_exec(plinkx, 'Plink')
test_exec(args.rscript_ex, 'Rscript')
test_exec(args.reap_ex, 'REAP')

# pca file
if plot_pca:
    assert os.path.isfile(args.plot_admix_pca), "PCA file does not exist (%r)" % args.plot_admix_pca
    assert '/' not in args.target_bfile, "--plot-admix-pca must specify only a file, not a path"

# verify bfiles are files, not paths
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

# link plink files (with verification)
if run_admix or args.use_exemplars:
    link(str(wd+'/'+args.unrel_bfile+'.bed'), str(args.unrel_bfile+'.bed'), 'bed file for unrelated individuals')
    link(str(wd+'/'+args.unrel_bfile+'.bim'), str(args.unrel_bfile+'.bim'), 'bim file for unrelated individuals')
    link(str(wd+'/'+args.unrel_bfile+'.fam'), str(args.unrel_bfile+'.fam'), 'fam file for unrelated individuals')

link(str(wd+'/'+args.target_bfile+'.bed'), str(args.target_bfile+'.bed'), 'bed file for target individuals')
link(str(wd+'/'+args.target_bfile+'.bim'), str(args.target_bfile+'.bim'), 'bim file for target individuals')
link(str(wd+'/'+args.target_bfile+'.fam'), str(args.target_bfile+'.fam'), 'fam file for target individuals')

# link pca file, if provided
if not (args.plot_admix_pca==None or args.plot_admix_pca=="None"):

    link(str(wd+'/'+args.plot_admix_pca), str(args.plot_admix_pca), 'PCA file')



#############
print '\n...Running Admixture on unrelated dataset...'
#############

if run_admix:
    admix_call = [args.admixture_ex,
                  str(args.unrel_bfile+'.bed'),
                  str(args.npops),
                  '-j'+str(args.multithread_cores)]
    admix_unrel_log = open(str('admix_'+args.out+'_unrel.log'), 'w')
    
    print str(' '.join(admix_call))
    print 'Logging to ' + admix_unrel_log.name + '\n'
    subprocess.check_call(admix_call, stdout=admix_unrel_log)
    
    admix_unrel_log.close()


if args.use_exemplars:

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
    
    if run_admix:
        admix_pops_file = str(args.unrel_bfile+'.'+str(args.npops)+'.Q')
    else:
        admix_pops_file = args.admix_q
    
    
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
    


# no exemplars, using projection instead
else:
    
    #############
    print '\n...Projecting admixture analysis to target data...'
    #############
    
    ref_p_name = str(args.target_bfile)+'.'+str(args.npops)+'.P.in'
    if run_admix:
        link(str(args.unrel_bfile)+'.'+str(args.npops)+'.P', ref_p_name, 'admixture allele freqs')
    else:
        ref_p_in = str(args.admix_p)
        link(wd+'/'+ref_p_in, ref_p_name,'input admixture allele freqs')

    
    admix_project_call = [args.admixture_ex,
                          '-P', str(args.target_bfile+'.bed'),
                        str(args.npops),
                        '-j'+str(args.multithread_cores)]
    admix_target_log = open(str('admix_'+args.out+'_target.log'), 'w')
    
    print str(' '.join(admix_project_call))
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
                       '--allow-no-sex',
                       '--recode12',
                       '--output-missing-genotype', '0',
                       '--transpose',
                       '--out', reap_tped])
                       
### attach FID/IIDs to mixture proportions

# verify files are same length
target_Qfile_nam = str(args.target_bfile + '.' + str(args.npops) + '.Q')                     
target_fam_nam = str(args.target_bfile + '.fam')

if not (file_len(target_Qfile_nam) == file_len(target_fam_nam)):
    raise ValueError('Length of admixture proportions output (%s) does not match fam file (%s). ' + \
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

# create plot directory
os.makedirs("plots")

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


### Optional plots on PCA  
if plot_pca:
    ###
    # - create plotinfo files for each ancestry population
    #    - parse admixture results
    #    - get population proportion
    #    - check if exemplar, in (imus) reference set
    # - create legend files
    # - loop 1-npops ancestries for plotting
    #    - plot admixture proportion on pcs 1-3
    #       - split to 8 bins for color gradient
    #    - plot exemplar status on pcs 1-3
    ###

    # setup file streams for plotinfo files
    pop_info_files = []
    if args.use_exemplars:
        exemp_info_files = []
    for i in xrange(args.npops):
        pop_info_files.append( open(str(args.target_bfile) + '.' + popnames[i] + '.admixture.plotinfo.txt', 'w') )
        pop_info_files[i].write('FID IID col pch layer\n')
        
        if args.use_exemplars:
            exemp_info_files.append( open(str(args.target_bfile) + '.' + popnames[i] + '.exemplar.plotinfo.txt', 'w') )
            exemp_info_files[i].write('FID IID col pch layer\n')
        
    # parse admixture proportions
    reap_mix_props = open(str(args.target_bfile + '.props.tmp.txt'), 'r')
    for line in reap_mix_props:
        splitline = line.split()
        
        # pull off fid/iid in first 2 columns
        fid = str(splitline.pop(0))
        iid = str(splitline.pop(0))
        joinid = fid + ':' + iid
        
        # loop populations
        for i in xrange(args.npops):
            # bin the population proportion
            # Note: .tolist()[0] need to get back from np.array to scalar
            prop_bins = [float(x) / 10.0 for x in range(0,11)]
            in_bin = digitize([splitline[i]], prop_bins).tolist()[0]
            bin_col = col_gradient[in_bin-1]
            
            # admix proportion info file: FID, IID, col, pch, layer
            pop_info_files[i].write(' '.join([fid, iid, bin_col, str(1), str(in_bin)])+'\n')
            
            # exemplar info file: FID, IID, col, pch, layer
            if args.use_exemplars:
                if joinid in pop_dict:
                    if pop_dict[joinid] == popnames[i]:
                        exemp_info_files[i].write(' '.join([fid, iid, '\"'+str(exemplar_color)+'\"', str(exemplar_pch), str(3)]) + '\n')
                    else:
                        exemp_info_files[i].write(' '.join([fid, iid, '\"'+str(ref_color)+'\"', str(ref_pch), str(2)]) + '\n')
                else:
                    exemp_info_files[i].write(' '.join([fid, iid, '\"'+str(other_color)+'\"', str(other_pch), str(1)]) + '\n')

    # close plotinfo files
    for i in xrange(args.npops):
        pop_info_files[i].close()
        if args.use_exemplars:
            exemp_info_files[i].close()
    
    # create legend files: col, pch, fill, text (either col/pch or fill should be NA)
    if args.use_exemplars:
        exem_legend = open(str(args.target_bfile) + '.exemplar.legend.txt', 'w')
        exem_legend.write('col pch fill text\n')
        exem_legend.write(str(exemplar_color) + ' ' + str(exemplar_pch) + ' NA ' + '\"Population exemplar\"\n')
        exem_legend.write(str(ref_color) + ' ' + str(ref_pch) + ' NA ' + '\"Reference set\"\n')
        exem_legend.write(str(other_color) + ' ' + str(other_pch) + ' NA ' + '\"Non-reference set\"\n')        
        exem_legend.close()
    
    prop_legend = open(str(args.target_bfile) + '.admixture.legend.txt', 'w')
    prop_legend.write('col pch fill text\n')
    for i in xrange(len(col_gradient)):    
        prop_legend.write('NA NA ' + str(col_gradient[i]) + ' ' + str(prop_bins[i])+'-'+str(prop_bins[i+1])+ '\n')       
    prop_legend.close()

    ### generate plots
    for i in xrange(args.npops):
        if args.use_exemplars:
            r_pca_ex_log = open(str(args.out) + '.' + popnames[i] + '.plot_exemplars.log', 'w')
            subprocess.check_call([Rplotpcax,
                                   str(args.plot_admix_pca),
                                   str(args.target_bfile) + '.' + popnames[i] + '.exemplar.plotinfo.txt',
                                   str(args.target_bfile) + '.exemplar.legend.txt',
                                   str(3),
                                   str(args.out) + '.' + popnames[i] + '.exemplars'],
                                   stderr=subprocess.STDOUT,
                                   stdout=r_pca_ex_log)    
            r_pca_ex_log.close()

        r_pca_admix_log = open(str(args.out) + '.' + popnames[i] + '.plot_admixture.log', 'w')
        subprocess.check_call([Rplotpcax,
                               str(args.plot_admix_pca),
                               str(args.target_bfile) + '.' + popnames[i] + '.admixture.plotinfo.txt',
                               str(args.target_bfile) + '.admixture.legend.txt',
                               str(3),
                               str(args.out) + '.' + popnames[i] + '.admixture'],
                               stderr=subprocess.STDOUT,
                               stdout=r_pca_admix_log)    
        r_pca_admix_log.close()
        print 'PCA plots for %s: %s, %s (completed %d/%d populations)' % (popnames[i], str(args.out)+'.'+popnames[i]+'.exemplars.pca.pairs.png', str(args.out)+'.'+popnames[i]+'.exemplars.pca.pc##_pc##.png', i+1, args.npops)



# final cleanup
if not args.no_cleanup:
    
    #############
    print '\n...Cleaning up output files...'
    #############
    
    if plot_pca:
        ###
        print '\nZipping PCA plot files to ' + args.out + '.pca_plot_files.tar.gz:'
        ###        
        subprocess.check_call(["tar", "-zcvf", 
                               str(args.out+'.plot_pca_files.tar.gz')] + \
                               glob(args.target_bfile+".*.admixture.plotinfo.txt") + \
                               [str(args.target_bfile)+".admixture.legend.txt"] + \
                               glob(args.out+".*.plot_admixture.log"))
                               
        subprocess.check_call(["tar", "-zcvf",
                               str(args.out+'.plot_exemplar_files.tar.gz')] + \
                               glob(args.target_bfile+".*.exemplar.plotinfo.txt") + \
                               [str(args.target_bfile)+".exemplar.legend.txt"] + \
                               glob(args.out+".*.plot_exemplars.log")  )
        
        # remove files after successfully compressed
        subprocess.check_call(['rm'] + glob(args.target_bfile+".*.admixture.plotinfo.txt"))
        subprocess.check_call(['rm'] + glob(args.target_bfile+".admixture.legend.txt"))
        subprocess.check_call(['rm'] + glob(args.out+".*.plot_admixture.log"))
        
        if args.use_exemplars:
            subprocess.check_call(['rm'] + glob(args.target_bfile+".*.exemplar.plotinfo.txt"))
            subprocess.check_call(['rm'] + glob(args.target_bfile+".exemplar.legend.txt"))
            subprocess.check_call(['rm'] + glob(args.out+".*.plot_exemplars.log"))


    ###
    print '\nZipping Admixture output files:'
    ###
    
    gz_confirm(str(args.target_bfile)+'.'+str(args.npops)+'.P', 
               str(args.target_bfile)+'.'+str(args.npops)+'.P.gz', force=False)
    gz_confirm(str(args.target_bfile)+'.'+str(args.npops)+'.Q', 
               str(args.target_bfile)+'.'+str(args.npops)+'.Q.gz', force=False)
    
    if run_admix:
        gz_confirm(str(args.unrel_bfile)+'.'+str(args.npops)+'.P', 
                   str(args.unrel_bfile)+'.'+str(args.npops)+'.P.gz', force=False)
        gz_confirm(str(args.unrel_bfile)+'.'+str(args.npops)+'.Q', 
                   str(args.unrel_bfile)+'.'+str(args.npops)+'.Q.gz', force=False)

    
    ###
    print '\nZipping REAP output files:'
    ###
    gz_confirm('REAP_Inbreed.txt', 
               str(args.out)+'.REAP_Inbreed.txt.gz', force=False)
    gz_confirm('REAP_Individual_Index.txt', 
               str(args.out)+'.REAP_Individual_Index.txt.gz', force=False)
    gz_confirm('REAP_pairs_relatedness.txt', 
               str(args.out)+'.REAP_pairs_relatedness.txt.gz', force=False)
                     
    
    ###
    print '\nRemoving temporary files:'
    ###
    subprocess.check_call(['rm', '-v',
                           str(args.target_bfile)+'.tmp_recode.tped',
                           str(args.target_bfile)+'.tmp_recode.tfam'])
    
    if plot_pca:
        subprocess.check_call(['rm', '-v', 
                               str(args.target_bfile)+'.props.tmp.txt'])

    ###
    print '\nRemove if exist:'
    ###
    subprocess.call(['rm', '-v', str(args.target_bfile)+'.tmp_recode.nosex'])
    subprocess.call(['rm', '-v', str(args.target_bfile)+'.tmp_recode.hh'])
    

# finish
print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)
