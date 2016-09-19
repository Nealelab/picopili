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
import argparse
from glob import glob
from py_helpers import read_conf, unbuffer_stdout, test_exec, find_from_path
from args_pca import *
unbuffer_stdout()

#############
if not (('-h' in sys.argv) or ('--help' in sys.argv)):
    print '\n...Parsing arguments...' 
#############

### init vars that may be set as functions of others
pcadir = ""

### parse arguments
parser = argparse.ArgumentParser(prog='imus_pca.py',
                                 formatter_class=lambda prog:
                                 argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=40),
                                 parents=[parserbase, parserpca])
args = parser.parse_args()


# set remaining defaults
if args.pcadir == None or args.pcadir == "None":
    pcadir = args.out + '_imus_pca'
else:
    pcadir = args.pcadir
    
wd = os.getcwd()

# plot settings
default_nplot = 6
case_color = "orange"
case_pch = "3"
con_color = "blue"
con_pch = "1"
miss_color = "green"
miss_pch = "2"
other_color = "black"
other_pch = "4"


# print settings
print 'Using settings:'
print '--bfile '+args.bfile
print '--out '+args.out
print '--rel-th '+str(args.rel_th)
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

# find required files
if args.rscript_ex == None or args.rscript_ex == "None":
    args.rscript_ex = find_from_path("Rscript", 'Rscript')

Rplotpcax = find_from_path("plot_pca.Rscript", 'PCA plotting script')


# test executables
test_exec(args.primus_ex, 'PRIMUS')
test_exec(plinkx, 'Plink')
test_exec(smartpcax, 'Eigensoft smartpca')
test_exec(args.rscript_ex, 'Rscript')



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

primelog = open(str('primus_' + args.out + '_imus.log'), 'w')
subprocess.check_call([args.primus_ex,
                       "--file", args.bfile,
                       "--genome",
                       "--rel_threshold", str(args.rel_th),
                       "--no_PR",
                       "--plink_ex", plinkx,
                       "--smartpca_ex", smartpcax,
                       "--output_dir", str(args.out+'_primus')],
                       stderr=subprocess.STDOUT,
                       stdout=primelog)
                       
primelog.close()

# verify successful output
primedir = os.getcwd() + '/' + args.out + '_primus'
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
os.symlink(imus_dirfile, imus_file)
os.symlink(wd+'/'+args.bfile+'.bed', args.bfile+'.bed')
os.symlink(wd+'/'+args.bfile+'.bim', args.bfile+'.bim')
os.symlink(wd+'/'+args.bfile+'.fam', args.bfile+'.fam')

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
                       "--silent",
                       "--memory", str(2000),
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

#############
print '\n...Preparing files for PCA...'
#############

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
id_conv = open(str(args.bfile)+'.pca.pedind.ids.txt', 'w')
id_conv.write('FID IID pca_id' + '\n')

# init pedind file
pedind = open(str(args.bfile)+'.pca.pedind', 'w')

# init dict for converting back
id_dict = {}

# track sex to add to pca output
sex_dict = {}

# process fam file by line
bfile_fam = open(str(args.bfile+'.fam'), 'r')
idnum=0
for line in bfile_fam:
    # iterate line counter, used for short id
    idnum += 1
    
    # read
    (longfid, longiid, pat, mat, sex, phen) = line.split()

    # assign and record short identifier
    shortfid = str(idnum)
    shortiid = str(idnum)
    pca_id = shortfid +':'+ shortiid
    id_conv.write(longfid +' '+ longiid +' '+ pca_id + '\n')

    # get FID:IID identifier to compare to imus
    bfile_id = longfid +':'+ longiid
    
    # write pedind with "IMUS" if match, "RELATED" if not
    if any(bfile_id == refid for refid in imus_ids):
        pedind.write(shortfid +' '+ shortiid +' 0 0 U IMUS\n')
    else:
        pedind.write(shortfid +' '+ shortiid +' 0 0 U RELATEDS\n')
    
    # record id pair in dict for converting back
    id_dict[pca_id] = bfile_id
    sex_dict[pca_id] = sex

pedind.close()
id_conv.close()


### process bim file to remove any centi-morgan distances
bim_nocm = open(str(args.bfile+'.nocm.bim'), 'w')

with open(str(args.bfile+'.bim'), 'r') as bim_in:
    for bimline in bim_in:
        (chrom, snp, cm, bp, a1, a2) = bimline.split()
        bim_nocm.write(' '.join([chrom, snp, str(0), bp, a1, a2]) + '\n')
        
bim_nocm.close()


### create par file
par = open(str(args.bfile + '.pca.par'), 'w')

par.write('genotypename:     '+str(args.bfile+'.bed')+'\n')
par.write('snpname:          '+str(args.bfile+'.nocm.bim')+'\n')
par.write('indivname:        '+str(args.bfile+'.pca.pedind')+'\n')
par.write('poplistname:      '+str(args.bfile+'.pca.refpoplist.txt')+'\n')
par.write('evecoutname:      '+str(args.bfile+'.pca.raw.txt')+'\n')
par.write('evaloutname:      '+str(args.bfile+'.pca.eval.txt')+'\n')
par.write('snpweightoutname: '+str(args.bfile+'.pca.snpw.txt')+'\n')
par.write('fastmode:         '+'YES'+'\n')
par.write('altnormstyle:     '+'NO'+'\n')
par.write('numoutevec:       '+str(args.npcs)+'\n')
par.write('numoutlieriter:   '+str(0)+'\n')

par.close()


### create poplist file
poplist = open(str(args.bfile+'.pca.refpoplist.txt'), 'w')
poplist.write("IMUS\n")
poplist.close()


#############
print '\n...Running PCA...'
#############

### run smartpca
pcalog = open(str('smartpca_'+args.out+'_imusproj.log'), 'w')
subprocess.check_call([smartpcax, 
                       "-p", str(args.bfile + '.pca.par')],
                       stderr=subprocess.STDOUT,
                       stdout=pcalog)

pcalog.close()


#############
print '\n...Processing PCA results files...'
#############

# final formatted pca results
pc_out = open(str(args.out+'.pca.txt'), 'w')

# init header
# ouput is FID, IID, PCs 1-npcs; all space sep.
pc_out.write('FID IID ' + ' '.join( [str('PC'+str(i)) for i in xrange(1,args.npcs+1)] ) + ' SEX\n')

# read in smartpca output, strip header, convert back fid:iid, remove pop designation
with open(str(args.bfile+'.pca.raw.txt'), 'r') as evec:
    dumphead = evec.readline()
    for line in evec:
        pc_in = line.split()

        # verify num fields matches fid:iid + PCs + pop label
        if len(pc_in) != (int(args.npcs)+2):
            raise ValueError("Incorrect number of fields in %s (expected %d, found %d)" % (evec.name, int(args.npcs)+2, len(pc_in)))

        # convert back fid:iid
        matched_id = id_dict[pc_in[0]]
        out_id = matched_id.split(':')
        
        # get sex
        matched_sex = str(sex_dict[pc_in[0]])
        
        # verify nothing weird with resulting match
        if len(out_id) != 2:
            raise ValueError("Problem parsing fid:iid for %r, check %s?" % (pc_in[0], str(args.bfile+'.pca.pedind.ids.txt')))
        
        # print
        pc_out.write( out_id[0] +' '+ out_id[1] +' '+ ' '.join(pc_in[1:(int(args.npcs)+1) ]) +' '+ matched_sex + '\n' )

pc_out.close()


####################################
# Plot projected PCs
# a) create temp files of PCA results with plotting instructions
# b) plot with R
####################################

#############
print '\n...Preparing to plot PCs...'
#############

### Create input files for plot_pca.Rscript
# - Plotting info file, with columns: FID, IID, col, pch, layer
# - Legend file, with columns: col, pch, fill, text (either pch/col or fill = NA)

# plot info, from fam
plotinfo = open(str(args.out+'.pca.plotinfo.txt'), 'w')

# header
plotinfo.write('FID IID col pch layer\n')

# track if need legend for missing, "other"
anymiss = False
anyother = False

with open(str(args.bfile+'.fam'), 'r') as fam:
    for line in fam:
        (fid, iid, pat, mat, sex, phen) = line.split()
        if int(phen)==1:
            plotinfo.write(' '.join([fid, iid, con_color, con_pch, str(1)]) +'\n')
        elif int(phen)==2:
            plotinfo.write(' '.join([fid, iid, case_color, case_pch, str(2)]) +'\n')
        elif int(phen)==0 or int(phen)==-9:
            plotinfo.write(' '.join([fid, iid, miss_color, miss_pch, str(0)]) +'\n')
            anymiss = True
        else:
            plotinfo.write(' '.join([fid, iid, other_color, other_pch, str(-1)]) +'\n')
            anyother = True

plotinfo.close()

# legend
legend = open(str(args.out)+'.pca.legend.txt', 'w')

legend.write('col pch fill text\n') # header
legend.write(con_color+' '+con_pch+' NA \"control, '+str(args.bfile+'.fam')+'\"\n')
legend.write(case_color+' '+case_pch+' NA \"case, '+str(args.bfile+'.fam')+'\"\n')
if anymiss:
    legend.write(miss_color+' '+miss_pch+' NA \"missing, '+str(args.bfile+'.fam')+'\"\n')
if anyother:
    legend.write(other_color+' '+other_pch+' NA \"other (quant?), '+str(args.bfile+'.fam')+'\"\n')

legend.close()

# number of PCs to plot
if args.plot_all or (int(args.npcs) <= 6):
    nplot = int(args.npcs)
else:
    nplot = default_nplot


#############
print '\n...Plotting PCA results...'
#############
# Args for plot_pca.Rscript are:
# - PCA results file name
# - Plotting info file name
# - Legend file name
# - number of PCs to plot
# - output name stem

if not os.path.exists("plots"):
    os.makedirs("plots")

rplotlog = open(str('rplot_'+args.out+'_pca.log'), 'w')
subprocess.check_call([Rplotpcax,
                       str(args.out+'.pca.txt'),
                       str(args.out+'.pca.plotinfo.txt'),
                       str(args.out+'.pca.legend.txt'),
                       str(nplot),
                       str(args.out)],
                       stderr=subprocess.STDOUT,
                       stdout=rplotlog)

rplotlog.close()


####################################
# Clean up files
####################################

os.chdir(wd)

if not args.no_cleanup:
    
    #############
    print '\n...Cleaning up output files...'
    #############

    #############
    print 'Zipping ' + pcadir + '/' + args.out + '.pca_files.tar.gz:'
    #############
    os.chdir(pcadir)
    subprocess.check_call(["tar", "-zcvf",
                           args.out + '.pca_files.tar.gz',
                           args.bfile + '.pca.par',
                           args.bfile + '.pca.eval.txt',
                           args.bfile + '.pca.snpw.txt',
                           args.bfile + '.pca.raw.txt',
                           args.bfile + '.pca.refpoplist.txt',
                           args.bfile + '.nocm.bim',
                           args.bfile + '.pca.pedind',
                           args.bfile + '.pca.pedind.ids.txt'])

    # remove successfully zipped files
    subprocess.check_call(["rm",
                           args.bfile + '.pca.par',
                           args.bfile + '.pca.eval.txt',
                           args.bfile + '.pca.snpw.txt',
                           args.bfile + '.pca.raw.txt',
                           args.bfile + '.pca.refpoplist.txt',
                           args.bfile + '.nocm.bim',
                           args.bfile + '.pca.pedind',
                           args.bfile + '.pca.pedind.ids.txt'])
                           
    os.chdir(wd)

    
    #############
    print '\nCompressing PCA results:'
    #############
    os.chdir(pcadir)
    subprocess.check_call(["gzip", "-fv", str(args.out + '.pca.txt')])
    os.chdir(wd)


    #############
    print '\nLinking useful PRIMUS plots to ./' + pcadir + '/plots/'
    #############
    os.chdir(pcadir)
    os.symlink(wd+'/'+args.out+'_primus/'+args.bfile+'_prePRIMUS/'+args.bfile+'_noDups_autosomal_unrelateds_PCV1vPCV2.pdf', './plots/'+args.bfile+'_noDups_autosomal_unrelateds_PCV1vPCV2.pdf')
    os.symlink(wd+'/'+args.out+'_primus/'+args.bfile+'_prePRIMUS/'+args.bfile+'_noDups_autosomal_unrelateds_merged_KDE_contours.pdf', './plots/'+args.bfile+'_noDups_autosomal_unrelateds_merged_KDE_contours.pdf')
    os.symlink(wd+'/'+args.out+'_primus/'+args.bfile+'_prePRIMUS/'+args.bfile+'_noDups_autosomal_unrelateds_merged_PCV1vPCV2.pdf', './plots/'+args.bfile+'_noDups_autosomal_unrelateds_merged_PCV1vPCV2.pdf')
    os.symlink(wd+'/'+args.out+'_primus/'+args.bfile+'_prePRIMUS/'+args.bfile+'_cleaned.genome_IBD0_vs_IBD1.jpeg', './plots/'+args.bfile+'_cleaned.genome_IBD0_vs_IBD1.jpeg')
    os.chdir(wd)


    #############
    print '\nRemoving interim files:'
    #############
    os.chdir(primedir)
    print 'Files in ./' + args.out + '_primus/:'
    subprocess.check_call(["rm", "-v"] + glob(args.bfile+"_cleaned.genome_network*.genome"))
    subprocess.check_call(["rm", "-v"] + glob(args.bfile+'_prePRIMUS/'+args.bfile+"_*.bed~"))
    subprocess.check_call(["rm", "-v"] + glob(args.bfile+'_prePRIMUS/'+args.bfile+"_*.bim~"))
    subprocess.check_call(["rm", "-v"] + glob(args.bfile+'_prePRIMUS/'+args.bfile+"_*.fam~"))
    subprocess.check_call(["rm", "-v"] + glob(args.bfile+'_prePRIMUS/'+args.bfile+"*.log"))
    subprocess.check_call(["rm", "-v"] + glob(args.bfile+'_prePRIMUS/'+args.bfile+"_cleaned.genome"))
    os.chdir(wd)
    
    os.chdir(pcadir)
    print '\nFiles in ./' + pcadir + '/:'
    subprocess.check_call(["rm", "-v",
                           args.out + '.pca.legend.txt',
                           args.out + '.pca.plotinfo.txt'])  
    os.chdir(wd)
   
    #############
    print '\nRemove if exist:'
    #############
    # allowing failure, since files may or may not exists    
    # additional primus files to remove
    os.chdir(primedir)
    print 'Files in ./' + args.out + '_primus/:'
    subprocess.call(["rm", "-r",
                     str(args.bfile+'_prePRIMUS/'+args.bfile+'_noDups_autosomal_IMUS')])
    subprocess.call(["rm", "-v"] + glob(args.bfile+'_prePRIMUS/'+args.bfile+"*.nosex"))
    subprocess.call(["rm", "-v"] + glob(args.bfile+'_prePRIMUS/'+args.bfile+"_*.het"))
    subprocess.call(["rm", "-v"] + glob(args.bfile+'_prePRIMUS/'+args.bfile+"_*.*miss"))
    subprocess.call(["rm", "-v"] + glob(args.bfile+'_prePRIMUS/'+args.bfile+".fam_temp"))
    os.chdir(wd)    
    
    # additional pca files to remove
    os.chdir(pcadir)
    print '\nFiles in ./' + pcadir + '/'
    subprocess.call(["rm", "-v",
                     str(bfile_imus +'.nosex'),
                     str(bfile_imus +'.hh')])
    os.chdir(wd)




print '\n############'
print '\n'
print 'SUCCESS!\n'
exit(0)
