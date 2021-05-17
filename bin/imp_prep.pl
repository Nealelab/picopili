#!/usr/bin/env perl
use strict;

####
# imp_prel.pl
# closely adapted from ricopili's impute_dirsub_57
# Original by Stephan Ripke
# Adapted by Raymond Walters
#
# 
# Changes in this adaptation:
# - stop before prephasing/imputation (buigue/chepos/chefli only)
# - remove old version history
# - remove code for fam file prep/chunking/refind and beyond
# - remove commented out blocks, unused subroutines
# - remove --phase (now using --reffiles)
# - remove unused arguments
# - remove non-"readref" alignment
# - remove checks for unused scripts
#
# TODO: refine output format
# 
#### 


#############################
# load utility functions
#############################

use FindBin;
use lib "$FindBin::Bin";
use rp_perl::Utils qw(trans);

my $version = "1.0.24";
my $progname = $0;
$progname =~ s!^.*/!!;
my $command_line = "$progname @ARGV";

use Cwd;
use File::Path;
my $rootdir = &Cwd::cwd();
my $sjainfotxt = "$rootdir\t$command_line";

my $jnum = 7; ### number of imputation job per node

my $multithread1 = 4;
my $multithread2 = 8;

use Sys::Hostname;
my $host = hostname;

#my $broad = 1 if ($ENV{"DOMAINNAME"} =~ /broadinstitute/);


#############################
# read config file
#############################

my $ploc = &trans("p2loc");
my $qloc = &trans("cluster");
my $email = &trans("email");

my $email_on = 0;
if ($email =~ m/\@/){
	$email_on = 1;
}

###############################################

my $iname = "" ;

my $suminfo = "infosum_pos";
my $suminfo_s = "NAN";

my $job_bn_th = 1000;

my $sjamem_incr = 0;
my $sec_freq = .2;  ## secure freq distance around 50%
my $fth_th = 0.15;
my $popname = "eur";



##### help message
my $usage = "
Usage : $progname [options] --bim FILE --reffiles STRING --outname OUTNAME

version: $version


 --help            print this text and exits
 
 --bim FILE        bim file of dataset to prepare (in current directory)

 --reffiles STRING   full path of reference directory (mandatory)
                            can use \"###\" as chr placeholder 

 --outname STRING  identifier for imputation run (mandatory)


##### alignment to reference:

  --popname STRING    important for freq-checks, either 
                           eur (default), 
                           eas (east asian),
						   sas (south asian),
                           amr (america), 
                           afr (africa), 
                           all (combined populations)

  --sfh FLOAT         secure frequency around 50% (default: $sec_freq)
                                for checkflip (compare to reference),
                                only applied to AT/CG SNPs
  --fth FLOAT         frequency-diff to exclude SNPs, default $fth_th
                                for checkflip (compare to reference)


### technical options

  --sjamem_incr INT   increase all memore requests by INT Mb in steps of 1000 (1000 is 1Gb)

  --force1            do not exit if same fail, but do this only once

  --sleep INT         sleep for INT seconds, try this if you think there is a race condition
                       (hints: stops a different steps, --serial works)
                       try using 5 seconds first.

  --serial            no sending jobs to queue all in one run
                          -> usually only used for testing  


### remarks 

  --bim, --reffiles, and --outname are mandatory!!

";



use Getopt::Long;
GetOptions( 

    "help"=> \my $help,
	"bim=s"=> \my $bim,
	"reffiles=s"=> \my $reffile_struct,
	"outname=s"=> \my $outname,

	"popname=s"=> \$popname,
    "sfh=f"=> \$sec_freq,
    "fth=f"=> \$fth_th,
	
    "sjamem_incr=i"=> \$sjamem_incr,
    "force1"=> \my $force1,
    "sleep=i"=> \my $sleep_sw,
    "serial"=> \my $serial,
    
    );

my $popname_uc = uc($popname);



if ($sleep_sw) {
  print "sleeping for $sleep_sw seconds (only use if suspect of race condition)\n";
  sleep ($sleep_sw);
}




############################################################
## testing binaries
##############################################################
my @test_scripts;

my $readref_script = "readref_pico.pl";    ### my.pipeline_tar
my $readrefsum_script = "readrefsum_pico.pl";  ### my.pipeline_tar
my $buigue_script = "buigue_pico.pl";              ### my.pipeline_tar
my $lift_script = "lift_to_hg19.pl";
my $checkpos_script = "checkpos_pico.pl";         ### my.pipeline_tar
my $checkflip_script = "checkflip_pico.pl";       ### my.pipeline_tar
my $mutt_script = "mutt";                  ### my.pipeline_tar
my $blue_script = "blueprint.py";             ### my.pipeline_tar

push @test_scripts, $readref_script;
push @test_scripts, $readrefsum_script;
push @test_scripts, $buigue_script;
push @test_scripts, $lift_script;
push @test_scripts, $checkpos_script;
push @test_scripts, $checkflip_script;
push @test_scripts, $blue_script;

#push @test_scripts, $mutt_script ;



print ".......testing necessary binaries....\n";
my @miss_scripts;


#my $err_scr = 0;
foreach my $scr_name (@test_scripts) {
    my $scr_path = '';
    
    for my $path ( split /:/, $ENV{PATH} ) {
	if ( -f "$path/$scr_name" && -x _ ) {
	    print "$scr_name\tfound in $path\n";
	    $scr_path = "$path/$scr_name";
	    last;
	}
    }
    if ( $scr_path eq  '') {
	push @miss_scripts, "$scr_name\n";
	print "!!Error!! : No $scr_name command available\n" ;
    }
 
}



if (@miss_scripts > 0) {

#  if (-e "get_scripts_on_broad.txt") {
#    print "please remove this file and restart: get_scripts_on_broad.txt\n";
#  }
  die $! unless open FILE1, "> missing_picopili_scripts.txt";
  foreach (@miss_scripts) {
    print FILE1 "$_";
  }
  close FILE1;

  die "Missing required scripts. See missing_picopili_scripts.txt\n";

#  print "exiting now -> have a look at get_scripts_on_broad.txt\n";
#  exit;

}




if($email_on){

	print ".......testing email program....\n";

	my $err_scr = 0;
	{
	    my $scr_path = '';
    
	    for my $path ( split /:/, $ENV{PATH} ) {
		if ( -f "$path/$mutt_script" && -x _ ) {
		    print "$mutt_script\tfound in $path\n";
		    $scr_path = "$path/$mutt_script";
		    last;
		}
	    }
	    unless ( $scr_path ) {

		print "!!Warning!! : No $mutt_script command available, trying mail\n" ;

		$mutt_script = "mail";
		for my $path ( split /:/, $ENV{PATH} ) {
		    if ( -f "$path/$mutt_script" && -x _ ) {
			print "$mutt_script\tfound in $path\n";
			$scr_path = "$path/$mutt_script";
			last;
		    }
		}
		unless ( $scr_path ) {
		    $err_scr = 1;
		    print "!!Error!! : No $mutt_script command available\n" ;
		}
	    }
 
	}
	die if $err_scr == 1;

}

print "....all necessary binaries found....\n";
print "------------------------------------\n";
#push @scripts,"id_tager_3";



#####################################
# "testing environment variable rp_perlpackages
####################################

# print "testing environment variable rp_perlpackages....\n";
# unless (exists $ENV{rp_perlpackages}) {
#     print "Error: no environment variable for perl-packages, please re-install ricopili and make sure to follow all instructions\n";
#     print "------------------------------------\n";
#     exit;
# }
print "....all set....\n";
print "------------------------------------\n";




die $usage if $help;

die $usage unless $bim;
die $usage unless $reffile_struct;
die $usage unless $outname;


#########
## check files exist for readref
#########
foreach my $chrloc(1..22) {
    my $reffi = $reffile_struct;
	$reffi =~ s/###/$chrloc/g;
    unless (-e $reffi) {
		die "Error: $reffi not found\n";
    }
}

unless (-e $bim) {
	die "Error: $bim not found\n";
}

my $impute_dir = "pi_sub";


#my $postimp_dir = "$impute_dir/postimp_data";






###################################################
###  system call with test if successfull
###################################################
sub mysystem(){
    my ($systemstr)="@_";
    system($systemstr);
    my $status = ($? >> 8);
    die "$systemstr\n->system call failed: $status" if ($status != 0);
}


##########################################
# subroutine to split a plink-output-line
##########################################

sub split_line {
    my ($line)=@_;
    chomp($line);
    $line =~ s/^[\s]+//g;
    my @cols=  split /\s+/, $line;
}


#####################################
# print array to file with newline
####################################

sub a2filenew {
    my ($file, @lines)=@_;
    die $! unless open FILE, "> $file";
    foreach (@lines){
	print FILE "$_\n";
    }
    close FILE;
}


#####################################
# append array to file with newline
####################################

sub a2filenew_app {
    my ($file, @lines)=@_;
    die "$!: $file" unless open FILE, ">> $file";
    foreach (@lines){
	print FILE "$_\n";
    }
    close FILE;
}


#####################################
# send jobs to cluster and also send navi again
#####################################

my $sjadir = "";
my $sjaweek = 0;
my $sjaname = "";
my $sjarow = "";
my @sjaarray;
my $sjamem = 0;
my $sjamaxpar = 0;

my $sjatime = -1;
my $sjamaxjobs = 30000;


if ($qloc eq "qsub") {
    $sjamaxjobs = 8000;
}


my $sjainfofile = "$rootdir/imp_prep_job_info.log";
unless (-e $sjainfofile) {
	&mysystem ("touch $sjainfofile");
}
my $sjamulti = 0;


sub send_jobarray {

    die "send_jobarray with undefined variables, dir" if ($sjadir eq "");
    die "send_jobarray with undefined variables, name" if ($sjaname eq "");
    die "send_jobarray with undefined variables, array" if (@sjaarray == 0);
    die "send_jobarray with undefined variables, mem" if ($sjamem == 0);
    die "send_jobarray with undefined variables, time" if ($sjatime < 0);
    die "send_jobarray with undefined variables, info" if ($sjainfotxt eq "");


    my $now = localtime time;
    $now =~ s/ /_/g;


#    if ($sjaname eq "finished") {
#
#	my $fini_message ;
#	$fini_message .= "\n\n##################################################################\n";
#	$fini_message .= "##### CONGRATULATIONS: \n";
#	$fini_message .= "##### rp_pipeline finished successfully:\n";
#	$fini_message .= "##### $sjainfotxt\n";
#	$fini_message .= "##### now start with PCA (see README in subdir pcaer_sub/)\n";
#	$fini_message .= "##### or directly with postimputation analysis\n";
#	$fini_message .= "##### have a look at the wiki page\n"; 
#	$fini_message .= "##### https://sites.google.com/a/broadinstitute.org/ricopili/\n";
#	$fini_message .= "##################################################################\n";
#	print "$fini_message\n";
#
#	
#	die $! unless open SUC, "> success_file";
#	print SUC $fini_message."\n";
#	close SUC;
#
#	if($email_on){
#		&mysystem ('cat success_file | '.$mutt_script.' -s RP_pipeline_finished '.$email) ;
#	}
#
#	my $sjarow      = $sjainfotxt."\t$sjaname\t$now";
#	&a2filenew_app("$sjainfofile",$sjarow);
#
#
#	exit;
#
#    }


    chdir ($sjadir);
    my $jobfile = "$sjaname.job_list";
    while (-e $jobfile) {
	$jobfile .= ".s";
	if (@sjaarray < 200) {
	    $sjatime = 4;
	}
    }


    
    &a2filenew ($jobfile, @sjaarray);

#    print "$jobfile\n";
#    exit;

    my $nsja = @sjaarray;


    
    my $nsja_loc = $nsja;
    if ($nsja_loc > $sjamaxjobs) {
	$nsja_loc = $sjamaxjobs;
    }
    
    

    my $multi_txt = "";
    if ($sjamulti > 0) {
	$multi_txt = "--multi $nsja_loc,$sjamulti";
    }

    ### with array
    $sjamem = $sjamem + $sjamem_incr;



    my $sja_week_str = "";
    if ($sjaweek > 0) {
	$sja_week_str = "--week 1";
    }



    
    if ($serial) {
	print "starting step $sjaname with ".@sjaarray." jobs\n";
	print "please be patient.\n";
	my $jc = 1;
	foreach (@sjaarray) {
	    print "running job $jc...\n";
	    &mysystem($_);
	    $jc++;
	    
	}
    }
    else { 
	my $sys_loc = "$blue_script $sja_week_str --maxpar $sjamaxpar --noerr --njob $nsja_loc --array $jobfile --wa $sjatime --mem $sjamem --j --na $sjaname.$outname $multi_txt";
#    print "$sys_loc\n";
#    exit;

	
	&mysystem ($sys_loc);
    }
#    exit;


    my $old_cmd = `tail -1 $sjainfofile | head -1`;

    my $nsja_txt = sprintf "%06d",$nsja;

    my $sjacontent = "$sjaname.".$nsja_txt;

    my $sjarow_part = $sjainfotxt."\t$sjacontent";
    my $sjarow      = $sjainfotxt."\t$sjacontent\t$now";
	
    &a2filenew_app("$sjainfofile",$sjarow);

    if ($old_cmd =~ /$sjarow_part/){
	unless ($force1 ){
	    my $err_message ;
	    $err_message .= "##################################################################\n";
	    $err_message .= "##### Error: \n";
	    $err_message .= "##### step $sjaname has been done repeatedly without any progress\n";
	    $err_message .= "##### imputation pipeline stopped: $command_line\n";
	    $err_message .= "##### $sjainfotxt\n";
#	    $err_message .= "##### if reason does not appear obvious\n";
#	    $err_message .= "##### have a look at the wiki page\n"; 
#	    $err_message .= "##### https://sites.google.com/a/broadinstitute.org/ricopili/\n";
#	    $err_message .= "##### or contact the developers\n";
	    $err_message .= "##################################################################\n";
	    print "$err_message\n";

	    die $! unless open ERR, "> error_file";
	    print ERR $err_message."\n";
	    close ERR;

		if($email_on){
			&mysystem ('cat error_file | '.$mutt_script.' -s Picopili_pipeline_error '.$email) ;
		}

	    unless ($serial) {
		exit;
	    }

	}

    }


    $command_line =~ s/--force1//;


#    my $wt_file = "$sjadir/blueprint_joblist_file-$sjaname.$outname";
    chdir "$rootdir" or die "something strange";
#    if ($qloc eq "bsub") {
#	$wt_file =~ s/.*blueprint_joblist_file-//;
#    }
#
#    if ($qloc eq "slurm") {
#	$wt_file = "$sjadir/$jobfile.script.id";
#    }
#
#    if ($qloc eq "qsub") {
#	$wt_file = "$sjadir/j.$sjaname.$outname.id";
#    }
#    if ($qloc eq "qsub_c") {
#	$wt_file = "$sjadir/j.$sjaname.$outname.id";
#    }
#    if ($qloc eq "qsub_b") {
#	$wt_file = "$sjadir/j.$sjaname.$outname.id";
#    }
    
    my $wt_name = "$sjaname.$outname";

    if ($serial) {
	my $sys_re = "$command_line";
	&mysystem ($sys_re);
	exit;
    }
    else {
	my $sys_re = "$blue_script --njob $job_bn_th -b \"$command_line\" --wa 2 --di -j --wait-name $wt_name --na _if_$outname";
	&mysystem ($sys_re);
    }



    print "------------------------------------------------------------\n";
    print "$nsja jobs successfully submitted\n";
    print "please see tail of $sjainfofile for regular updates\n";
    print "also check bjobs -w for running jobs\n";
    print "possibly differnt command on different computer cluster: e.g. qstat -u USER\n";
    print "you will be informed via email if errors or successes occur\n";
    print "------------------------------------------------------------\n";

    exit;


}




##############################################
##############################################
#############  BEGIN
##############################################
##############################################


unless (-e $impute_dir){
    print "impute_dir is not existing, create one for you\n";
    my @created = mkpath(   ## $created ?
			    $impute_dir,
			    {verbose => 0, mode => 0750},
	);
}


#####################################
## file management
###################################

my @files = ();
opendir(DIR, ".") || die "can't opendir .: $!";
@files = readdir(DIR);
closedir DIR;

my @pi_files = ();

unless (-e "$rootdir/puting_done") {
    opendir(DIR, "$impute_dir") || die "can't opendir .: $!";
    @pi_files = readdir(DIR);
    closedir DIR;
}


### read bim-files
my @bim_files = ($bim);
#print "@bim_files\n";

foreach (@bim_files) {
    if ($_ =~ /.hg19.ch.fl.bim$/){
	print "wrong filename, will rename:\n";
	my $obfile = $_;
	$obfile =~ s/.bim$//;
	my $nbfile = $obfile;
	$nbfile =~ s/.hg19.ch.fl/.bf/;
	print "mv $obfile.bed/bim/fam $nbfile.bed/bim/fam\n";
	&mysystem ("mv $obfile.fam $nbfile.fam");
	&mysystem ("mv $obfile.bed $nbfile.bed");
	&mysystem ("mv $obfile.bim $nbfile.bim");
	print "to redo\n";
	print "mv $nbfile.bim $obfile.bim\n";
	print "mv $nbfile.fam $obfile.fam\n";
	print "mv $nbfile.bed $obfile.bed\n";

	exit;
    }
}

#print "sleep\n";
#sleep(10);

my @bimfli_files = grep {/.ch.fl.bim$/} @pi_files;
my @bimpos_files = grep {/.ch.bim$/} @pi_files;
my @bimref_files = grep {/.bim.ref/} @pi_files;
my @bimhg19_files = grep {/.hg19.bim$/} @pi_files;
if (-e "$rootdir/puting_done") {
#if (@bimfli_files == 0) {
    foreach (@bim_files) {
	my $bitemp = $_;
	$bitemp =~ s/.bim$//;
	$bitemp .= ".hg19.ch.fl.bim";
	push @bimfli_files,$bitemp;
    }
}

#print @bimfli_files."\n";
print "@bimhg19_files\n";
#print "debug\n";
#sleep(10);
#exit;


### read flipped bim-files
my %bimfli_array = ();
foreach (@bimfli_files) {
    $bimfli_array{$_} = 1;
}

### read flipped bim-files
my %bimpos_array = ();
foreach (@bimpos_files) {
    $bimpos_array{$_} = 1;
}

### read flipped bimref-files
my %bimref_array = ();
foreach (@bimref_files) {
    $bimref_array{$_} = 1;
}

### read flipped bim-files
my %bimhg19_array = ();
foreach (@bimhg19_files) {
    $bimhg19_array{$_} = 1;
}



## name for log-files
$iname = $bimfli_files[0];
$iname = $bim_files[0] if ($iname eq "");
$iname =~ s/.bim$//;
$iname =~ s/qc2report_//;


#####################################
# prepare pi_subdir
#####################################

chdir ($impute_dir);

unless (-e "$rootdir/puting_done") {
    foreach (@bim_files) {
	my $bfile = $_;
	$bfile =~ s/.bim$//;
	&mysystem("ln -s $rootdir/$bfile.bim .") unless (-e "$bfile.bim");
	&mysystem("ln -s $rootdir/$bfile.bed .") unless (-e "$bfile.bed");
	&mysystem("ln -s $rootdir/$bfile.fam .") unless (-e "$bfile.fam");
    }
}


###################################
### GUESS BUILD
###################################

my @buigue_arr = ();
my $buigue_fini = 0;

unless (-e "$rootdir/buigue_done") {
    unless (-e "$rootdir/posing_done") {
		foreach (@bim_files) {
		    my $bfile = $_;
		    $bfile =~ s/.bim$//;
		    my $accfli ="$bfile".".hg19.bim";
	#	    print "he: $accfli\n";
	#	    exit;
		    unless (exists $bimhg19_array{$accfli}) {
				push @buigue_arr, "$buigue_script --lift19 $bfile.bim" ;#
		    }
		    else {
				$buigue_fini++;
		    }
		}

		if (@buigue_arr > 0) {
	    
		    $sjadir = $impute_dir;
		    $sjaname = "buigue";
	    
		    $sjatime = 2;
	#	    $sjatime = 4 if ($buigue_fini > 0);
	    
		    $sjamem = 3000;
		    @sjaarray = @buigue_arr;
	    
		    &send_jobarray;
		}
		else {
		    &mysystem ("touch $rootdir/buigue_done");
		    print "build_guess done\n";
		}
    }
}


###################################
### READREF
###################################



my @readref_arr = ();
my $readref_fini = 0;

unless (-e "$rootdir/readref_done") {
	foreach (@bim_files) {
	    my $bimfile = $_;
	    my $bfile = $bimfile;
	    $bfile =~ s/.bim$//;
	    my $accfli ="$bfile".".hg19.bim";
	    
	    
	    foreach my $chrloc(1..22) {
			my $bimref ="$accfli".".ref.chr$chrloc";
		    my $reffi = $reffile_struct;
			$reffi =~ s/###/$chrloc/g;
			unless (exists $bimref_array{$bimref}) {
				push @readref_arr, "$readref_script --ref $reffi --refheads id,NULL,position,a1,a0,$popname_uc --chr $chrloc $accfli" ;#
				#		print "$readref_script --chr $chrloc --ref $reffi $bimfile\n" ;#
			}
			else {
				$readref_fini++;
			}
	    }
	}
#    exit;
	if (@readref_arr > 0) {
	    
	    $sjadir = $impute_dir;
	    $sjaname = "readref";
	    $sjatime = 2;
	    $sjatime = 4 if ($readref_fini > 0);
	    
	    $sjamem = 3000;
	    $sjamaxpar = 100;
	    @sjaarray = @readref_arr;
	    


	    &send_jobarray;
	}
	else {
	    &mysystem ("touch $rootdir/readref_done");
	    print "readref done\n";
	}
}


###################################
### sum readref
###################################

unless (-e "$rootdir/readrefsum_done") {

	my @readrefsum_arr = ();
	my $readrefsum_fini = 0;
	
	unless (-e "$rootdir/readrefsum_done") {
	    foreach (@bim_files) {
		my $bimfile = $_;
		my $bfile = $bimfile;
		$bfile =~ s/.bim$//;
		my $accfli ="$bfile".".hg19.bim";
		my $bimref_done ="$accfli".".ref.sum.done";
#		print "looking for $bimref_done\n";
		unless (exists $bimref_array{$bimref_done}) {
		    push @readrefsum_arr, "$readrefsum_script $accfli" ;#
		}
		else {
		    $readrefsum_fini++;
		}
	    }
	    
	    if (@readrefsum_arr > 0) {

#		print "stragne\n";
#		exit;
		
		$sjadir = $impute_dir;
		$sjaname = "reresum";
		$sjatime = 2;
#		$sjatime = 4 if ($readrefsum_fini > 0);

		$sjamem = 1000;
		@sjaarray = @readrefsum_arr;
		
		&send_jobarray;
	    }
	    else {
		&mysystem ("touch $rootdir/readrefsum_done");
		print "readrefsum done\n";
	    }
	}
}


###################################
### CHECKPOS
##################################################
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### checkpos6 needs the var_chr_renaming (see my.readref) -> done
###################################

my @chepos_arr = ();
my $chepos_fini = 0;

#print "???\n";
#exit;

unless (-e "$rootdir/posing_done") {
    foreach (@bim_files) {
	my $bfile = $_;
	$bfile =~ s/.bim$//;
	my $locref = $bfile.".hg19.bim.ref.sum";
	my $accfli ="$bfile".".hg19.ch.bim";

	if (-e $locref) {
	    print "locref $locref is existing! safes some time\n";
	    unless (exists $bimpos_array{$accfli}) {
		push @chepos_arr, "$checkpos_script --dbcol 1,2,3 --dbsnp $rootdir/$impute_dir/$locref $bfile.hg19.bim" ;#
#		print "$checkpos_script --dbcol 1,2,3 --dbsnp $rootdir/$impute_dir/$locref $bfile.hg19.bim\n" ;#

	    }
	    else {
			$chepos_fini++;
	    }
	}
	else {
		die "$locref does not exist. Readref may have failed.\n"
#	    print "locref $locref is not existing! would be better if it did\n";
#	    unless (exists $bimpos_array{$accfli}) {
#		push @chepos_arr, "$checkpos_script --dbcol 1,8,9 --dbsnp $refdir/$suminfo_s $bfile.hg19.bim" ;#
#	    }
#	    else {
#		$chepos_fini++;
#	    }
	}

    }
#    exit;

   
    if (@chepos_arr > 0) {
	
	$sjadir = $impute_dir;
	$sjaname = "chepos";
	$sjatime = 2;
#	$sjatime = 4 if ($chepos_fini > 0);
	$sjamem = 3000;
	@sjaarray = @chepos_arr;
	
	&send_jobarray;
    }
    else {
		&mysystem ("touch $rootdir/posing_done");
		print "checkpos done\n";
    }
}

#exit;

###################################
### CHECKFLIP
###################################

#print "checkflip3?\n";
#exit;

my @chefli_arr = ();
my $chefli_fini = 0;

unless (-e "$rootdir/flipping_done") {
	foreach (@bim_files) {
	    my $bfile = $_;
	    $bfile =~ s/.bim$//;
	    my $accfli ="$bfile".".hg19.ch.fl.bim";
	    my $locref = $bfile.".hg19.bim.ref.sum";

	    if (-e $locref) {
			print "locref $locref is existing! safes some time\n";
			unless (exists $bimfli_array{$accfli}) {
				my $systmp = "$checkflip_script --dbcol 0,3,4,5 --fth $fth_th --sfh $sec_freq --info $rootdir/$impute_dir/$locref $bfile.hg19.ch.bim" ;
				push @chefli_arr, $systmp ;
#			    print "$systmp\n";
#			    exit;
#			    push @chepos_arr, "$checkpos_script --dbcol 1,2,3 --dbsnp $rootdir/$impute_dir/$locref $bfile.hg19.bim" ;#
		    
			}
			else {
		    	$chefli_fini++;
			}
	    }
	    else {
			die "$locref not found. This should not be possible; check readref and checkpos for errors.\n"
#			unless (exists $bimfli_array{$accfli}) {
#		    	push @chefli_arr, "$checkflip_script --fth $fth_th --sfh $sec_freq --info $refdir/$suminfo_s $bfile.hg19.ch.bim" ;
#			}
#			else {
#		    	$chefli_fini++;
#			}
	    }


	}

#	exit;

	if (@chefli_arr > 0) {
	    
	    $sjadir = $impute_dir;
	    $sjaname = "chefli";
	    $sjatime = 2;
#	    $sjatime = 4 if ($chefli_fini > 0);
	    $sjamem = 3000;
	    @sjaarray = @chefli_arr;
	    
	    &send_jobarray;
	    
	}
	else {
	    &mysystem ("touch $rootdir/flipping_done");
	    print "checkflip done\n";
	}
}




print "\n\n";
print "###############";
print "\n";
print "FINISHED!\n";
print "\n\n";

exit;
