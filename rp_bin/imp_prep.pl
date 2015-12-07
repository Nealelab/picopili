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
# - stop before prephasing/imputation (buigue/chepos/chefli/prep fam file only)
# - remove old version history
# - remove code for chunking/refind and beyond
# - remove commented out blocks
# - remove trioset, phase_txt, pcaer_dir
#
# TODO: remove unused args, dependencies
# - refine output format
# 
#### 


my $version = "1.0.24";
my $progname = $0;
$progname =~ s!^.*/!!;
my $command_line = "$progname @ARGV";



my $jnum = 7; ### number of imputation job per node

my $spliha_n = 1500; ## split haplotypes with N individuals

my $best_lahunt = 5;

my $multithread1 = 4;
my $multithread2 = 8;


my $phas = -1;




my $info_txt = "";
#my $homedir = "/home/gwas";
my $rootdir = "";

my $iname = "" ;


my $suminfo = "infosum_pos";
my $suminfo_n = "$suminfo.nsnps";
my $suminfo_r = "$suminfo.reffiles";
#my $suminfo_s = "$suminfo.sorted";
my $suminfo_s = "NAN";

my $job_bn_th = 1000;

my @ref_coll = ();



#my $hapmap_ref_root = "/home/gwas/pgc-samples/hapmap_ref/";


my $fth_th = 0.15;


use Sys::Hostname;
my $host = hostname;

#my $broad = 1 if ($ENV{"DOMAINNAME"} =~ /broadinstitute/);


#############################
# read config file
#############################

my $conf_file = $ENV{HOME}."/ricopili.conf";
my %conf = ();

die $!."($conf_file)" unless open FILE, "< $conf_file";
while (my $line = <FILE>){
    my @cells = split /\s+/, $line;
    $conf{$cells[0]} = $cells[1];
}
close FILE;

sub trans {
    my ($expr)=@_;
    unless (exists $conf{$expr}) {
	die "config file without entry: $expr\n";
    }
    $conf{$expr};
}

my $ploc = &trans("p2loc");
my $shloc = &trans("shloc"); # shapeit
my $hapmap_ref_root = &trans("hmloc");
my $homedir = &trans("home");
my $qloc = &trans("queue");
my $i2loc = &trans("i2loc");
my $liloc = &trans("liloc");
my $email = &trans("email");
my $loloc = &trans("loloc");


###############################################

$ref_coll[5] = "$hapmap_ref_root"."1KG/phased/subchr" ;
$ref_coll[6] = "$hapmap_ref_root"."1KG_june10/hapmap3_r2_plus_1000g_jun2010_b36_ceu/bgl/subchr_5" ;
$ref_coll[7] = "$hapmap_ref_root"."1KG_aug10/subchr" ;
$ref_coll[8] = "$hapmap_ref_root"."1KG_aug10_nodup/subchr" ;
$ref_coll[8883] = "$hapmap_ref_root"."1KG_aug10_nodup/mhc_window" ;
$ref_coll[8882] = "$hapmap_ref_root"."mars_window/1KG" ;

$ref_coll[88] = "$hapmap_ref_root"."1KG_phas1_umich/ref_0611/subchr" ;
$ref_coll[881] = "$hapmap_ref_root"."1KG_phas1_umich/ref_0611_eur/subchr" ;
$ref_coll[882] = "$hapmap_ref_root"."1KG_phas1_umich/ref_0611_eur_3Mb/subchr" ;
$ref_coll[8821] = "$hapmap_ref_root"."1KG_phas1_umich/ref_0611_eur_3Mb/subchr/test2" ;
$ref_coll[8822] = "$hapmap_ref_root"."1KG_phas1_umich/ref_0611_eur_3Mb/subchr/NOD2" ;
$ref_coll[88222] = "$hapmap_ref_root"."1KG_phas1_umich/ref_0611_eur_3Mb/subchr/best_basam" ;
$ref_coll[882222] = "$hapmap_ref_root"."1KG_phas1_umich/ref_0611_eur_3Mb/subchr/chr3_itih" ;
$ref_coll[882223] = "$hapmap_ref_root"."1KG_phas1_umich/ref_0611_eur_3Mb/subchr/chr19_ncan" ;

$ref_coll[9999] = "$hapmap_ref_root"."impute2_ref/1KG_Mar12/ALL_1000G_phase1integrated_feb2012_impute/test/subchr" ;
$ref_coll[9] = "$hapmap_ref_root"."impute2_ref/1KG_Mar12/ALL_1000G_phase1integrated_feb2012_impute/subchr" ;
$ref_coll[91] = "$hapmap_ref_root"."impute2_ref/1KG_Aug12/ALL_1000G_phase1integrated_v3_impute_macGT1/subchr" ;
$ref_coll[9111] = "$hapmap_ref_root"."impute2_ref/1KG_Aug12/ALL_1000G_phase1integrated_v3_impute_macGT1/subchr/test" ;
$ref_coll[9123] = "$hapmap_ref_root"."impute2_ref/1KG_Aug12/ALL_1000G_phase1integrated_v3_impute_macGT1/chr23/subchr" ;
$ref_coll[91231] = "$hapmap_ref_root"."impute2_ref/1KG_Aug12/ALL_1000G_phase1integrated_v3_impute_macGT1/chr23/subchr/test" ;


$ref_coll[923] = "$hapmap_ref_root"."impute2_ref/1KG_Mar12/ALL_1000G_phase1integrated_feb2012_impute/chr23/subchr" ;
$ref_coll[9231] = "$hapmap_ref_root"."impute2_ref/1KG_Mar12/ALL_1000G_phase1integrated_feb2012_impute/chr23_pseudo/subchr" ;


$ref_coll[555] = "$hapmap_ref_root"."1KG/phased/subchr/test";
$ref_coll[3] = "$hapmap_ref_root"."subchr";
$ref_coll[323] = "$hapmap_ref_root"."subchr/23";
$ref_coll[333] = "$hapmap_ref_root"."subchr/test" ;

$ref_coll[334] = "$hapmap_ref_root"."subchr/HLA" ;
$ref_coll[335] = "$hapmap_ref_root"."subchr/chr2_test" ;
$ref_coll[3333] = "$hapmap_ref_root"."subchr/test/local" ;
$ref_coll[3331] = "$hapmap_ref_root"."subchr/SDCCA" ;
$ref_coll[3332] = "$hapmap_ref_root"."mars_window" ;
$ref_coll[39] = "$hapmap_ref_root"."impute2_ref/HM3/hapmap3_r2_b36/subchr" ;
$ref_coll[399] = "$hapmap_ref_root"."impute2_ref/HM3/hapmap3_r2_b36/subchr.3mb" ;

$ref_coll[15] ="$hapmap_ref_root"."hla_t1d/subchr";
$ref_coll[152] ="$hapmap_ref_root"."impute2_ref/HLA_0813/orig/subchr";
$ref_coll[1521] ="$hapmap_ref_root"."impute2_ref/HLA_0813/orig/hg19";
$ref_coll[1522] ="$hapmap_ref_root"."mhc";
$ref_coll[515] = "$hapmap_ref_root"."1KG/phased/subchr/hla" ;
$ref_coll[511] = "$hapmap_ref_root"."1KG/phased/subchr/mir" ;
$ref_coll[512] = "$hapmap_ref_root"."1KG/phased/subchr/mir708" ;
$ref_coll[513] = "$hapmap_ref_root"."1KG/phased/subchr/HBII-108" ;

$ref_coll[514] = "$hapmap_ref_root"."1KG/phased/subchr/chr10_11" ;
$ref_coll[515] = "$hapmap_ref_root"."1KG/phased/subchr/chr8_9" ;
$ref_coll[516] = "$hapmap_ref_root"."1KG/phased/subchr/cacna1c" ;
$ref_coll[517] = "$hapmap_ref_root"."1KG/phased/subchr/csmd1" ;
$ref_coll[518] = "$hapmap_ref_root"."1KG/phased/subchr/tcf4" ;
$ref_coll[519] = "$hapmap_ref_root"."1KG/phased/subchr/chr2_20" ;
$ref_coll[520] = "$hapmap_ref_root"."1KG/phased/subchr/top_scz" ;
$ref_coll[521] = "$hapmap_ref_root"."1KG/phased/subchr/chr11_13" ;

$ref_coll[5555] ="$hapmap_ref_root"."impute2_ref/lboettger/chr16";

$ref_coll[4] = "$hapmap_ref_root"."CNV/subchr";
$ref_coll[444] = "$hapmap_ref_root"."CNV/subchr_test";

$ref_coll[2] = "$hapmap_ref_root"."phas2/subchr/outdir";
$ref_coll[222] = "$hapmap_ref_root"."phas2/subchr/outdir/chr12_1";
$ref_coll[311] = "$hapmap_ref_root"."hm3_ww/subchr";
$ref_coll[322] = "$hapmap_ref_root"."finref/fineur/refdan/subchr";


my $info_th = 0.1;
my $freq_th = 0.005;
my $bg_th = 0.8;

my $sjamem_incr = 0;
#exit;

my $sec_freq = .2;  ## secure freq distance around 50%






my $popname = "eur";

##### help message
my $usage = "
Usage : $progname [options] --phase PHASE --outname OUTNAME

version: $version


 --help            print this text and exits

 --phase INT       impute with HM - Phase INT as ref., no default; 
                       (mandatory if --refdir is not specified)

 --refdir STRING   full path of reference directory, overwrites --phase

 --outname STRING  identifier for imputation run (mandatory)



#### for trio datasets

 --triset STRING  for subset of trio datasets (can contain bimfiles)

 --spliha INT      split haplotypes with N individuals




##### alignment to reference:

  --popname STRING    important for freq-checks, either 
                           eur (default), 
                           asn (asian), 
                           amr (america), 
                           afr (africa), 
                           asw (african american)

  --sfh FLOAT         secure frequency around 50% (default: $sec_freq)
                                for checkflip (compare to reference),
                                only applied to AT/CG SNPs
  --fth FLOAT         frequency-diff to exclude SNPs, default $fth_th
                                for checkflip (compare to reference)



 
#### post - imputation cleaning

  --info_th FLOAT  threshold for infoscore, default = $info_th

  --freq_th FLOAT  threshold for frequence (cases and controls), default = $freq_th

  --bg_th FLOAT    threshold for frequence (cases and controls), default = $bg_th


### technical options

  --refiex            file containing refinds to exclude

  --sjamem_incr INT   increase all memore requests by INT Mb in steps of 1000 (1000 is 1Gb)

  --noclean           do not clean up intermediate files at the very end

  --force1            do not exit if same fail, but do this only once

  --sleep INT         sleep for INT seconds, try this if you think there is a race condition
                       (hints: stops a different steps, --serial works)
                       try using 5 seconds first.

  --serial            no sending jobs to queue all in one run
                          -> usually only used for testing  


### remarks 

  --phase is mandatory!! use --phase 0 for list of options

  --outname is mandatory!!




 created by Stephan Ripke 2009 at MGH, Boston, MA

";



use Getopt::Long;
GetOptions( 


    "sjamem_incr=i"=> \$sjamem_incr,
    "info_th=f"=> \$info_th,
    "freq_th=f"=> \$freq_th,
    "bg_th=f"=> \$bg_th,
    "triset=s"=> \my $trioset_file,

    "help"=> \my $help,
    "serial"=> \my $serial,
    "sleep=i"=> \my $sleep_sw,



    "outname=s"=> \my $outname,
    "refdir=s"=> \my $refdir_str,
    "phase=i"=> \ $phas,

    "sfh=f"=> \$sec_freq,
    "fth=f"=> \$fth_th,

    "spliha_n=i"=> \$spliha_n,
    "noclean"=> \my $noclean,
    "force1"=> \my $force1,


    "popname=s"=> \$popname,
    "refiex=s"=> \my $refiex_file,

    );





if ($sleep_sw) {
  print "sleeping for $sleep_sw seconds (only use if suspect of race condition)\n";
  sleep ($sleep_sw);
}




############################################################
## testing binaries
##############################################################
my @test_scripts;


my $readref_script = "my.readref";         ### my.pipeline_tar
my $readrefsum_script = "my.readref_sum";  ### my.pipeline_tar
my $buigue_script = "buigue";              ### my.pipeline_tar
my $checkpos_script = "checkpos6";         ### my.pipeline_tar
my $checkflip_script = "checkflip4";       ### my.pipeline_tar
my $chuck_script = "my.chuck2";            ### my.pipeline_tar
my $preph_script = "my.preph";             ### my.pipeline_tar
my $imp2_script = "my.imp2.3";             ### my.pipeline_tar
my $dos_script = "haps2dos4";              ### my.pipeline_tar
my $impprob_script = "impprob_to_2dos";    ### my.pipeline_tar
my $dabg_script = "daner_bg3";             ### my.pipeline_tar
my $cobg_script = "bcomb_5_p2";            ### my.pipeline_tar
my $cobg_gw_script = "comb_bg_dir_1";      ### my.pipeline_tar
my $prune_script = "my.prune";             ### my.pipeline_tar
my $merge_script = "my.merge";             ### my.pipeline_tar
my $pseudo_script = "haps2pseudo2";        ### my.pipeline_tar
my $lift_script = "lift18219";             ### my.pipeline_tar
my $trisha_script = "trio2shape";          ### my.pipeline_tar
my $splithap_script = "splithap_1";        ### my.pipeline_tar
my $cleandir_script = "my.cleandir";       ### my.pipeline_tar
my $cleanerrandout_script = "my.cleanerrandout";  ### my.pipeline_tar
my $pdflatex_script = "pdflatex";          ### my.pipeline_tar
my $mystart_script = "my.start_job";       ### my.pipeline_tar
my $mutt_script = "mutt";                  ### my.pipeline_tar
my $du_script = "my.du";                  ### my.pipeline_tar
my $blue_script = "blueprint";         ### my.pipeline_tar



push @test_scripts, $readref_script;
push @test_scripts, $impprob_script;
push @test_scripts, $readrefsum_script ;
push @test_scripts, $buigue_script ;
push @test_scripts, $checkpos_script ;
push @test_scripts, $checkflip_script ;
push @test_scripts, $chuck_script ;
push @test_scripts, $preph_script ;
push @test_scripts, $imp2_script ;
push @test_scripts, $dos_script ;
push @test_scripts, $dabg_script ;
push @test_scripts, $cobg_script ;
push @test_scripts, $cobg_gw_script ;
push @test_scripts, $prune_script ;
push @test_scripts, $merge_script ;
push @test_scripts, $pseudo_script ;
push @test_scripts, $lift_script ;
push @test_scripts, $trisha_script ;
push @test_scripts, $pdflatex_script ;
push @test_scripts, $splithap_script ;
push @test_scripts, $cleandir_script ;
push @test_scripts, $cleanerrandout_script ;
push @test_scripts, $du_script ;
push @test_scripts,  $mystart_script;
push @test_scripts,  $blue_script;

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
	push @miss_scripts, "cp /home/unix/sripke/bin/$scr_name ./\n";
	print "!!Error!! : No $scr_name command available\n" ;
    }
 
}



if (@miss_scripts > 0) {
  if (-e "get_scripts_on_broad.txt") {
    print "please remove this file and restart: get_scripts_on_broad.txt\n";
  }
  die $! unless open FILE1, "> get_scripts_on_broad.txt";
  foreach (@miss_scripts) {
    print FILE1 "$_";
  }
  close FILE1;


  print "exiting now -> have a look at get_scripts_on_broad.txt\n";
  exit;

}






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


print "....all necessary binaries found....\n";
print "------------------------------------\n";
#push @scripts,"id_tager_3";



#####################################
# "testing environment variable rp_perlpackages
####################################

print "testing environment variable rp_perlpackages....\n";
unless (exists $ENV{rp_perlpackages}) {
    print "Error: no environment variable for perl-packages, please re-install ricopili and make sure to follow all instructions\n";
    print "------------------------------------\n";
    exit;
}
print "....all set....\n";
print "------------------------------------\n";















my $nomega_sw = 1;



my $nomega = 0;
$nomega = 1 if ($nomega_sw);




die $usage if $help;

die $usage unless $outname;
if ($phas == -1) {
    unless ($refdir_str) {
	print "$usage\n";
	exit;
    }
}
if ($phas == 0) {
    unless ($refdir_str) {
	print "$usage\n";
	exit;
    }
}



#my ($xsnp,$xchr,$xbeg,$xend);
#($xsnp,$xchr,$xbeg,$xend)= split ',', $xareastr if ($xareastr);


if ($phas == 9) {
  print "please do not use old reference any more\n";
  exit;
}

my $p2_txt = "";
if ($phas == 2 ){
    $p2_txt = "--phase2";
}


my $refdir = "";

if ($refdir_str) {
    $refdir = $refdir_str;
}
else {
    $refdir = $ref_coll[$phas];
}

unless (-d $refdir) {
    print "reference directory ($refdir) is not existing\n";
    exit;
}


my $impute_dir = "pi_sub";


#my $postimp_dir = "$impute_dir/postimp_data";







sub fisher_yates_shuffle {
    my $deck = shift;  # $deck is a reference to an array
    my $i = @$deck;
    while ($i--) {
	my $j = int rand ($i+1);
	@$deck[$i,$j] = @$deck[$j,$i];
    }
}


#####################################
# print array to file
####################################

sub a2file {
    my ($file, @lines)=@_;
    die $! unless open FILE, "> $file";
    foreach (@lines){
	print FILE $_;
    }
    close FILE;
}


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

##########################################
# subroutine to split a plink-output-line with references
##########################################

sub split_line_ref {
    my ($line)=${$_[0]};
    chomp($line);
    $line =~ s/^[\s]+//g;
    my @cols=  split /\s+/, $line;
    \@cols;
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
# subroutine to count lines of a file
#####################################

sub count_lines {
    my ($file)=@_;
    my $lc=0;
    die "$file: ".$! unless open FILE, "< $file";
    while (<FILE>){
	$lc++;
    }
    close FILE;
    $lc;
}



#####################################
# subroutine to re-invoke this script
#####################################

sub reinvo_b {
    my ($message, $wt_file)=@_;
    my $now = localtime time;
    my $old_cmd = `tail -3 $loloc/impute_dir_info | head -1`;

    my $message_part = $info_txt."\t$message";
    $message = $info_txt."\t$message\t$now";

    &a2filenew_app("$loloc/impute_dir_info",$message);
    die "2 times already" if ($old_cmd =~ /$message_part/);
    chdir "$rootdir" or die "something strange";
    if ($qloc eq "bsub") {
	$wt_file =~ s/.*blueprint_joblist_file-//;
    }

    my $sys_re = "$blue_script --njob $job_bn_th -b \"$command_line\" --wa 4 --di -j --fwt $wt_file --na _if_$outname";
#    print "$sys_re\n";
    &mysystem ($sys_re);
    exit;

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


my $sjainfofile = "$loloc/impute_dir_info";
unless (-e $sjainfofile) {
    print "log-file ($sjainfofile) is not existing\n";
    print "please check loloc in ~/ricopili.conf\n";
    exit;
}
#my $sjainfofile = "$homedir/impute_dir_info_35_test";
my $sjainfotxt = "";
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


    if ($sjaname eq "finished") {

	my $fini_message ;
	$fini_message .= "\n\n##################################################################\n";
	$fini_message .= "##### CONGRATULATIONS: \n";
	$fini_message .= "##### rp_pipeline finished successfully:\n";
	$fini_message .= "##### $sjainfotxt\n";
	$fini_message .= "##### now start with PCA (see README in subdir pcaer_sub/)\n";
	$fini_message .= "##### or directly with postimputation analysis\n";
	$fini_message .= "##### have a look at the wiki page\n"; 
	$fini_message .= "##### https://sites.google.com/a/broadinstitute.org/ricopili/\n";
	$fini_message .= "##################################################################\n";
	print "$fini_message\n";

	
	die $! unless open SUC, "> success_file";
	print SUC $fini_message."\n";
	close SUC;

	&mysystem ('cat success_file | '.$mutt_script.' -s RP_pipeline_finished '.$email) ;

	my $sjarow      = $sjainfotxt."\t$sjaname\t$now";
	&a2filenew_app("$sjainfofile",$sjarow);


	exit;

    }


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
#    $message = $info_txt."\t$message\t$now";

    &a2filenew_app("$sjainfofile",$sjarow);

    if ($old_cmd =~ /$sjarow_part/){
	unless ($force1 ){
	    my $err_message ;
	    $err_message .= "##################################################################\n";
	    $err_message .= "##### Error: \n";
	    $err_message .= "##### step $sjaname has been done repeatedly without any progress\n";
	    $err_message .= "##### imputation pipeline stopped: $command_line\n";
	    $err_message .= "##### $sjainfotxt\n";
	    $err_message .= "##### if reason does not appear obvious\n";
	    $err_message .= "##### have a look at the wiki page\n"; 
	    $err_message .= "##### https://sites.google.com/a/broadinstitute.org/ricopili/\n";
	    $err_message .= "##### or contact the developers\n";
	    $err_message .= "##################################################################\n";
	    print "$err_message\n";

	    die $! unless open ERR, "> error_file";
	    print ERR $err_message."\n";
	    close ERR;


	    &mysystem ('cat error_file | '.$mutt_script.' -s RP_pipeline_error '.$email) ;

	    unless ($serial) {
		exit;
	    }

	}

    }


    $command_line =~ s/--force1//;


    my $wt_file = "$sjadir/blueprint_joblist_file-$sjaname.$outname";
    chdir "$rootdir" or die "something strange";
    if ($qloc eq "bsub") {
	$wt_file =~ s/.*blueprint_joblist_file-//;
    }

    if ($qloc eq "slurm") {
	$wt_file = "$sjadir/$jobfile.script.id";
    }

    if ($qloc eq "qsub") {
	$wt_file = "$sjadir/j.$sjaname.$outname.id";
    }
    if ($qloc eq "qsub_c") {
	$wt_file = "$sjadir/j.$sjaname.$outname.id";
    }
    if ($qloc eq "qsub_b") {
	$wt_file = "$sjadir/j.$sjaname.$outname.id";
    }
    


    if ($serial) {
	my $sys_re = "$command_line";
	&mysystem ($sys_re);
	exit;
    }
    else {
	my $sys_re = "$blue_script --njob $job_bn_th -b \"$command_line\" --wa 2 --di -j --fwt $wt_file --na _if_$outname";
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




#####################################
# subroutine to re-invoke this script
#####################################

sub reinvo_b_week {
    my ($message, $wt_file)=@_;
    my $now = localtime time;
    my $old_cmd = `tail -3 $loloc/impute_dir_info | head -1`;

    my $message_part = $info_txt."\t$message";
    $message = $info_txt."\t$message\t$now";

    &a2filenew_app("$loloc/impute_dir_info",$message);
    die "2 times already" if ($old_cmd =~ /$message_part/);
    chdir "$rootdir" or die "something strange";
    if ($qloc eq "bsub") {
	$wt_file =~ s/.*blueprint_joblist_file-//;
    }

    &mysystem ("$blue_script --week 1 --njob $job_bn_th -b \"$command_line\" --wa 10 --di -j --fwt $wt_file --na _if_$outname");
    exit;

}



##############################################
##############################################
#############  BEGIN
##############################################
##############################################


use Cwd;
use File::Path;
$rootdir = &Cwd::cwd();
$sjainfotxt = "$rootdir\t$command_line";






my $archive_dir = "/archive/gwas/scz/archive/$outname";

#exit;

unless (-e $impute_dir){
    print "impute_dir is not existing, create one for you\n";
    my @created = mkpath(   ## $created ?
			    $impute_dir,
			    {verbose => 0, mode => 0750},
	);
}



#####################################
## if new frequency file is existing
###################################
if (-e "$refdir/sumfrq.$popname") {
    $suminfo_s = "sumfrq.$popname";
}
else {

    my $popname_uc = uc($popname);
    if (-e "$refdir/sumfrq.$popname_uc") {
	$suminfo_s = "sumfrq.$popname_uc";
    }
    else {
	print "$refdir/sumfrq.$popname_uc in refdir is not existing!!!\n";
	die;
#	sleep(10);
    }
}

#print $refdir."\n";
#print $suminfo_s."\n";
#exit;


unless (-e "$refdir/$suminfo_n"){
    print "ERROR: $refdir/$suminfo_n not existing\n";

#    chdir ($refdir);
#    &mysystem ("wc -l *.info_pos > $suminfo_n");
#    chdir ($rootdir);
}


unless (-e "$refdir/$suminfo_r"){
    print "ERROR: $refdir/$suminfo_r not existing\n";
    die;
#    chdir ($refdir);
#    &mysystem ("ls sc_*.bgl > $suminfo_r");
#    chdir ($rootdir);
}

#my @refallfiles;

#    opendir(DIR, "$refdir") || die "can't opendir .: $!";
#    @refallfiles = readdir(DIR);
#    closedir DIR;
#}
my $cc=0;



my %refiex;
if ($refiex_file) {
    print "read $refiex_file\n";
    die $!." <$refiex_file>" unless open IN, "< $refiex_file";
    while (my $line = <IN>){
	chomp($line);
	$refiex{$line} = 1;
	print "$line\n";
    }
    close IN;

}







my @reffiles;
print "read $refdir/$suminfo_n....\n";
die $!." <$refdir/$suminfo_n>" unless open IN, "< $refdir/$suminfo_n";
while (my $line = <IN>){
    my @cells = &split_line($line);
    die "problem with $refdir/$suminfo_n" if (@cells < 2);
    my $bgl_file = $cells[1];
    $bgl_file =~ s/.info_pos$//;

    if ($refiex_file) {

	my $refind = $bgl_file;
	if ($refind =~ /chr[0-9]*_[0-9]*_[0-9]*/){
	    $refind =~ s/.*(chr[0-9]*_[0-9]*_[0-9]*).*/\1/;
	}
	else {
	    $refind =~ s/.*(chr[0-9]*_[0-9]*).*/\1/;
	}

	if (exists $refiex{$refind}){
	    print "exclude: $bgl_file\n" ;
	    next;
	}

    }

    next if ($bgl_file eq "total");

#    print "$bgl_file\n";
    push @reffiles, $bgl_file;
    $cc++;
}
close IN;
#print "finished reading $refdir/$suminfo_n\n";
die "reference directory <$refdir> empty (no sc.*bgl)" if (@reffiles == 0);



#exit;
#print "sleep\n";
#sleep(5);


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
my @bim_files = grep {/bim$/} @files;
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


################################################
### set info text
####################################################


$info_txt = "command:\t\"$command_line\"\tdir:\t$rootdir";



#####################################################
## check readref
########################################################


my $readref_sw = 1;

foreach my $chrloc(1..22) {
    my $reffi ="$refdir/$suminfo_s.$chrloc.gz";
    unless (-e $reffi) {
	print "Warning: $reffi is not existing, it's ok if using older reference\n";
	$readref_sw = 0;
	last;
    }
}

if ($readref_sw == 1) {
    print "efficient reference alignment switched on\n";
}
else {
    print "efficient reference alignment switched off, please check refdir, will continue in 3 sec...\n";
    sleep(3);
}
#print "exit;\n";
#exit;



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

if ($readref_sw == 1) {
    unless (-e "$rootdir/readref_done") {
	foreach (@bim_files) {
	    my $bimfile = $_;
	    my $bfile = $bimfile;
	    $bfile =~ s/.bim$//;
	    my $accfli ="$bfile".".hg19.bim";

	    
	    foreach my $chrloc(1..22) {
		my $bimref ="$accfli".".ref.chr$chrloc";
		my $reffi ="$refdir/$suminfo_s.$chrloc.gz";
		unless (exists $bimref_array{$bimref}) {
		    push @readref_arr, "$readref_script --chr $chrloc --ref $reffi $accfli" ;#
#		print "$readref_script --chr $chrloc --ref $reffi $bimfile\n" ;#
		}
		else {
		    $readref_fini;
		}
	    }
	}
#    exit;
	if (@readref_arr > 0) {
	    
	    $sjadir = $impute_dir;
	    $sjaname = "readref";
	    $sjatime = 2;
#	    $sjatime = 4 if ($readref_fini > 0);
	    
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
}


#exit;



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
	    print "locref $locref is not existing! would be better if it did\n";
	    unless (exists $bimpos_array{$accfli}) {
		push @chepos_arr, "$checkpos_script --dbcol 1,8,9 --dbsnp $refdir/$suminfo_s $bfile.hg19.bim" ;#
	    }
	    else {
		$chepos_fini++;
	    }
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
	    my $locref = $bfile.".bim.ref.sum";

	    if (-e $locref) {
		print "locref $locref is existing! safes some time\n";
		unless (exists $bimfli_array{$accfli}) {
		    my $systmp = "$checkflip_script --dbcol 0,3,4,5 --fth $fth_th --sfh $sec_freq --info $rootdir/$impute_dir/$locref $bfile.hg19.ch.bim" ;
		    push @chefli_arr, $systmp ;
#		    print "$systmp\n";
#		    exit;
#		    push @chepos_arr, "$checkpos_script --dbcol 1,2,3 --dbsnp $rootdir/$impute_dir/$locref $bfile.hg19.bim" ;#
		    
		}
		else {
		    $chefli_fini++;
		}
	    }
	    else {
		unless (exists $bimfli_array{$accfli}) {
		    push @chefli_arr, "$checkflip_script --fth $fth_th --sfh $sec_freq --info $refdir/$suminfo_s $bfile.hg19.ch.bim" ;
		}
		else {
		    $chefli_fini++;
		}
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



# exit;


###########################
#### here preparation of famfiles for shapeit
############################


print "prepare famfiles for shapeit\n";

unless (-e "$rootdir/puting_done") {
    foreach (@bim_files) {

	my $bprefix = $_;
	$bprefix =~ s/.bim$//;
	my %sex_hash = ();

	### include sex check for chrX
	if ($phas == 923 || $phas == 9123 || $phas == 91231) {
	    if (-e "$bprefix.hg19.ch.fl.fam") {
		unless (-e "$bprefix.hg19.ch.fl.sexcheck") {
		    my $sx = "$ploc/plink --memory 2000  --bfile $bprefix.hg19.ch.fl --check-sex --out $bprefix.hg19.ch.fl";
		    &mysystem ($sx);
		}
		
		
		die $! unless open SI, "< $bprefix.hg19.ch.fl.sexcheck";
		while (my $line = <SI>){
		    my @cells = &split_line($line);
		    if ($cells[5] < .5) {
			$sex_hash{"$cells[0]\t$cells[1]"} = 2;
		    }
		    else {
			$sex_hash{"$cells[0]\t$cells[1]"} = 1;
		    }
		}
		close SI;
	    }
	}
	
	if (-e "$bprefix.hg19.ch.fl.fam"){ 
	    unless (-e "$bprefix.hg19.ch.fl.fam.idnum") {
		die $! unless open FI, "< $bprefix.fam";
		die $! unless open FO, "> $bprefix.hg19.ch.fl.fam.idnum.tmp";
		die $! unless open FT, "> $bprefix.hg19.ch.fl.fam.transl";
		my $cc = 1;
		while (my $line = <FI>){
		    my @cells = &split_line($line);
		    print FO "$cc $cc"; 

		    print FO " 0"; 
		    print FO " 0"; 

#		    print FO " $cells[2]"; 
#		    print FO " $cells[3]"; 
		    if (exists $sex_hash{"$cells[0]\t$cells[1]"}){
			print FO " ".$sex_hash{"$cells[0]\t$cells[1]"}; 
		    }
		    else {
			print FO " $cells[4]"; 
			if ($phas == 923 || $phas == 9123 || $phas == 91231) {
			    print "Error: no sex-check on X-chr?\n";
			    die;
			}
		    }
		    
		    print FO " $cells[5]\n";
		    print FT "$cc"; 
		    print FT " $cells[0]";
		    print FT " $cells[1]\n";
		    $cc++;
		}
		close FI;
		close FO;
		close FT;
		my $nloc = $cc -1;
		die $! unless open FN, "> $bprefix.hg19.ch.fl.fam.n";
		print FN $nloc."\n";
		close FN;
		&mysystem ("mv $bprefix.hg19.ch.fl.fam.idnum.tmp $bprefix.hg19.ch.fl.fam.idnum");
	    }
	}
    }
}

print "\n\n"
print "###############"
print "\n";
print "FINISHED!\n"
print "\n\n"

exit;
