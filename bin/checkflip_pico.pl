#!/usr/bin/env perl
use strict;



###########################################################################################################
#
#
#    checkerflip
#
#          created by Stephan Ripke, Broadinstitute, sripke@broadinstitute.org
#
#                                  01/14/10
#
#          Adapted for Picopili by Raymond Walters, rwalters(at)broadinstitute.org
#
#
#
#    checks alelles of a bim-file (plink-binary-dataset) with reference-info (created with refinfo)
#
#
#
#I converted all datasets to the same format (name, frequency, etc.) as the reference (here HM3):
#- first I check, if the same rs-name has the same position-information. if bp-position differs, I correct it, but I exclude with differing chromosome-information. I exclude SNPs not found in reference (around 20K on Chr.X) or multi-occurent in dataset.
#- then I check alleles: all unambigous SNPs are easy, either flipped, not-flipped, or to-be-excluded (e.g. AC on reference and AG in dataset)
#- for ambiguous SNPs (AT or CG) I check on frequency and exclude all SNPs with MAF > .4. SNPs with MAF <= .4 are matched on the rarer allele.
#- I exclude all SNPs with a frequency missmatch of more than 15%
#
#
##########################################################################################################
# v4: with columns from reference file

#############################

my $conf_file = $ENV{HOME}."/picopili.conf";
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

my $sloc = &trans("sloc");
my $p2loc = &trans("p2loc");


#######################################

my $version = "1.0.0";
my $progname = $0;
$progname =~ s!^.*/!!;

my $dscol = 0;  ## snp-col in reference
my $da1col = 2;  ## chr-col in reference
my $da2col = 4;  ## chr-col in reference
my $dfcol = 3;  ## chr-col in reference


my $info_file = "";

my $refdir = "";

my $frq_th = .15;
my $subdir = "flip_subdir";

my $sec_freq = .2;  ## secure freq distance around 50%


my $usage = "
Usage : $progname [options] bim-file1 bim-file2 ......

version: $version

  --refdir STRING     location of reference-directory, default $refdir
  --ploc STRING       location of plink-binary (default is found from picopili.conf)
                       default: $p2loc
  --info STRING       other info-file (absolute path) -> overwrites --refdir
  --subdir STRING     subdir, to put end-dataset into, default: $subdir
  --pos               first translate SNP names into positions
  --fth INT           frequency difference threshold to exclude SNPs, default: $frq_th
  --help              print this message and exit

  --dbcol INT,INT,INT,INT snp-col,a1-col,a2-col,Freq_a1_col in dbsnp-file: default: $dscol,$da1col,$da2col,$dfcol 

#chr10_100002841_I	758	I2	0.065	D	0.935	48	10	100002841
#0                         1      2        3      4       

#0 2 4 3


  --sfh FLOAT         secure frequency around 50% (default: $sec_freq)
                        -> .5 +- .5 * sfh

#  --replace           replace old dataset with new one

";

use Cwd;
use File::Path;
use Getopt::Long;
GetOptions( 
    "refdir=s"=> \$refdir,
    "ploc=s"=> \$p2loc,
    "info=s"=> \$info_file,
    "fth=f"=> \$frq_th,
    "sfh=f"=> \$sec_freq,
    "help"=> \my $help,
    "pos"=> \my $pos_tr,
    "dbcol=s"=> \my $dcolstr,
#    "replace"=> \my $replace,
    "subdir=s"=> \$subdir,
    );

die "$usage\n" if ($help);

if ($info_file eq "") {
    die "$usage\n";
}
else {
    die "couldn't find info-file $info_file" unless (-e $info_file);
}

if ($dcolstr) {
    ($dscol,$da1col,$da2col,$dfcol) = split(/,/, $dcolstr);
    die "*****\nwrong col-usage\n****\n$usage\n" if ($dfcol eq "");
}


#############################
# test, if plink is present
#############################


unless (-e "$p2loc/plink" ){
    print "\n***** Error: couldn't find the following:\n";
    print "please check --ploc or $conf_file\n";
    exit;
}


my $bim_file=@ARGV[0];



my $bim_xal=$bim_file.".xal";
my $bim_inde=$bim_file.".inde";
my $bim_gf=$bim_file.".gf";
my $bim_bf=$bim_file.".bf";
my $bim_uif=$bim_file.".uif";
my $bim_usf=$bim_file.".usf";
my $bim_fli=$bim_file.".fli";
my $bim_nal=$bim_file.".nal";


my $bim_excl=$bim_file.".excl";
my $bim_log=$bim_file.".report";



my $bfile = $bim_file;
die "please bim_file ($bfile)" unless ($bfile =~ /.bim$/);
$bfile =~ s/.bim$//;


die "Error: bim_file ($bfile) does not exist, please check command" unless (-e "$bfile.bim");
die "Error: bed_file ($bfile) does not exist, please check command" unless (-e "$bfile.bed");
die "Error: fam_file ($bfile) does not exist, please check command" unless (-e "$bfile.fam");




my $bfile_flipped = $bfile.".fl";
my $frq_file = "$bfile.frq";
my $frq_sorted = "$bfile.frq_sorted";
my $frq_joined = "$bfile.frq_joined";



###################################################
###  system call with test if successfull
###################################################
my @cmd_collect;

sub mysystem(){
    my ($systemstr)="@_";
    system($systemstr);
    my $status = ($? >> 8);
    die "$systemstr\n->system call failed: $status" if ($status != 0);
    push @cmd_collect, $systemstr;

}



##########################################
# split a whitespace line
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

###############################
### bring alleles in alphabetical order genotype

sub al2gt {
    my ($a1,$a2) = @_;
    my $gt;
    if ($a1 lt $a2) {
	$gt = $a1.$a2;
    }
    else {
	$gt = $a2.$a1;
    }
}

sub al2gt_str {
    my ($a1,$a2) = @_;
    my $gt;
    if ($a1 lt $a2) {
	$gt = $a1."/".$a2;
    }
    else {
	$gt = $a2."/".$a1;
    }
}

###############################
### flip SNP
my %tr =();
$tr{"A"}="T";
$tr{"C"}="G";
$tr{"G"}="C";
$tr{"T"}="A";

sub flip {
    my ($a1,$a2,$fr) = @_;
    my $fro = $fr;
    my $a1o = $tr{$a1};
    my $a2o = $tr{$a2};
    ($a1o,$a2o,$fro) ;
}




###############################
### turn SNP

sub turn {
    my ($a1,$a2,$fr) = @_;
    my $fro = 1 - $fr;
    ($a2,$a1,$fro) ;
}

###############################
###   BEGIN
######################

my @log_lines;
push @log_lines,  "\nsummary checkflip, version: $version\n\n";

use File::Copy;
use File::Path;
use Cwd;


my $rootdir = &Cwd::cwd();


my $workdir = "$sloc/flip_$bim_file";

while (-e $workdir) {
    $workdir .= ".f";
}

print "workdir:\t$workdir\n";

my @created = mkpath(   ## $created ?
			$workdir,
			{verbose => 0, mode => 0750},
    );

chdir ($workdir);

my $bfile = $bim_file;
$bfile =~ s/.bim$//;

&mysystem("ln -sf $rootdir/$bfile.bim .") unless (-e "$bfile.bim");
&mysystem("ln -sf $rootdir/$bfile.bed .") unless (-e "$bfile.bed");
&mysystem("ln -sf $rootdir/$bfile.fam .") unless (-e "$bfile.fam");



###############################
### turn on pos_tr, if needed, now switched off, this is outdated

die $!." <$info_file>" unless open REFI, "< $info_file";
my $line = <REFI>;
$line = <REFI>;
my @cells = &split_line($line);
#$pos_tr = 1 if ($cells[0] =~ /[0-9]*:[0-9]*/);
close REFI;



if ($pos_tr) {
    &mysystem ("cp $bim_file $bim_file.sich") unless (-e "$bim_file.sich");
    die $!." <$bim_file.sich>" unless open BIMI, "< $bim_file.sich";
    die $!." <$bim_file>" unless open BIMO, "> $bim_file";

    while (my $line = <BIMI>) {
	my @cells = &split_line($line);
	$cells[1] = "$cells[0]:$cells[3]";
	print BIMO "@cells\n";
    }
    close BIMI;
    close BIMO;
}

#exit;




####### create freq-file
unless (-e $frq_file) {
    print "creating frq-file\n";
    &mysystem ("$p2loc/plink --memory 2000 --bfile $bfile --out $bfile --freq --nonfounders> /dev/null");
}



### read bim-file
#print "checker, version $version\n";
print "read bim-file\n";
my %pos;
my %gt;
my %a1;
my %a2;
my %gt_name;
my %frqa1;


if (0) {
die $!." <$bim_file>" unless open BIMI, "< $bim_file";

while (my $line = <BIMI>) {
    my @cells = &split_line($line);
    $pos{$cells[1]} = "$cells[0]:$cells[3]";
#    $pos{$cells[1]} = "$cells[1]" if ($pos);
    my $gt;
    my $a1 = $cells[4];
    my $a2 = $cells[5];
    if ($a1 lt $a2) {
	$gt = $a1.$a2;
    }
    else {
	$gt = $a2.$a1;
    }
    $gt{$pos{$cells[1]}} = $gt;
    $gt_name{$pos{$cells[1]}} = $cells[1];
}
close BIMI;
}

### read frq-file


print "read frq-file\n";
my %frq_info;

die "$frq_file".$! unless open FILE , "< $frq_file";
while (my $line = <FILE>){
    my @cells = @{&split_line_ref(\$line)};
    $frq_info{$cells[1]} = " $cells[2] $cells[3] $cells[4]";
}
close FILE;


print "read reference-file\n";
die "$info_file".$! unless open FILE , "< $info_file";
die "$frq_joined".$! unless open OF , "> $frq_joined";
while (my $line = <FILE>){
    my @cells = @{&split_line_ref(\$line)};
    if (exists $frq_info{$cells[$dscol]}) {
#	my $info_loc = " $cells[2] $cells[4] $cells[3]";
	my $info_loc = " $cells[$da1col] $cells[$da2col] $cells[$dfcol]";
	print OF $cells[$dscol].$frq_info{$cells[0]}.$info_loc."\n";
	delete ($frq_info{$cells[0]});
    }
}
close FILE;

print "print remainder\n";
foreach my $snp_loc (keys %frq_info) {
    print OF $snp_loc.$frq_info{$snp_loc}."\n";
}
close OF;






#### depreciated
if (0) {
print "sort frq-file\n";
&mysystem("sort -b -k 2,2 $frq_file > $frq_sorted");
print "join with reference-db\n";
&mysystem("join -o 0,1.3,1.4,1.5,2.3,2.5,2.4 -a 1 -1 2 $frq_sorted $info_file > $frq_joined");

}







if (0) {
    die $!." <$frq_file>" unless open FRQ, "< $frq_file";
    my $header = <FRQ>;
    while (my $line = <FRQ>) {
	my @cells = &split_line($line);
#	$pos{$cells[1]} = "$cells[0]:$cells[3]";
#	$pos{$cells[1]} = "$cells[1]" if ($nopos);
	$pos{$cells[1]} = $cells[1];
	$a1{$cells[1]} = $cells[2];
	$a2{$cells[1]} = $cells[3];
	$frqa1{$cells[1]} = $cells[4];

    }
    close FRQ;
}


#exit;

###################################
### compare to allele_reference
##################################
print "compare to strand reference\n";
#### flipped genotypes
my %trab =();
$trab{"AC"}="GT";
$trab{"AG"}="CT";
$trab{"CT"}="AG";
$trab{"GT"}="AC";

#### unflippable genotypes
my %xtr =();
$xtr{"AT"}=1;
$xtr{"CG"}=1;

die $!." <$frq_joined>" unless open FJ, "< $frq_joined";

my @snp_bf;  ## snps with good alleles, bad freq
push @snp_bf, "snp\ta1_bim\ta2_bim\ta1_ref\ta2_ref\tf_bim\tf_ref\tf_diff\n";
my @snp_gf;  ## snps with good alleles, good freq

my @snp_fg;  ## snps with flippable alleles, matched
my @snp_ff;  ## snps with flippable alleles, flipped

my @snp_uif;  ## snps with unflippable alleles, insecure freq
my @snp_usf;  ## snps with unflippable alleles, secure freq
my @snp_xal;  ## snps with bad alleles
my @snp_inde;  ## snps with insertion-deletion translation
my @snp_fli_am;  ## snps to be flipped, ambigous
my @snp_fli_un;  ## snps to be flipped, unambigous


my @snp_nan;  ## snps not annotated

my $sfh_lo = .5 - ($sec_freq / 2);  ### half of sec_freq
my $sfh_hi = .5 + ($sec_freq / 2);  ### half of sec_freq

my $cc=0;
while (my $line = <FJ>) {
    $cc++;
    print "$cc entries read\n" if ($cc % 1000000 == 0);
    my @cells = &split_line($line);

    my $ncells = @cells;

    next if ($cells[0] eq "SNP"); ## skip header

    if ($ncells < 7){
	push @snp_nan, "$line";
	next;
    }

    ### process reference information
    my $snp = $cells[0];
    my $a1_bim = $cells[1];
    my $a2_bim = $cells[2];
    my $fr_bim = $cells[3];
    my $a1_ref = $cells[4];
    my $a2_ref = $cells[5];
    my $fr_ref = $cells[6];



    if ($a1_bim eq "I" || $a2_bim eq "I") {

	my $a1_bim_orig = $a1_bim;
	my $a2_bim_orig = $a2_bim;

	my $la1r = length($a1_ref);
	my $la2r = length($a2_ref);

	if ($a1_bim eq "I") {
	    if ($la1r > $la2r) {
		$a1_bim = $a1_ref;
		$a2_bim = $a2_ref;
	    }
	    else {
		$a1_bim = $a2_ref;
		$a2_bim = $a1_ref;
	    }
	}
	if ($a2_bim eq "I") {
	    if ($la1r > $la2r) {
		$a2_bim = $a1_ref;
		$a1_bim = $a2_ref;
	    }
	    else {
		$a2_bim = $a2_ref;
		$a1_bim = $a1_ref;
	    }
	}
	
	
	push @snp_inde, "$snp\t$a1_bim_orig\t$a2_bim_orig\t$a1_bim\t$a2_bim\t$fr_bim\t$fr_ref\n";
    }


   if ($a1_bim eq "-" || $a2_bim eq "-") {

	my $a1_bim_orig = $a1_bim;
	my $a2_bim_orig = $a2_bim;

	my $la1r = length($a1_ref);
	my $la2r = length($a2_ref);

	### the other way round, - is the deletion not the insertion, so changed ep into ne
	if ($a1_bim ne "-") {
	    if ($la1r > $la2r) {
		$a1_bim = $a1_ref;
		$a2_bim = $a2_ref;
	    }
	    else {
		$a1_bim = $a2_ref;
		$a2_bim = $a1_ref;
	    }
	}
	if ($a2_bim ne "-") {
	    if ($la1r > $la2r) {
		$a2_bim = $a1_ref;
		$a1_bim = $a2_ref;
	    }
	    else {
		$a2_bim = $a2_ref;
		$a1_bim = $a1_ref;
	    }
	}
	
	
	push @snp_inde, "$snp\t$a1_bim_orig\t$a2_bim_orig\t$a1_bim\t$a2_bim\t$fr_bim\t$fr_ref\n";
#	print "$snp\t$a1_bim_orig\t$a2_bim_orig\t$a1_bim\t$a2_bim\t$fr_bim\t$fr_ref\n";


    }





    

    my $gt_ref = &al2gt ($a1_ref,$a2_ref);
    my $gt_ref_str = &al2gt_str ($a1_ref,$a2_ref);
    my $gt_bim = &al2gt ($a1_bim,$a2_bim);
    my $gt_bim_str = &al2gt_str ($a1_bim,$a2_bim);


    my $frq_diff = 0;

    my $nomatch = 0;


    ## unflippable
    if (exists $xtr{$gt_ref} || exists $xtr{$gt_bim}){
	## must be the same, otherwise non-matching
	if ( $gt_ref ne $gt_bim){
	    if (exists $xtr{$gt_ref} && exists $xtr{$gt_bim}){
		push @snp_xal, "$snp\t$gt_bim_str\t$gt_ref_str\tboth_are_ambiguous\n";
		$nomatch = 1;
	    }
	    else {
		push @snp_xal, "$snp\t$gt_bim_str\t$gt_ref_str\tone_is_ambiguous\n";
		$nomatch = 1;
	    }
#	    print "$snp\t$gt_bim_str\t$gt_ref_str\n";
	}
 	else {
	    ## flip, if Flipdiff AND same allele
	    if (abs ($fr_bim - $fr_ref) > abs ((1-$fr_bim)-$fr_ref)){
		## if match, then flip
		if ($a1_bim eq $a1_ref && $a2_bim eq $a2_ref) {
		    push @snp_fli_am, "$snp\t$a1_bim\t$a2_bim\t$a1_ref\t$a2_ref\t$fr_bim\t$fr_ref\n";
		    ($a1_bim,$a2_bim,$fr_bim) = &turn($a1_bim,$a2_bim,$fr_bim);
		}
		else {
		    ($a1_bim,$a2_bim,$fr_bim) = &turn($a1_bim,$a2_bim,$fr_bim);
		}
	    }
	    ## flip, if no Flipdiff AND other allele
	    else {
		## if match, then flip
		if ($a1_bim eq $a2_ref && $a2_bim eq $a1_ref) {
		    push @snp_fli_am, "$snp\t$a1_bim\t$a2_bim\t$a1_ref\t$a2_ref\t$fr_bim\t$fr_ref\n";
		    ($a1_bim,$a2_bim,$fr_bim) = &flip($a1_bim,$a2_bim,$fr_bim);
		}
	    }


	    ## insecure freq
	    if (($fr_bim > $sfh_lo && $fr_bim < $sfh_hi )|| ($fr_ref > $sfh_lo && $fr_ref < $sfh_hi)){
		push @snp_uif, "$snp\t$a1_bim\t$a2_bim\t$a1_ref\t$a2_ref\t$fr_bim\t$fr_ref\n";
	    }
	    ## secure freq
	    else {
		push @snp_usf, "$snp\t$a1_bim\t$a2_bim\t$a1_ref\t$a2_ref\t$fr_bim\t$fr_ref\n";
	    }

	}
    }
    else {

	## flip if necessary
	if ($trab{$gt_bim} eq $gt_ref){
#	    printf "$snp\t$a1_bim\t$a2_bim\t$gt_bim\t$gt_ref\t%.4f\t%.4f\n",$fr_bim,$fr_ref;
	    push @snp_fli_un, "$snp\t$a1_bim\t$a2_bim\t$a1_ref\t$a2_ref\t$fr_bim\t$fr_ref\n";
	    ($a1_bim,$a2_bim,$fr_bim) = &flip ($a1_bim,$a2_bim,$fr_bim);
#	    printf "$snp\t$a1_bim\t$a2_bim\t$gt_bim\t$gt_ref\t%.4f\t%.4f\n",$fr_bim,$fr_ref;
	}

	## if match
	if ($a1_bim eq $a1_ref && $a2_bim eq $a2_ref) {
	    push @snp_fg, "$snp\t$gt_bim\t$gt_ref\t$fr_bim\t$fr_ref\n";
	}
	## turn if necessary
	elsif ($a1_bim eq $a2_ref && $a2_bim eq $a1_ref) {
#	    printf "$snp\t$a1_bim\t$a2_bim\t$gt_bim\t$gt_ref\t%.4f\t%.4f\n",$fr_bim,$fr_ref;
	    ($a1_bim,$a2_bim,$fr_bim) = &turn ($a1_bim,$a2_bim,$fr_bim);
	    push @snp_fg, "$snp\t$a1_bim\t$a2_bim\t$a1_ref\t$a2_ref\t$fr_bim\t$fr_ref\n";
#	    printf "$snp\t$a1_bim\t$a2_bim\t$gt_bim\t$gt_ref\t%.4f\t%.4f\n",$fr_bim,$fr_ref;
	}
	
	## no match
	else {
#	    push @snp_xal, "$snp\t$gt_bim\t$gt_ref\n";
	    $nomatch =1;
	    push @snp_xal, "$snp\t$gt_bim_str\t$gt_ref_str\tunambiguous\n";
	}
#	printf "$snp\t$gt_bim\t$gt_ref\t%.4f\t%.4f\n",$fr_bim,$fr_ref;
    }

    ## report freq_violation
    $frq_diff = abs ($fr_bim - $fr_ref);
    if ($frq_diff < $frq_th){
	if ($nomatch == 0) {
	    push @snp_gf, sprintf "$snp\t$a1_bim\t$a2_bim\t$a1_ref\t$a2_ref\t%.4f\t%.4f\t%.4f\n",$fr_bim,$fr_ref,$frq_diff;
	}
    }
    else {
	push @snp_bf, sprintf "$snp\t$a1_bim\t$a2_bim\t$a1_ref\t$a2_ref\t%.4f\t%.4f\t%.4f\n",$fr_bim,$fr_ref,$frq_diff;
    }




 
}
close FJ;


my $nsnp_bf = @snp_bf -1;

&a2file($bim_xal,@snp_xal);
&a2file($bim_inde,@snp_inde);
&a2file($bim_gf,@snp_gf);
&a2file($bim_bf,@snp_bf);
&a2file($bim_fli,@snp_fli_am,@snp_fli_un);
&a2file($bim_uif,@snp_uif);
&a2file($bim_usf,@snp_usf);
&a2file($bim_nal,@snp_nan);


push @log_lines,  "\ncomparsion to allele-db:\n";
push @log_lines,  "\nkeep:\n";
push @log_lines,  "allele match:              \t".@snp_gf."\t-> $bim_gf\n";
push @log_lines,  "flipped, ambiguous:        \t".@snp_fli_am."\t-> $bim_fli\n";
push @log_lines,  "flipped, unambiguous:      \t".@snp_fli_un."\t-> $bim_fli\n";
push @log_lines,  "ambiguous, secure freq:    \t".@snp_usf."\t-> $bim_usf\n\n";


push @log_lines,  "\nexclude:\n";
push @log_lines,  "no match in alleles:       \t".@snp_xal."\t-> $bim_xal\n";
push @log_lines,  "not found:                 \t".@snp_nan."   \t-> $bim_nal\n";
push @log_lines,  "ambiguous, insecure freq:  \t".@snp_uif."\t-> $bim_uif\n";
push @log_lines,  "freq-diff > $frq_th   :       \t".$nsnp_bf."\t-> $bim_bf\n";

push @log_lines,  "\nexclude:\n";
push @log_lines,  "translated indels:       \t".@snp_inde."\t-> $bim_inde\n";

push @log_lines,  "\n\n\nI would exclude unmatched, unfound, different and insecure unflippable SNPs:\n";
push @log_lines,  "\nthese are roughly ".@snp_nan + @snp_xal + @snp_uif + $nsnp_bf." SNPs (multiple entries possible):\n";

my $cmd1 = "cat $bim_xal $bim_uif $bim_bf $bim_nal > $bim_excl\n";
push @log_lines,  $cmd1;
push @log_lines,  "use plink\n";
#my $cmd2 = "$p2loc/plink --bfile $bfile --exclude $bim_excl --flip  $bim_fli --out $bfile_flipped --make-bed > /dev/null\n";
my $cmd2 = "$p2loc/plink --memory 2000 --update-alleles $bim_inde --bfile $bfile --exclude $bim_excl --flip  $bim_fli --out $bfile_flipped --make-bed > /dev/null\n";
push @log_lines,  $cmd2;


&a2file($bim_log,@log_lines);




#print "create flipped dataset in this directory: $subdir\n";
&mysystem ($cmd1);
&mysystem ($cmd2);

&mysystem ("mv $bfile_flipped.bed $rootdir/");
&mysystem ("mv $bfile_flipped.bim $rootdir/");
&mysystem ("mv $bfile_flipped.fam $rootdir/");

&mysystem ("tar -cvzf $bfile_flipped.tar.gz $bim_inde $bim_xal $bim_gf $bim_bf $bim_fli $bim_uif $bim_usf $bim_nal");

&mysystem ("mv $bfile_flipped.tar.gz $rootdir/");


&mysystem ("mv $bim_log $rootdir/$bfile_flipped.report");

chdir ($rootdir);
&mysystem ("rm -r $workdir");


die $!."($bim_file).chefli.cmd" unless open BC, "> $bim_file.chefli.cmd";
foreach (@cmd_collect) {
    print BC "$_\n";
}
close BC;


exit;


&mysystem ("mkdir -p $subdir");
&mysystem ("cp -fl $bfile_flipped.bed $subdir/");
&mysystem ("cp -fl $bfile_flipped.bim $subdir/");
&mysystem ("cp -fl $bfile_flipped.fam $subdir/");

&mysystem ("cp -fl $bfile_flipped.bed ../");
&mysystem ("cp -fl $bfile_flipped.bim ../");
&mysystem ("cp -fl $bfile_flipped.fam ../");

#if ($replace){

#    &mysystem ("mv $bfile.bed $bfile.bed.orig");
#    &mysystem ("mv $bfile.bim $bfile.bim.orig");
#    &mysystem ("mv $bfile.fam $bfile.fam.orig");
#    &mysystem ("mv $bfile_flipped.bed $bfile.bed");
#    &mysystem ("mv $bfile_flipped.bim $bfile.bim");
#    &mysystem ("mv $bfile_flipped.fam $bfile.fam");


#}



print "checkflip completed, please read this report: $bim_log\n\n";



