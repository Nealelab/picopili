#!/usr/bin/env perl
use strict;



###########################################################################################################
#
#
#    checkpos6
#
#          Created by Stephan Ripke, Broadinstitute, sripke@broadinstitute.org
#
#                                  01/14/10
#
#          Adapted for Picopili by Raymond Walters, rwalters(at)broadinstitute.org
#
#
#
#
#    checks position bim-file (plink-binary-dataset) with dbsnp reference
#
#   checkpos3 looks also for positions and translates the SNP name if necessary
#   checkpos4 with plink2 support
#   checkpos5 with snp-translation if rs-name is buried, also positional translation
#   checkpos6 needs var_chr_renaming
##########################################################################################################

################################
## check multi-positions


#############################
# read config file
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

my $scol = 2;  ## snp-col in bim-file
my $ccol = 1;  ## chr-col in bim-file
my $kcol = 4;  ## kb-col in bim-file

my $dscol = 1;  ## snp-col in reference
my $dccol = 2;  ## chr-col in reference
my $dkcol = 3;  ## kb-col in reference

my $version = "1.0.0";
my $progname = $0;
$progname =~ s!^.*/!!;

my @new_arg = @ARGV;
pop (@new_arg);
my $cmd_line = "$progname @new_arg";

my $subdir = "dbsnp_subdir";

my $home_dir = "$ENV{HOME}";
my $dbsnp_file = "";


my $usage = "
Usage : $progname [options] bim-file

version: $version

  --dbsnp STRING      dbSNP reference file (created by readref)
  --ploc STRING       location of plink-binary (default read from picopili.conf)
                       default: $p2loc
  --col INT,INT,INT   snp-col,chr-col,kb-col in bim-file: default: $scol,$ccol,$kcol 
  --dbcol INT,INT,INT snp-col,chr-col,kb-col in dbsnp-file: default: $dscol,$dccol,$dkcol 

  --nocreate          do not create clean datasets, only analyse
  --exmulti           exclude all multiple annotated
  --nokeep            keepnon-annoted
  --subdir STRING     subdir, to put end-dataset into, default: $subdir
  --help              print this message and exit

  --exmulti, --nokeep and --subdir are in effect, as long as --ncreate is not switched

";

use Getopt::Long;
GetOptions( 
    "dbsnp=s"=> \$dbsnp_file,
    "ploc=s"=> \$p2loc,
    "col=s"=> \my $colstr,
    "dbcol=s"=> \my $dcolstr,
    "help"=> \my $help,
    "nocreate"=> \my $nocreate,
    "nokeep"=> \my $nokeep,
    "exmulti"=> \my $exmulti,
    "subdir=s"=> \$subdir,
    );

die "$usage\n" if (@ARGV != 1 || $help);

if ($colstr) {
    ($scol,$ccol,$kcol) = split(/,/, $colstr);
    die "*****\nwrong col-usage\n****\n$usage\n" if ($kcol eq "");
}
if ($dcolstr) {
    ($dscol,$dccol,$dkcol) = split(/,/, $dcolstr);
    die "*****\nwrong col-usage\n****\n$usage\n" if ($dkcol eq "");
}

unless (-e $dbsnp_file) {
    die "*** Error, dbSNP file not found\n";
}


#############################
# test, if plink is present
#############################


unless (-e "$p2loc/plink" ){
    print "\n***** Error: couldn't find the following:\n";
    print "please check --p2loc or $conf_file\n";
    exit;
}



my $bim_file=$ARGV[0];
my $bim_sorted=$bim_file.".sorted";
my $bim_joined=$bim_file.".joined";
my $bim_addpos=$bim_file.".addpos";
my $bim_xpos=$bim_file.".xpos";
my $bim_xchr=$bim_file.".xchr";
my $bim_xkb=$bim_file.".xkb";
my $bim_xdup=$bim_file.".xdup";
my $bim_npos=$bim_file.".npos";
my $bim_rpos=$bim_file.".rpos";
my $bim_uchr=$bim_file.".uchr";
my $bim_ukb=$bim_file.".ukb";
my $bim_extr=$bim_file.".extract";



my $bim_updated=$bim_file.".updated";

my $bfile = $bim_file;
print "Warning, no bim_file no create please\n" unless ($bfile =~ /.bim$/);
$bfile =~ s/.bim$//;

my $bfile_chr = $bfile."_chr";
my $bfile_kb = $bfile."_kb";
my $bfile_dbsnp = $bfile.".ch";

my $bim_log=$bfile_dbsnp.".report";

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
    die "$!: $file" unless open FILE, "> $file";
    foreach (@lines){
	print FILE $_;
    }
    close FILE;
}

###############################
###   BEGIN
######################


use File::Copy;
use File::Path;
use Cwd;


my $rootdir = &Cwd::cwd();
my $work_dir = "$sloc/checkpos_$bim_file";
while (-e $work_dir) {
    $work_dir .= ".p";
}
print "work_dir: $work_dir\n";
#&mysystem ("rm -rf $work_dir");

my @created = mkpath(
    $work_dir,
    {verbose => 0, mode => 0750},
    );

chdir $work_dir or die $!;

&mysystem("ln -sf $rootdir/$bfile.bim .") unless (-e "$bfile.bim");
&mysystem("ln -sf $rootdir/$bfile.bed .") unless (-e "$bfile.bed");
&mysystem("ln -sf $rootdir/$bfile.fam .") unless (-e "$bfile.fam");


my @log_lines;
push @log_lines,  "\nsummary checkpos, version: $version\n\n";



print "read bim-file\n";


################################################
### read all snp names
########################################################

my %bim_hash;
die "$bfile.bim".$! unless open FI , "< $bfile.bim";
while (my $line = <FI>){

    my @cells = @{&split_line_ref(\$line)};
    my $snp = $cells[$scol-1];
    if (exists $bim_hash{$snp}) {
	print "double snp name entry: $snp\n";
	exit;
    }
    $bim_hash{$snp} = 1;

}
close FI;




################################################
### work on snp names
########################################################

my $name_rs = 0;
my $name_nrs = 0;
my $name_pos = 0;
my $name_war = 0;
die "$bfile.bim".$! unless open FI , "< $bfile.bim";
die "$bfile.bim.ow".$! unless open FO , "> $bfile.bim.ow";
die "$bfile.bim.ow.det".$! unless open FD , "> $bfile.bim.ow.det";
while (my $line = <FI>){

    my @cells = @{&split_line_ref(\$line)};
    my $snp = $cells[$scol-1];
    my $chr = $cells[$ccol-1] ;
    my $pos = $cells[$kcol-1] ;



    if ($snp !~ /^rs[0-9]*$/) {
#	print "$snp\t";

	if ($snp =~ /rs[0-9]/) {
	    my $csnp = $snp;
	    $csnp =~ s/.*(rs[0-9]*).*/\1/;
	    unless (exists $bim_hash {$csnp}) {
		$cells[$scol-1] = $csnp;
		print FD "rs_name extracted from snp_name: $snp $csnp\n";
		$name_rs++;
		$bim_hash {$csnp} = 1;
	    }
	    else {
		print FD "rs_name not extracted from snp_name, since already present: $snp $csnp\n";
		$name_nrs++;
	    }
	}
	elsif  ($snp =~ /chr[0-9]*_[0-9]*/) {

	    my $psnp = $snp;
	    $psnp =~ s/.*(chr[0-9]*_[0-9]*).*/\1/;

	    my @tarr = split '_', $psnp;
	    my $cname = $tarr[0];
	    my $pname = $tarr[1];
	    $cname =~ s/chr//;
	    if ($chr == 0 && $pos == 0) {
		$cells[$ccol-1] = $cname;
		$cells[$kcol-1] = $pname;
		print FD "position extracted from snp_name: $snp $chr $pos $cname $pname\n";
		$name_pos++;
	    }
	    else {
		if ($cells[$ccol-1] != $cname || $cells[$kcol-1] != $pname){
		    print FD "Warning: pos_name doesnt match positions: $snp $cells[$ccol-1] $cells[$kcol-1] $cname $pname\n";
		    $name_war++;	    
		}
	    }
	}


	elsif  ($snp =~ /var_[0-9]*_[0-9]*/) {

	    my $psnp = $snp;
	    $psnp =~ s/.*(var_[0-9]*_[0-9]*).*/\1/;
	    $psnp =~ s/var_//;

	    my @tarr = split '_', $psnp;
	    my $cname = $tarr[0];
	    my $pname = $tarr[1];
	    $cname =~ s/chr//;
#	    my $in = 0;
	    if ($chr == 0 && $pos == 0) {
#		if ($chr == $chrstr) {
		    $cells[$ccol-1] = $cname;
		    $cells[$kcol-1] = $pname;
		    print FD "position extracted from snp_name: $snp $chr $pos $cname $pname\n";
#		}
#		$in =1;
	    }
	    else {
		if ($chr != $cname || $pos != $pname){

			print FD "Warning: pos_name doesnt match positions (still translated): $snp $chr $pos $cname $pname\n";
			$cells[$ccol-1] = $cname;
			$cells[$kcol-1] = $pname;

		}
	    }
#	    if ($in ==1) {
#		my $chrpos_loc = "$cname\t$pname";
#		if ($cname==$chrstr) {
#		    $lochash{$chrpos_loc} = $snp;
#		}
#	    }
	}




    }
    print FO "@cells\n";
}
close FI;
close FO;
close FD;


print "number of rs_name extractions: $name_rs\n";
print "number of no rs_name extractions: $name_nrs\n";
print "number of position extractions: $name_pos\n";
print "number of warning extractions: $name_war\n";

#exit;

################################################
### read bim file into hash
########################################################

#my @bim_lines;
my %bim_info; # position for snp_name
my %loc_info; # also read all positions, safes SNP name for position in bim-file

die "$bfile.bim.ow".$! unless open FILE , "< $bfile.bim.ow";
while (my $line = <FILE>){
    my @cells = @{&split_line_ref(\$line)};
    my $snp = $cells[$scol-1];
    $bim_info{$snp} = " $cells[$ccol-1] $cells[$kcol-1]";


#    if ($snp eq "rs493347") {
#	print "bim-found: rs493347\n";
#	sleep(10);
 #   }

    if ($cells[$ccol-1] != 0 || $cells[$kcol-1] != 0) {
	my $posloc = " $cells[$ccol-1] $cells[$kcol-1]";
#	if (exists $loc_info{$posloc}){
#	    print "double entry on position: $posloc\n";
#	    exit;
#	}
	$loc_info{$posloc} = $cells[$scol-1];
    }
    else {
	## check the presumably right positions
	if ($cells[$kcol-1] !=0) {
	    print  "chr position equals zero, but basepair ($cells[$kcol-1]) not: $snp\n";
	    exit;
	}
	if ($cells[$ccol-1] !=0) {
	    print  "basepair position equals zero, but chr position ($cells[$ccol-1]) not: $snp\n";
	    exit;
	}
    }

#    push @bim_lines,$line;
}
close FILE;

#if (exists $bim_info{"rs493347"}) {
#    print "rs493347 exists in bim file\n";
#}


######################################################################
## read reference file
#####################################################################
print "read reference-file: $dbsnp_file\n";
my %bim_in_ref;  ## safe snp names that are found in reference
die "$dbsnp_file".$! unless open FILE , "< $dbsnp_file";
die "$bim_joined".$! unless open OF , "> $bim_joined";
die "$bim_addpos".$! unless open AP , "> $bim_addpos";
while (my $line = <FILE>){
    my @cells = @{&split_line_ref(\$line)};
    my $info_loc = " $cells[$dccol-1] $cells[$dkcol-1]";


#    if ($cells[$dscol-1] eq "rs493347") {
#	print "found in ref: rs493347\n";
#	sleep(10);
#    }

    ###########################################
    ### this is crucual:
    ###  only if the SNP is not present on the bim file, there is a position lookup
    ##############################################

    if (exists $bim_info{$cells[$dscol-1]}) {

#	if ($cells[$dscol-1] eq "rs493347") {
#	    print "found still in bim: rs493347\n";
#	    sleep(10);
#	}

	$bim_in_ref{$cells[$dscol-1]} = 1;

	print OF $cells[$dscol-1].$bim_info{$cells[$dscol-1]}.$info_loc."\n";
	delete ($bim_info{$cells[$dscol-1]});
    }
    elsif (exists $loc_info{$info_loc}) {

	## loc_info: snp_name for position in bim fie
	## $cells[$dscol-1]: snp name in reference
	## $info_loc: position
	## $bim_info{$cells[$dscol-1]}: position of reference SNP in bim-file

#	if ($cells[$dscol-1] eq "chr10_116123724_D") {
#	    print $cells[$dscol-1]."\n";
#	    print $loc_info{$info_loc}."\n";
#	    print $bim_info{$cells[$dscol-1]}."\n";
#	    print $info_loc."\n";
#	    print "all\n";
#	}

#if ()

#	$bim_in_ref{$cells[$dscol-1]} = 1;
	
	print AP $cells[$dscol-1]." ".$loc_info{$info_loc}." ".$info_loc."\n";
#	delete ($bim_info{$loc_info{$info_loc}});
    }
}
close FILE;

print "print remainder\n";
foreach my $snp_loc (keys %bim_info) {
    print OF $snp_loc.$bim_info{$snp_loc}."\n";
}
close OF;
close AP;

#exit;





### read bim-file

print "read bim-file\n";
my %mupos; ## multiple positions in reference file
my %pos_count;
my %multi_ref;




###################################
### compare to position_reference
##################################
my @snp_xpos;  ## snps with bad position (chr and kb)
my @snp_xchr;  ## snps with bad position (chr)
my @snp_xkb;  ## snps with bad position (kb)
my @snp_kb_update;  ## snps with updateinformation
my @snp_chr_update;  ## snps with updateinformation
my @snp_name_update;  ## snps with updateinformation
my @snp_xdup;  ## snps with bad position (chr and kb)
my @snp_rpos; ## snps with good position (sanity check)
my @non_pos;  ## snps without reference position

my %ukb;  ## update position
my %uchr;  ## update chromosome
my %usnp;  ## update snpname
my $usnp;  ## count
my $unsnp = 0;  ## non update snpname count

print "analyse joined file\n";
die $!." <$bim_joined>" unless open POS, "< $bim_joined";
my $cc = 0;
my $old_snp = "";
while (my $line = <POS>) {
    $cc++;
    print "$cc entries analysed\n" if ($cc % 100000 == 0);
    my @cells = @{&split_line_ref(\$line)};


    my $snp = $cells[0];

    my $crefcol = 3;
    my $krefcol = 4;
    my $chr_ref = "$cells[$crefcol]";
    my $kb_ref = "$cells[$krefcol]";
    my $chr_bim = "$cells[1]";
    my $kb_bim = "$cells[2]";

#   print "$line\n";
#    print "$snp\t$chr_bim\t$kb_bim\t$chr_ref\t$kb_ref\n";
#    exit;
	
    if ($chr_ref ne "") {

	if (($chr_bim eq $chr_ref) && ($kb_bim eq $kb_ref)){
	    push @snp_rpos, $snp."\t"."$chr_ref:$kb_ref"."\n";
	}
	else {
#	    if (($chr_bim ne $chr_ref) && ($kb_bim ne $kb_ref)){
#		push @snp_xpos, $snp."\t"."$chr_bim:$kb_bim"."\t"."$chr_ref:$kb_ref"."\n";
#	    }
	    if ($chr_bim ne $chr_ref){
		push @snp_xchr, $snp."\t".$chr_bim."\t".$chr_ref."\n";
		push @snp_chr_update, $snp."\t".$chr_ref."\n";
		$uchr {$snp} = $chr_ref;
	    }
	    if ($kb_bim ne $kb_ref){
		my $dist = $kb_ref - $kb_bim;
		push @snp_xkb, $snp."\t".$chr_bim."\t".$kb_bim."\t".$kb_ref."\t".$dist."\n";
		push @snp_kb_update, $snp."\t".$kb_ref."\n";
		$ukb {$snp} = $kb_ref;
	    }
	}

	$pos_count{"$snp\t$chr_bim\t$kb_bim"}++ ;  ## count occurences
	$multi_ref{"$snp\t$chr_bim\t$kb_bim"} .= "\t$chr_ref:$kb_ref";  ## count occurences
	$mupos{$snp} = 1 if ($pos_count{"$snp\t$chr_bim\t$kb_bim"} > 1);
    }
    else {
	my $info_loc = " $chr_bim $kb_bim";
	unless (exists $loc_info{$info_loc}) {
	    push @non_pos, $snp."\t$chr_bim:$kb_bim\n";
	}
    }
	
#    delete ($pos_chr{$snp});
}
close POS;


## print out multiple-annotated SNPs
foreach (keys %pos_count) {
    if ($pos_count{$_} > 1){
	push @snp_xdup, $_."$multi_ref{$_}\n";
    }
}



print "read addpos file for translating SNP names\n";
die $!." <$bim_addpos>" unless open AP, "< $bim_addpos";
die $!." <$bim_addpos.det>" unless open AD, "> $bim_addpos.det";
while (my $line = <AP>) {
    $cc++;
    print "$cc entries analysed\n" if ($cc % 100000 == 0);
    my @cells = @{&split_line_ref(\$line)};
    unless (exists $bim_in_ref{$cells[1]}) {
	print AD "$cells[1] $cells[0] $cells[2] $cells[3] replace\n";
	$usnp {$cells[1]} = $cells[0];
	$usnp++;
    }
    else {
	print AD "$cells[1] $cells[0] $cells[2] $cells[3] non_replace\n";
	$unsnp++;
    }
}
close AP;
close AD;

##########################################
# exclude multiple SNPs from update lists.
##########################################

sub excl_multi {
    my @list=@_;
    my @tmp;
    foreach (@list) {
	my @c = split "\t", $_;
	push @tmp, $_ unless (exists $mupos{$c[0]})
    }
    @tmp;
}

@snp_chr_update = &excl_multi (@snp_chr_update);
@snp_kb_update = &excl_multi (@snp_kb_update);
@snp_xchr = &excl_multi (@snp_xchr);
@snp_xkb = &excl_multi (@snp_xkb);
@snp_xpos = &excl_multi (@snp_xpos);

print "write files...\n";
&a2file($bim_xpos,@snp_xpos);
&a2file($bim_xchr,@snp_xchr);
&a2file($bim_xkb,@snp_xkb);


&a2file($bim_rpos,@snp_rpos);
&a2file($bim_npos,@non_pos);
&a2file($bim_xdup,@snp_xdup);
&a2file($bim_ukb,@snp_kb_update);
&a2file($bim_uchr,@snp_chr_update);


push @log_lines,  "\nsnp-name extractions:\n";
push @log_lines,  "rs_name extraction           :\t".$name_rs."\t-> $bfile.bim.ow.det\n";
push @log_lines,  "rs_name extraction (not done):\t".$name_nrs."\t-> $bfile.bim.ow.det\n";
push @log_lines,  "position extraction          :\t".$name_pos."\t-> $bfile.bim.ow.det\n";
push @log_lines,  "warning extraction           :\t".$name_war."\t-> $bfile.bim.ow.det\n";
push @log_lines,  "\nsnp-names taken from reference based on position:\n";
push @log_lines,  "as noted above               :\t".$usnp."\t-> $bim_addpos.det\n";
push @log_lines,  "not done since snpname exists:\t".$unsnp."\t-> $bim_addpos.det\n";
push @log_lines,  "\npositions of rs-name in db-snp:\n";
push @log_lines,  "wrong (chr and kb)           :\t".@snp_xpos."\t-> $bim_xpos\n";
push @log_lines,  "wrong (chr), update possible :\t".@snp_xchr."\t-> $bim_xchr\n";
push @log_lines,  "wrong (kb), update possible  :\t".@snp_xkb."\t-> $bim_xkb\n";
push @log_lines,  "not found                    :\t".@non_pos."\t-> $bim_npos\n";
push @log_lines,  "multiple annot in dbSNP      :\t".@snp_xdup."\t-> $bim_xdup\n";
push @log_lines,  "\nsanity check               :\n";
push @log_lines,  "right position               :\t".@snp_rpos."\t-> $bim_rpos\n";
push @log_lines,  "N - SNPs in bim-file\n\n";
push @log_lines,  "** multiple Annotation contains possibly SNPs from the wrong-section\n\n";



push @log_lines,  "\n\n -> have a look at the ouput files\n";


print "update bim...\n";
#foreach my $line (@bim_lines) {
die "$bfile.bim.ow".$! unless open FILE , "< $bfile.bim.ow";
die "$bim_updated".$! unless open OF , "> $bim_updated";

while (my $line = <FILE>){
    chomp ($line);
    my @cells = @{&split_line_ref(\$line)};
    if (exists $usnp {$cells[1]}) {
	$cells[1] = $usnp {$cells[1]};
    }
    if (exists $uchr {$cells[1]}) {
	$cells[0] = $uchr {$cells[1]};
    }
    if (exists $ukb {$cells[1]}) {
	$cells[3] = $ukb {$cells[1]};
    }
    print OF "$cells[0]";
    print OF "\t$cells[1]";
    print OF "\t$cells[2]";
    print OF "\t$cells[3]";
    print OF "\t$cells[4]";
    print OF "\t$cells[5]\n";
}

close FILE;
close OF;

print "copy files...\n";
#&mysystem ("mv $bim_updated $bim_file");
&mysystem ("cat $bim_rpos $bim_ukb $bim_uchr $bim_addpos > $bim_extr");


my $exclude_txt = ""; 
$exclude_txt = "--exclude $bim_xdup" if ($exmulti);
my $cmd3 = "$p2loc/plink --memory 2000 --bed $bfile.bed --fam $bfile.fam --bim $bim_updated --out $bfile_dbsnp --extract $bim_extr $exclude_txt --make-bed \n";

print "make new bed...\n";
&mysystem ($cmd3) unless (-e "$bfile_dbsnp.bed");


push @log_lines, $cmd3;
&a2file($bim_log,@log_lines);

&mysystem ("mv $bfile_dbsnp.bed $rootdir/");
&mysystem ("mv $bfile_dbsnp.bim $rootdir/");
&mysystem ("mv $bfile_dbsnp.fam $rootdir/");


&mysystem ("tar -cvzf $bfile_dbsnp.tar.gz $bim_xpos $bim_xchr $bim_xkb $bim_rpos $bim_npos $bim_xdup $bim_ukb $bim_uchr $bim_addpos.det $bfile.bim.ow.det");

&mysystem ("mv $bfile_dbsnp.tar.gz $rootdir/");
&mysystem ("mv $bim_log $rootdir/");
&mysystem ("mv $bim_npos $rootdir/");
&mysystem ("mv $bim_ukb $rootdir/");
&mysystem ("mv $bfile.bim.ow.det $rootdir/");

die $!."($bim_file).chepos.cmd" unless open BC, "> $bim_file.chepos.cmd";
foreach (@cmd_collect) {
    print BC "$_\n";
}
close BC;

&mysystem ("mv  $bim_file.chepos.cmd $rootdir/");
chdir ($rootdir);
sleep(1);
#exit;
&mysystem ("rm -r $work_dir");


exit;


my $cmd1 = "$p2loc/ploc --memory 2000 --bfile $bfile --update-map $bim_uchr --update-chr --out $bfile_chr --make-bed > /dev/null\n";
my $cmd2 = "$p2loc/ploc --memory 2000 --bfile $bfile --update-map $bim_ukb --out $bfile_kb --make-bed > /dev/null\n";

my $nokeep_txt = ""; 
$nokeep_txt = "cat $bim_npos >> $bim_rpos; " if ($nokeep);
$cmd3 = $nokeep_txt.$cmd3;





unless ($nocreate){
    if (@snp_xchr > 0) {
	push @log_lines, $cmd1;
	&a2file($bim_log,@log_lines);
	print "update chr-information, have a look here: $bim_log\n";
	exit;
	&mysystem ($cmd1);
	&mysystem ("$cmd_line $bfile_chr.bim");
	exit;
    }
}

unless ($nocreate){
    if (@snp_xkb > 0) {
	push @log_lines, $cmd2;
	&a2file($bim_log,@log_lines);
	print "update kb-information, have a look here: $bim_log\n";
	&mysystem ($cmd2);
	&mysystem ("$cmd_line $bfile_kb.bim");
	exit;
    }
}

push @log_lines, $cmd3;
&a2file($bim_log,@log_lines);

unless ($nocreate){
    print "create dbsnp dataset, have a look here: $bim_log\n";
    &mysystem ($cmd3);
    &mysystem ("mv $bfile_dbsnp.bed $rootdir/");
    &mysystem ("mv $bfile_dbsnp.bim $rootdir/");
    &mysystem ("mv $bfile_dbsnp.fam $rootdir/");

#    &mysystem ("$cmd_line --nocreate $bfile_dbsnp.bim");
    exit;
}

chdir ($rootdir);


&mysystem ("rm -r $work_dir");

print "have a look at this file: $bim_log\n";






