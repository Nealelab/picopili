#!/usr/bin/env perl
use strict;

my $version = "1.0.0";
my $progname = $0;


###################################################
## creates a subfile of reference based on bimfile
##
##
## created by Stephan Ripke 2012 at MGH, Boston, MA
##
## modified by Raymond Walters 2016
##
######################################################



##### help message
my $usage = "
Usage : $progname [options] bimfile

version: $version

  --ref      gzipped reference file
  --refheads SNP,CHR,POS,A1,A2,FRQ_A1
                 column headers for required info in ref file (if no chr, 
                 use NULL to default to --chr value)
  --chr INT  chromosome
  --help     print this help message and exit
 
";



use Getopt::Long;
GetOptions( 

    "help"=> \my $help,
    "ref=s"=> \my $reffile,
	"refheads=s"=> \my $refhead_str,
    "chr=i"=> \my $chrstr,
    );

die $usage if $help;


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



###################################################
###  system call with test if successfull
###################################################
sub mysystem(){
    my ($systemstr)="@_";
    system($systemstr);
    my $status = ($? >> 8);
    die "$systemstr\n->system call failed: $status" if ($status != 0);
}


use Compress::Zlib ;



###################################################
###  validate args
###################################################

my $bimfile = $ARGV[0];
unless (-e $bimfile) {
    print "ERROR: $bimfile is not existing\n";
    exit;
}
unless (-e $reffile) {
    print "ERROR: $reffile is not existing\n";
    exit;
}
unless ($chrstr) {
    print "ERROR: please specify chromosome --chr\n";
    exit;
}
my $outfile = $bimfile.".ref.chr$chrstr"; 

print "$bimfile\n$reffile\t$chrstr\n";

## check reffile has correct columns
# parse argument
my $snphead;
my $chrhead;
my $bphead;
my $a1head;
my $a2head;
my $frqhead;
($snphead, $chrhead, $bphead, $a1head, $a2head, $frqhead) = split(/,/, $refhead_str);
my @colnames = ($snphead, $chrhead, $bphead, $a1head, $a2head, $frqhead);

if( scalar @colnames != 6 ){
	die "Incorrect number of names for --refheads. Should be 6 (SNP,CHR,POS,A1,A2,FRQ_A1).\n";
}
if( "" ~~ @colnames ){
	die "Some elements of --refheads are empty.\n";
}

# read header in ref
my $inz = gzopen("$reffile", "rb")  or die "Cannot open $reffile\n" ;
$inz->gzreadline(my $head);
my @refheader = split /\s+/, $head;

# find columns matching requested names
my @snpcol_tmp = grep { @refheader[$_] =~ /^$snphead$/ } 0..$#refheader;
my @chrcol_tmp = grep { @refheader[$_] =~ /^$chrhead$/ } 0..$#refheader;
my @bpcol_tmp = grep { @refheader[$_] =~ /^$bphead$/ } 0..$#refheader;
my @a1col_tmp = grep { @refheader[$_] =~ /^$a1head$/ } 0..$#refheader;
my @a2col_tmp = grep { @refheader[$_] =~ /^$a2head$/ } 0..$#refheader;
my @frqcol_tmp = grep { @refheader[$_] =~ /^$frqhead$/ } 0..$#refheader;
my $snpcol;
my $chrcol;
my $bpcol;
my $a1col;
my $a2col;
my $frqcol;

if( scalar @snpcol_tmp > 1 ){
	die "SNP column name ($snphead) matches multiple columns of reference file header.\n"
}elsif(scalar @snpcol_tmp < 1 ){
	die "SNP column name ($snphead) not found in reference file header.\n";
}else{
	$snpcol = @snpcol_tmp[0];
}

if( scalar @bpcol_tmp > 1 ){
	die "BP position column name ($bphead) matches multiple columns of reference file header.\n"
}elsif(scalar @bpcol_tmp < 1 ){
	die "BP position column name ($bphead) not found in reference file header.\n";
}else{
	$bpcol = @bpcol_tmp[0];
}

if( scalar @a1col_tmp > 1 ){
	die "Allele1 column name ($a1head) matches multiple columns of reference file header.\n"
}elsif(scalar @a1col_tmp < 1 ){
	die "Allele1 column name ($a1head) not found in reference file header.\n";
}else{
	$a1col = @a1col_tmp[0];
}

if( scalar @a2col_tmp > 1 ){
	die "Allele2 column name ($a2head) matches multiple columns of reference file header.\n"
}elsif(scalar @a2col_tmp < 1 ){
	die "Allele2 column name ($a2head) not found in reference file header.\n";
}else{
	$a2col = @a2col_tmp[0];
}

if( scalar @frqcol_tmp > 1 ){
	die "A1 allele frequency column name ($frqhead) matches multiple columns of reference file header.\n"
}elsif(scalar @frqcol_tmp < 1 ){
	die "A1 allele frequency column name ($frqhead) not found in reference file header.\n";
}else{
	$frqcol = @frqcol_tmp[0];
}

if( scalar @chrcol_tmp > 1 ){
	die "Chromosome column name ($chrhead) matches multiple columns of reference file header.\n"
}elsif(scalar @chrcol_tmp < 1 ){
	# allow missing chr column
	if( $chrhead == "NULL" ){
		$chrcol = -9;
	}else{
		die "Chromosome column name ($chrhead) not found in reference file header.\n";
	}
}else{
	$chrcol = @chrcol_tmp[0];
}

# print "snp=$snpcol, chr=$chrcol, bp=$bpcol, a1=$a1col, a2=$a2col, frq=$frqcol\n";
# exit;


###############################
# BEGIN
#####################################

#read bimfile
print "read bimfile into hash\n";
my %snphash;
my %lochash;
my %lochash_frq;
my %loctrans;
die $!."$bimfile" unless open FILE, "< $bimfile";
die "$outfile.renames".$! unless open FD , "> $outfile.renames";
while (my $line = <FILE>){
    my @cells = @{&split_line_ref(\$line)};

    my $snp = $cells[1];
    my $chr = $cells[0] ;
    my $pos = $cells[3] ;



    my $chrpos = "$chr\t$pos";


    $snphash{$snp} = 1;
    if ($chrstr == $chr) {
	$lochash{$chrpos} = $snp;
    }



    if ($snp !~ /^rs[0-9]*$/) {


	if ($snp =~ /rs[0-9]/) {
	    my $csnp = $snp;
	    $csnp =~ s/.*(rs[0-9]*).*/\1/;
	    unless (exists $snphash {$csnp}) {
		if ($chr == $chrstr) {
		    print FD "rs_name extracted from snp_name: $snp $csnp\n";
		}
		$snphash{$csnp} = 1;
	    }
	}
	elsif  ($snp =~ /chr[0-9]*_[0-9]*/) {

	    my $psnp = $snp;
	    $psnp =~ s/.*(chr[0-9]*_[0-9]*).*/\1/;

	    my @tarr = split '_', $psnp;
	    my $cname = $tarr[0];
	    my $pname = $tarr[1];
	    $cname =~ s/chr//;
	    my $in = 0;
	    if ($chr == 0 && $pos == 0) {
		if ($chr == $chrstr) {
		    print FD "position extracted from snp_name: $snp $chr $pos $cname $pname\n";
		}
		$in =1;
	    }
	    else {
		if ($chr != $cname || $pos != $pname){
		    if ($chr == $chrstr) {
			print FD "Warning: pos_name doesnt match positions: $snp $cname $pname $chr $pos\n";
		    }
		    $in =1;
		}
	    }
	    if ($in ==1) {
		my $chrpos_loc = "$cname\t$pname";
		if ($cname==$chrstr) {
		    $lochash{$chrpos_loc} = $snp;
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
	    my $in = 0;
	    if ($chr == 0 && $pos == 0) {
		if ($chr == $chrstr) {
		    print FD "position extracted from snp_name: $snp $chr $pos $cname $pname\n";
		}
		$in =1;
	    }
	    else {
		if ($chr != $cname || $pos != $pname){
		    if ($chr == $chrstr) {
			print FD "Warning: pos_name doesnt match positions: $snp $cname $pname $chr $pos\n";
		    }
		    $in =1;
		}
	    }
	    if ($in ==1) {
		my $chrpos_loc = "$cname\t$pname";
		if ($cname==$chrstr) {
		    $lochash{$chrpos_loc} = $snp;
		}
	    }
	}

    }

}
close FILE;
close FD;




## read reference
# note: ref file already open from checking header
my $line;
die $!."$outfile.tmp" unless open OF, "> $outfile.tmp";
print OF "SNP\tCHR\tPOS\tA1\tA2\tFRQ_A1\n";
# expects: SNP	CHR	POS	A1	A2	FA1
while ($inz->gzreadline($line)){
    my @cells = @{&split_line_ref(\$line)};

    my $sname = $cells[$snpcol];
	my $rsid = $sname;
	if ($rsid !~ /^rs[0-9]*$/) {
		if ($rsid =~ /rs[0-9]/) {
			$rsid =~ s/.*(rs[0-9]*).*/\1/;
		}elsif ($rsid =~ /^[0-9]*:[0-9]*:.*:.*/) {
			$rsid =~ s/^([0-9]*:[0-9]*):.*:.*/\1/;
		}
	}	
	
	my $sbp = $cells[$bpcol];
	
	my $schr;
	if( $chrcol = -9 ){
		$schr = $chrstr;
	}else{
    	$schr = $cells[$chrcol];
	}
	my $pos = "$schr\t$sbp";
#	print "$sname\t$rsid\t$pos\n";

    my $in = 0;
    if (exists $snphash{$sname}) {
		$in = 1;
		delete $snphash{$sname};
#	delete $lochash{$pos};
    }elsif ( exists $snphash{$rsid} ){
    	$in = 1;
		delete $snphash{$rsid};
		$sname = $rsid;
    } 

    if  (exists $lochash{$pos}) {
		$in = 1;
		delete $lochash{$pos};
    }

    if ($in) {
#	$snpin{$sname} = 1;
		print OF "$sname\t$schr\t$sbp\t$cells[$a1col]\t$cells[$a2col]\t$cells[$frqcol]\n";
    }

}
$inz -> gzclose();
close OF;

die $!."$outfile.leftloc" unless open OF, "> $outfile.leftloc";
foreach my $locpos (keys %lochash) {
 #   unless (exists $snpin{$lochash{$locpos}}) {
	print OF "$locpos\t$lochash{$locpos}\n";
  #  }
}
close OF;

&mysystem ("mv $outfile.tmp $outfile");

print "done\n";

