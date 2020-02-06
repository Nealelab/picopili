#!/usr/bin/env perl
use strict;

#############################
# load utility functions
#############################



use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
use Cwd 'abs_path';
use rp_perl::Utils qw(trans);


my $version = "1.0.0";
my $progname = $0;
$progname =~ s!^.*/!!;

my $picodir = dirname(dirname(abs_path($0)));

#############################
# read config file
#############################

my $liloc = &trans("liloc");
# my $liref = "$picodir/lib/buigue";
my $liref = "/stanley/genetics/users/nbaya/spark/array_May2019/spark_ricopili/imputation/refdir_chr23"; #added for X chrom imputation


my $perlpack;
BEGIN {
	$perlpack = &trans("perlpack");
}
use lib $perlpack;

#####################################################
# use lib $ENV{rp_perlpackages};



my @bu_files;
# push @bu_files, "$liref/snp.txt.pos.scz49.gz";
# push @bu_files, "$liref/snp125.txt.pos.scz49.gz";
# push @bu_files, "$liref/snp130.txt.pos.scz49.gz";
# push @bu_files, "$liref/snp138.txt.pos.scz49.gz";
push @bu_files, "$liref/snp.txt.pos.scz49.chr23.gz";
push @bu_files, "$liref/snp125.txt.pos.scz49.chr23.gz";
push @bu_files, "$liref/snp130.txt.pos.scz49.chr23.gz";
push @bu_files, "$liref/snp138.txt.pos.scz49.chr23.gz";

my $redhat_version = `cat /etc/redhat-release`; # check redhat version
print "\nbuigue_pico_chr23.pl: $redhat_version";

my @li_files;
push @li_files, "$liloc/hg16ToHg19.over.chain.gz";
push @li_files, "$liloc/hg17ToHg19.over.chain.gz";
push @li_files, "$liloc/hg18ToHg19.over.chain.gz";
push @li_files, "$liloc/hg19ToHg19.over.chain.gz"; ## fake

my $lift_script = "lift_to_hg19.pl";


##### help message
my $usage = "
Usage : $progname bim-file

version: $version

  -help            print this message and exit

  --lift19         lift dataset to hg19



 guesses the build of a bim file out of ucsc snp file

  find here the helping files:
    $liref

 created by Stephan Ripke 2014 at MGH, Boston, MA
 in the frame of the PGC

";

use Getopt::Long;
GetOptions( 
    "help"=> \my $help,
    "lift19"=> \my $lift19,
    );


die ("wrong: $ARGV\n$usage") if ($help || @ARGV != 1);



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
    my $redhat_version = `cat /etc/redhat-release`; # check redhat version
    print "redhat_buigue: $redhat_version";

    system($systemstr);
    my $status = ($? >> 8);
    die "$systemstr\n->system call failed: $status" if ($status != 0);
}




my $bim_file = $ARGV[0];

unless ($bim_file =~ /.bim$/) {
    $bim_file .= ".bim";
}

my $bfile = $bim_file;
$bfile =~ s/.bim$//;

my %bpos=();

use Compress::Zlib ;

## read bim-file
print "read bim-file positons\n";
die "$!:$bim_file" unless open BIM, "< $bim_file";   
while (<BIM>) {
    my @cells = @{&split_line_ref(\$_)};
    my $pos = $cells[0]." ".$cells[3];
    $bpos{$cells[1]} = $pos;
}
close BIM;

my $nbim=keys(%bpos);
print "size of bim-file: $nbim\n\n";


my $max_comp = 50000000;
my $min_noma = 100000000;
my $buig;

print "start comparison\n";


die $!."($bim_file).buigue.noma_comp" unless open NC, "> $bim_file.buigue.noma_comp";

print NC "build";
print NC "\tN_comp";
print NC "\tN_nocomp";
print NC "\tN_match";
print NC "\tN_nomatch";
print NC "\n";
my $bucount = 0;
my $licount = 0;
foreach my $bufile (@bu_files) {
    ## compare with snp-collection

    my $bush = $bufile;
    $bush =~ s/.*snp/snp/;
    $bush =~ s/.txt.pos.gz//;
#    print "$bush\n";
#    exit;
    print "comparing with $bufile\n";
    
    my $nmatch= 0;
    my $nnomatch= 0;
    my $ncomp= 0;
    my $nnocomp= 0;

    my $sc = gzopen($bufile, "rb")  or die "Cannot open $bufile: $gzerrno\n" ;
    die $!."($bim_file).$bush.noma" unless open NOMA, "> $bim_file.$bush.noma";
    while ($sc->gzreadline(my $line) > 0) {

	my @cells = @{&split_line_ref(\$line)};


	if (exists $bpos{$cells[0]}) {
	    my $pos = $cells[1]." ".$cells[3];

#	    if ($cells[0] eq "rs3131972") {
#		print "rs3131972:\n";
#		print "$pos\n";
#		print "$bpos{$cells[0]}\n";
#	    }
	    if ($pos eq $bpos{$cells[0]}) {
		$nmatch++;
	    }
	    else {
		$nnomatch++;
		print NOMA "$cells[0] $bpos{$cells[0]} $pos\n";
	    }
	    $ncomp++;
	}
	else {
	    $nnocomp++;
	}
	last if ($ncomp == $max_comp);
    }	    
    close NOMA;
    $sc->gzclose();

    print NC "$bufile";
    print NC "\t$ncomp";
    print NC "\t$nnocomp";
    print NC "\t$nmatch";
    print NC "\t$nnomatch";
    print NC "\n";
    if ($nnomatch < $min_noma) {
	$min_noma = $nnomatch;
	$buig = $bush;
	$licount = $bucount;
    }

    $bucount++;
}
close NC;


my $cmd_out = "no liftover necessary to get to hg19\n";
if ($lift19){
    if ($licount < 3) {
	my $sys_str = "$lift_script --noex --lilofile $li_files[$licount] $bim_file";
	#	print "$sys_str\n";
	$cmd_out = $sys_str;
	&mysystem($sys_str);
    }
    else {
	&mysystem("ln -s $bfile.bed $bfile.hg19.bed");
	&mysystem("ln -s $bfile.bim $bfile.hg19.bim");
	&mysystem("ln -s $bfile.fam $bfile.hg19.fam");
    }
}

die $!."($bim_file).buigue.liftover_script" unless open BU, "> $bim_file.buigue.liftover_script";
print BU "$cmd_out\n";
close BU;



die $!."($bim_file).buigue" unless open BU, "> $bim_file.buigue";
print BU "$buig\n";
close BU;
print "success: $bim_file.buigue\n";
