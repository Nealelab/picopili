#!/usr/bin/env perl
use strict;

my $version = "1.0.0";
my $progname = $0;


###################################################
## creates a subfile of reference based on bimfile
######################################################



##### help message
my $usage = "
Usage : $progname [options] bimfile

version: $version

  --help     print this help message and exit



 created by Stephan Ripke 2012 at MGH, Boston, MA
 
";



use Getopt::Long;
GetOptions( 

    "help"=> \my $help,

    );

die $usage if $help;

my $bimfile = $ARGV[0];
unless (-e $bimfile) {
    print "ERROR: $bimfile is not existing\n";
    exit;
}

my $outfile = $bimfile.".ref.sum"; 



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


###############################
# BEGIN
#####################################

#read bimfile
print "sum ref-files\n";

die "$outfile".$! unless open FD , "> $outfile";

# die $!."$bimfile.ref.chr1" unless open FILE, "< $bimfile.ref.chr1"; # old version for imputing across all autosomes
die $!."$bimfile.ref.chr23" unless open FILE, "< $bimfile.ref.chr23"; # for imputing chr 23
my $line = <FILE>;
print FD "$line";
close FILE;

my $chrloc = 23; # added for imputing chr 23
die $!."$bimfile.ref.chr$chrloc" unless open FILE, "< $bimfile.ref.chr$chrloc";
my $line = <FILE>;

while (my $line = <FILE>){
    print FD "$line";
}
close FILE;


close FD;


&mysystem ("tar -cvzf $bimfile.ref.addinfo.tar.gz $bimfile.ref.chr*.leftloc $bimfile.ref.chr*.renames");
#&mysystem ("gzip $outfile");


&mysystem ("touch $bimfile.ref.sum.done");
&mysystem ("rm $bimfile.ref.chr*");
print "done\n";

