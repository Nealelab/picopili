#!/usr/bin/env perl
use strict;

#############################
# load utility functions
#############################

use File::Basename;
use FindBin;
use Cwd 'abs_path';
use lib "$FindBin::Bin";
use rp_perl::Utils qw(trans);

my $version = "1.0.0";
my $progname = $0;
$progname =~ s!^.*/!!;

my $picodir = dirname(dirname(abs_path($0)));

#############################
# read config file
#############################

my $hmloc = "$picodir/lib/plague";
my $perlpack = &trans("perlpack");
use lib $perlpack;

#####################################################

my $sc_file = "$hmloc/snp_platform_collection.txt.new.0815.gz";
my $sc_file_0416 = "$hmloc/snp_platform_collection.txt.new.0416a.gz";

my $scol = 2;

##### help message
my $usage = "
Usage : $progname bim-file

version: $version

  --scf    STRING  SNP collection file
                       default: $sc_file
		       first checking this: $sc_file_0416
  --scol INT       column of SNPs, default = $scol
  --create STRING  create new entry with name STRING
  -help            print this message and exit


 guesses the platform of a bim file out of experience

 created by Stephan Ripke 2008 at MGH, Boston, MA
 in the frame of the PGC

";

use Getopt::Long;
GetOptions( 
    "not"=> \my $notrans_temp,
    "scf=s" =>\ $sc_file,
    "scol=i" =>\ $scol,
    "create=s" =>\my $create,
    "bgl"=> \my $bgl_tmp,
    "phase2"=> \my $phase2,
    "help"=> \my $help,
    );


die ("wrong: $ARGV\n$usage") if ($help || @ARGV != 1);



##########################################
# subroutine to split a plink-output-line
##########################################

sub split_line {
    my ($line)=@_;
    chomp($line);
    $line =~ s/^[\s]+//g;
    my @cols=  split /\s+/, $line;
}



my $bim_file = $ARGV[0];

unless ($bim_file =~ /.bim$/) {
    $bim_file .= ".bim";
}

my %bsnps=();

use Compress::Zlib ;

## read bim-file

die "$!:$bim_file" unless open BIM, "< $bim_file";   
while (<BIM>) {
    my @cells = &split_line($_);
    $bsnps{$cells[$scol-1]} = 1;
}
close BIM;

my $nbim=keys(%bsnps);
print "size of bim-file: $nbim\n\n";

my @out_lines = ();

## compare with snp-collection

if (-e  $sc_file_0416) {
    $sc_file =  $sc_file_0416;
}

unless (-e $sc_file) {
    $sc_file = "$hmloc/snp_platform_collection.txt.new.0114.gz";
    if (-e $sc_file) {
	print "Warning: using older platform reference $sc_file\n";
	sleep(4);
    }
    else {
	print "Error: no platform reference found\n";
	exit;
    }
}
print "platform reference: $sc_file\n";


my $sc = gzopen($sc_file, "rb")  or die "Cannot open $sc_file: $gzerrno\n" ;
while ($sc->gzreadline(my $line) > 0) {
    chomp($line);
    push @out_lines, $line if ($create);
    my @cells=  split /\t/, $line;
    my @snps=  split ',', $cells[1];
    my $pos = 0;
    foreach (@snps){
	if (exists $bsnps{$_}){
	    $pos++;
	}
    }
    print "size of\t$cells[0] is\t".@snps.",\t$pos positives. ";
    printf "\t %.2f %% of this sample, these are \t%.2f %% of bim-file\n",$pos/@snps*100, $pos/$nbim*100;
}	    
$sc->gzclose();

if ($create) {
    my $out_str = "$create\t" ;
    foreach my $snp (keys %bsnps){
	$out_str .= ",$snp";
    }
    push @out_lines, $out_str;


    my $outz = gzopen("snp_platform_collection_new.txt.gz", "wb")  or die "Cannot open: $gzerrno\n" ;
    foreach (@out_lines){
	$outz->gzwrite("$_\n") ;
    }
    $outz->gzclose();

    print "if you want to keep:\n";
    print "cp snp_platform_collection_new.txt.gz $sc_file\n";
}



#cut -f2 mdd_stard_eur_QC1B.bim | tr "\n" "," | sed 's/^/\nstard_A500\t/' > snp_platform_collection.txt
#cut -f2 mdd_gain_eur_QC1B.bim  | tr "\n" "," | sed 's/^/\ngain_P600\t/' >> snp_platform_collection.txt
#cut -f2 mdd_genred_eur_QC1B.bim| tr "\n" "," | sed 's/^/\ngenred_A6.0\t/' >> snp_platform_collection.txt
#cut -f2 mdd_gsk_eur_QC1B.bim   | tr "\n" "," | sed 's/^/\ngsk_I550\t/' >> snp_platform_collection.txt
#cut -f2 mdd_munich_eur_QC1B.bim| tr "\n" "," | sed 's/^/\nmunich_I317\t/' >> snp_platform_collection.txt
