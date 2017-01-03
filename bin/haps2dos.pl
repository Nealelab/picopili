#!/usr/bin/env perl
use strict;


####
# imp_prel.pl
# closely adapted from ricopili's haps2dos4
# Original by Stephan Ripke
# Adapted by Raymond Walters
#
####

my $version = "1.3.0";
my $progname = $0;
$progname =~ s!^.*/!!;

#############################
# load utility functions
#############################

use FindBin;
use lib "$FindBin::Bin";
use rp_perl::Utils qw(trans);


#############################
# read config file
#############################

my $ploc = &trans("p2loc");
my $perlpack = &trans("perlpack");
use lib $perlpack;
use Compress::Zlib ;

###############################################

my $plinkmem = 2000;

##### help message
my $usage = "
Usage : $progname haps-files (out of impute2)

version: $version

  --outname STRING    outdir, mandatory
  --outdir STRING     outname, mandatory
  --fam STRING        fam-file, mandatory
  --chr INT           chromosome
  --plinkmem INT      memory for plink --dosage (in MB), default $plinkmem
  --help              print this message and exit

 created by Stephan Ripke 2012 at MGH, Boston, MA
 in the frame of the PGC

";

use Getopt::Long;
GetOptions( 

    "help"=> \my $help,
    "outname=s"=> \my $outname,
    "outdir=s"=> \my $outdir,
    "fam=s"=> \my $famname,
	"idnum=s"=> \my $idnum,
    "chr=i"=> \my $chr,
    "plinkmem=i" => \$plinkmem,

    );

die ($usage) if $help;
die ($usage) unless $famname;
die ($usage) unless $outname;
die ($usage) unless $outdir;
die ($usage) unless $chr;


my $idnum_str = "";
if($idnum) {
	$idnum_str = "--update-parents $idnum";
}

# die "$usage" if (@ARGV != 1);


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

my @cmd_collect;

sub mysystem(){
    my ($systemstr)="@_";
    system($systemstr);
    my $status = ($? >> 8);
    die "$systemstr\n->system call failed: $status" if ($status != 0);
    push @cmd_collect, $systemstr;

}





############################################
######### BEGIN
########################################

my @haps_collection = @ARGV;
my $nf = @ARGV - 1;

my $haps_file1 = @haps_collection[0];
my $haps_file = "$haps_file1.combined";
my $haps_file_gz = "$haps_file1.combined.gz";




##########
# bring multiple files into one
########

my @filehandles_gz;


foreach my $infile (@haps_collection) {

    my $igz = gzopen("$infile", "rb")  or die "Cannot open file $infile: $gzerrno\n" ;
    push(@filehandles_gz, $igz);

}


print "open $haps_file_gz\n";

my $ogz = gzopen("$haps_file_gz", "wb")  or die "Cannot open file $haps_file_gz: $gzerrno\n" ;

my $igz = $filehandles_gz[0];


while ($igz->gzreadline(my $line)){
    chomp($line);
    my @cells = @{&split_line_ref(\$line)};

    if ($nf > 0) {
	foreach my $igzn (1..$nf) {
	    my $igz_add = $filehandles_gz[$igzn];
	    $igz_add->gzreadline(my $line);
	    chomp($line);
	    my @cells_add = @{&split_line_ref(\$line)};


	    foreach (1..5) {
		shift (@cells_add);
	    }

	    @cells = (@cells,@cells_add);
	    
	}
    }
    $ogz->gzwrite("@cells\n");

}

foreach my $igzn (1..$nf-1) {
    my $igz_add = $filehandles_gz[$igzn];
    $igz->gzclose();
}
$ogz->gzclose();

print "close $haps_file_gz\n";


##########
# process dosages
########

my $dosout = "$outdir/$outname";

# preprocess dosages to ensure probabilities sum to 1
my $twodos_tmp = "$haps_file1.combined.tmp2dos.gz";
my $twodos_tmp_success = "$twodos_tmp.fini";
my $imp_proc = "impprob_to_2dos $haps_file_gz $twodos_tmp";

unless (-e $twodos_tmp_success) {
    print "improb script: $imp_proc\n";
    &mysystem ($imp_proc);
}

# check success
# if successful, remove $haps_file_gz now to save space
if (-e $twodos_tmp_success) {
    &mysystem("rm $haps_file_gz")
} 
else {
    die "Failed to create $twodos_tmp from $haps_file_gz";
}

# convert to plink format
my $sys_loc = "$ploc/plink --threads 1 --memory $plinkmem --dosage $twodos_tmp noheader skip0=1 skip1=1 format=2  Zout --fam $famname $idnum_str --allow-no-sex --write-dosage --out $dosout";
print "$sys_loc\n";
&mysystem ($sys_loc);

# process corresponding info file
# my $mapout = "$outdir/$outname.out.dosage.map";
# my $info_file = $haps_file1;
# $info_file =~ s/.gz$//;
# $info_file .= "_info";
# if ($info_file) {
#    die $!."($info_file)" unless open IF, "< $info_file";
#    die $!."($mapout)" unless open MA, "> $mapout";
#    my $line = <IF>;
#    while (my $line = <IF>){
#	my @cells = @{&split_line_ref(\$line)};
#
#	my $snp = $cells[1];
#	my $pos = $cells[2];
#
#
#	my $bas_str = "$chr $snp 0 $pos";
#	print MA "$bas_str\n";
#    }
#    close IF;
#    close MA;
# }



my $outdosgz = "$dosout.out.dosage.gz";
my $finiout = "$outdir/$outname.out.dosage.fini";

if (-e $outdosgz) {
    &mysystem ("touch $finiout");
    &mysystem ("rm $twodos_tmp_success");
    &mysystem ("rm $twodos_tmp");
    &mysystem("rm $haps_file_gz");    
    print "done\n";
}else{
    die "Error: completed haps2dos.pl without creating $outdosgz. Something strange is wrong.\n";
}

