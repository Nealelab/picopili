#!/usr/bin/env perl
use strict;
use File::Basename;
use Cwd;
use Cwd 'abs_path';
use Data::Dumper;

### Script to configure settings for picopili pipeline
### Jackie Goldstein, Jan 2014

### Adapted for picopili by Raymond Walters, Sept 2016

my $version = "2.0.0";
my $progname = $0;

$progname =~ s!^.*/!!;

my $cdir = abs_path($0);
my $home = $ENV{HOME};
my $conf_file = $ENV{HOME}."/picopili.conf";
my $command_line = "$progname @ARGV";

print "\n";
print "##############################\n";
print "#\n";
print "# Creating config file for picopili\n";
print "# $conf_file\n";
print "#\n";
print "# Will index location of executables for other\n";
print "# programs (e.g. plink), reference files, and\n";
print "# job settings (e.g. email address for job logs).\n";
print "#\n";
print "# Default settings are available for clusters\n";
print "# with existing ricopili configurations.\n";
print "#\n";
print "##############################\n";
    


#############################
# Ask user what cluster they're using
#############################
#my %clusters = ("broad",0,"mssm",0,"genomedk",0,"lisa",0,"other",0);
#my %clusters = ("broad",0,"mssm",0,"genomedk",0,"lisa",0,"computerome",0,"other",0);
my %clusters = ("broad",0,"mssm",0,"genomedk",0,"lisa",0,"computerome",0,"co_ipsych",0,"other",0);
my @cluster_names = ("broad","mssm","genomedk","lisa","computerome","co_ipsych","other");
print "Please enter your cluster name from the following options:\n";
my $i = 1;
foreach (@cluster_names){
    print "\t($i) $_\n";
    $i += 1;
}
print "\n";
my $cluster = "other";
while (1) {    
    $cluster = lc <>;
    chomp $cluster;
    if (exists $clusters{$cluster}){$clusters{$cluster} = 1;last;}
    else {
	$cluster =~ s/(\)|\()//g;
	if ($cluster >= 1 && $cluster <= $i){$cluster -= 1; $cluster = $cluster_names[$cluster];$clusters{$cluster} = 1;last;}
	else {
	    print "Did not recognize option. Please enter a cluster name from the options below:\n";
	    my $i = 1;
	    foreach (@cluster_names){
		print "\t($i) $_\n";
		$i += 1;
	    }
	    print "\n";
	    my $cluster = "other";
	}
    }
}
print "\nUsing the following cluster: $cluster\n\n";


###################################################
###  system call with test if successful
###################################################
sub mysystem(){
    my ($systemstr)="@_";
    system($systemstr);
    my $status = ($? >> 8);
    die "$systemstr\n->system call failed: $status" if ($status != 0);
}

###################################################
### Make sure lapack is installed
### specific to genomedk, and unclear if needed?
###################################################
# if ($clusters{genomedk} == 1){
#    unless ($ENV{EXTRAS} =~ /lapack/) {
#        print "Run the following commands to add lapack to the default search path:\n";
#        print "echo \"source /com/extra/lapack/3.5.0/load.sh\" >> ~/.bashrc\n";
#        print "source /com/extra/lapack/3.5.0/load.sh\n";
#        print "./rp_config\n\n";
#        exit;
#    }
#    else { print "Detected lapack is installed.\n\n";}
#}


###################################################
### Check whether to overwrite existing config (if exists)
###################################################

my $ans_ow = "y";
if (-e $conf_file) {
    print "Configuration file already exists at $conf_file\n";
    print "Do you wish to overwrite this file? (y/n)\n";
    while (1) {
        $ans_ow = lc <>;
        chomp $ans_ow;
        if ($ans_ow eq "y") {
            print "Rewriting configuration file. Making a backup to $conf_file.copy\n\n";
            &mysystem("cp $conf_file $conf_file.copy");
            last;
        }
        elsif ($ans_ow eq "n") {print "Not overwriting $conf_file.\n";last;}
        else {print "Please answer with y or n.\n";}
    }
};

my $cd = cwd();
my $sloc = "";
my $initials = "";
my $email = "";
my @text = ();

if ($ans_ow eq "y"){
#############################
# make scratch directory
#############################
if ($clusters{broad} == 1) {
    my $user_name = basename($ENV{HOME});
    $sloc = "/broad/hptmp/$user_name/";
    print "Do you want to use the following default scratch directory? (y or n)\n";
    print "\t$sloc\n";
    while (1) {
        my $answer = lc <>;
        chomp $answer;
        if ($answer eq "y") {print "Using $sloc as the scratch directory.\n\n";last;}
        elsif ($answer eq "n") {print "Please enter a scratch directory to use:\n";
				$sloc = <>;
				chomp $sloc;
				$sloc =~ s/^~/$ENV{HOME}/g;
				$sloc =~ s/^\./$cd/g;
				last;}
        else {print "Please answer with y or n.\n";}
    }
}
elsif ($clusters{lisa} == 1) {
    $sloc = "/scratch/";
    print "Do you want to use the following default scratch directory? (y or n)\n";
    print "\t$sloc\n";
    while (1) {
        my $answer = lc <>;
        chomp $answer;
        if ($answer eq "y") {print "Using $sloc as the scratch directory.\n\n";last;}
        elsif ($answer eq "n") {print "Please enter a scratch directory to use:\n";
				$sloc = <>;
				chomp $sloc;
				$sloc =~ s/^~/$ENV{HOME}/g;
				$sloc =~ s/^\./$cd/g;        
				last;}
        else {print "Please answer with y or n.\n";}
    }
}

elsif ($clusters{computerome} == 1) {
    $sloc = "/scratch/";
    print "Do you want to use the following default scratch directory? (y or n)\n";
    print "\t$sloc\n";
    while (1) {
        my $answer = lc <>;
        chomp $answer;
        if ($answer eq "y") {print "Using $sloc as the scratch directory.\n\n";last;}
        elsif ($answer eq "n") {print "Please enter a scratch directory to use:\n";
				$sloc = <>;
				chomp $sloc;
				$sloc =~ s/^~/$ENV{HOME}/g;
				$sloc =~ s/^\./$cd/g;        
				last;}
        else {print "Please answer with y or n.\n";}
    }
}

elsif ($clusters{co_ipsych} == 1) {
    $sloc = "/data/scratch/";
    print "Do you want to use the following default scratch directory? (y or n)\n";
    print "\t$sloc\n";
    while (1) {
        my $answer = lc <>;
        chomp $answer;
        if ($answer eq "y") {print "Using $sloc as the scratch directory.\n\n";last;}
        elsif ($answer eq "n") {print "Please enter a scratch directory to use:\n";
				$sloc = <>;
				chomp $sloc;
				$sloc =~ s/^~/$ENV{HOME}/g;
				$sloc =~ s/^\./$cd/g;        
				last;}
        else {print "Please answer with y or n.\n";}
    }
}

elsif ($clusters{genomedk} == 1) {
    $sloc = "/project/ricopili/scratch_dir/";
    print "Do you want to use the following default scratch directory? (y or n)\n";
    print "\t$sloc\n";
    while (1) {
        my $answer = lc <>;
        chomp $answer;
        if ($answer eq "y") {print "Using $sloc as the scratch directory.\n\n";last;}
        elsif ($answer eq "n") {print "Please enter a scratch directory to use:\n";
				$sloc = <>;
				chomp $sloc;
				$sloc =~ s/^~/$ENV{HOME}/g;
				$sloc =~ s/^\./$cd/g;        
				last;}
        else {print "Please answer with y or n.\n";}
    }
}


elsif ($clusters{mssm} == 1) {
    my $user_name = $ENV{USER};
    $sloc = "/sc/orga/scratch/$user_name/";
    print "Do you want to use the following default scratch directory? (y or n)\n";
    print "\t$sloc\n";
    while (1) {
        my $answer = lc <>;
        chomp $answer;
        if ($answer eq "y") {print "Using $sloc as the scratch directory.\n\n";last;}
        elsif ($answer eq "n") {print "Please enter a scratch directory to use:\n";
				$sloc = <>;
				chomp $sloc;
				$sloc =~ s/^~/$ENV{HOME}/g;
				$sloc =~ s/^\./$cd/g;        
				last;}
        else {print "Please answer with y or n.\n";}
    }
}
else {
    print "Please enter a scratch directory to use:\n";
    $sloc = "$cd/tmp/";
    while (1) {
        my $answer = lc <>;
        chomp $answer;
        if ($answer eq "y") {print "Using $sloc as the scratch directory.\n\n";last;}
        elsif ($answer eq "n") {print "Please enter a scratch directory to use:\n";
				$sloc = <>;
				chomp $sloc;
				$sloc =~ s/^~/$ENV{HOME}/g;
				$sloc =~ s/^\./$cd/g;        
				last;}
        else {print "Please answer with y or n.\n";}
    }

}

unless (-d $sloc) {
    print "Making scratch directory: $sloc\n\n";
    &mysystem("mkdir $sloc");    
}
else {
    print "Scratch directory already exists at $sloc\n";
}
print "\n";

#############################
# analyst info
#############################
print "Please enter your initials (2 letters):\n";
while (1) {
    $initials = lc <>;
    chomp $initials;
    if (length($initials) == 2) {last;}
    else {print "Make sure initials are 2 letters!\n";}
}
print "\n";

print "Please enter your email address (optional, can enter \"None\"):\n";
my $email = <>;
chomp $email;
print "\n";


#############################
# allow default all remaining values on select platforms
#############################
my $defall = 0;

if ($clusters{lisa} == 1 || $clusters{broad} == 1) {
    print "Do you want to use default values for the rest of the installation process? (y or n)\n";
    while (1) {
	my $answer = lc <>;
	chomp $answer;
	if ($answer eq "y") {print "Using default values for the rest of the installation process\n\n"; $defall = 1;last;}
	elsif ($answer eq "n") {print "Not using default values for the rest of the installation process\n\n"; $defall = 0;last;}
	else {print "Please answer with y or n.\n";}
    }
}



print "\n";

my %longvar = ("p2loc","PLINK2",
	       "shloc","SHAPEIT",
	       "i2loc","IMPUTE2",
	       "liloc","Liftover",
	       "eloc","Eigenstrat",
		   "admloc","ADMIXTURE",
		   "reaploc","REAP",
		   "priloc","PRIMUS",
	       "rloc","R",
		   "rscloc","Rscript",
		   "rplloc","R-enabled Plink (e.g. v1.07, or a dev build of 1.90)",
	       "hmloc","HapMap reference",
	       "perlpack","Perl packages (for Compress::Zlib)",
    );
	       

my %variables = ("p2loc", "",
		  "shloc","",
		 "i2loc","",
		  "liloc","",
		  "eloc","",
		  "admloc","",
		  "reaploc","",
		  "priloc","",
		  "rloc","",
		  "rscloc","",
		  "rplloc","",
		  "hmloc","",
		  "perlpack","",
    );

    
if ($clusters{broad}){
    %variables = (
		  "p2loc","/home/unix/sripke/plink_src/plink_1.9_newest",
		  "shloc","/home/unix/sripke/shapeit/bin",
		  "i2loc", "/psych/genetics_data/ripke/references_from_debakkerscratch/impute_v2/impute_v2/impute_2.2.7_beta",
		  "liloc","/home/unix/sripke/liftover",
		  "eloc","/humgen/atgu1/fs03/shared_resources/shared_software/EIG6.0beta_noreq/bin",
		  "admloc","/humgen/atgu1/fs03/shared_resources/shared_software/admixture_linux-1.23",
		  "reaploc","/humgen/atgu1/fs03/shared_resources/shared_software/REAP",
		  "priloc","/humgen/atgu1/fs03/shared_resources/shared_software/PRIMUS_v1.8.0/bin",
		  "rloc","/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.0.1/bin",
		  "rscloc","/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.0.1/bin",
		  "rplloc","/home/unix/sripke/plink_src/src/",
		  "hmloc","/psych/genetics_data/ripke/references_outdated/hapmap_ref",
		  "perlpack","/home/unix/sripke/perl_modules",
	);
}

		  # "ploc", "/home/unix/sripke/plink_src/src/",
#		  "ldloc","/psych/genetics_data/ripke/ldsc/",
#		  "rloc","broadinstitute",
#		  "rpac","NA",
#		  "meloc","/psych/genetics_data/ripke/references_from_debakkerscratch/metal/",
#		  "hvloc","/home/radon01/sripke/bakker_ripke/haploview/",


elsif ($clusters{lisa}){
    %variables = (
		  "p2loc","/home/gwas/plink2/plink_1.9_newest",
		  "shloc","/home/gwas/shapeit/bin",
		  "i2loc","/home/gwas/bin_impute_v2/impute_v2.2.2_x86_64_static",
		  "liloc","/home/gwas/liftover",
		  "eloc","/home/gwas/eigensoft",
		  "admloc","",
		  "reaploc","",
		  "priloc","",
		  "rloc","/sara/sw/R-3.1.2/bin",
		  "rscloc","/sara/sw/R-3.1.2/bin",
		  "rplloc","/home/gwas/plink/1.08/src",
		  "hmloc","/home/gwas/pgc-samples/hapmap_ref",
		  "perlpack","/home/gwas/perl_modules",
	);
}

#		  "ploc", "/home/gwas/plink/1.08/src",
#		  "ldloc","/home/gwas/ldsc/",
#		  "rloc","/sara/sw/R-3.1.2/bin/",
#		  "rpac","NA",
#		  "meloc","/home/gwas/metal",
#		  "hvloc","./",


elsif ($clusters{computerome}){
    %variables = (
		  "p2loc","/home/people/sripke/rp_external_bins/plink_1.9_newest",
		  "shloc","/home/people/sripke/rp_external_bins/shapeit/bin",
		  "i2loc","/home/people/sripke/rp_external_bins/impute2",
		  "liloc","/home/people/sripke/rp_external_bins/liftover",
		  "eloc","/home/people/sripke/rp_external_bins/EIG6.0beta",
		  "admloc","",
		  "reaploc","",
		  "priloc","",
		  "rloc","/services/tools/R-3.1.2/bin",
		  "rscloc","/services/tools/R-3.1.2/bin",
		  "rplloc","/services/tools/R-3.1.2/bin",
		  "hmloc","/home/people/sripke/imputation_references",
		  "perlpack","/home/people/sripke/rp_external_bins/perl_packages",
	);
}

#		  "ploc", "/home/people/sripke/rp_external_bins/plink/",
#		  "ldloc","/home/people/sripke/rp_external_bins/ldsc/",
#		  "rloc","/services/tools/R-3.1.2/bin/",
#		  "rpac","/home/people/sripke/rp_external_bins/Rpackages/",
#		  "meloc","/home/people/sripke/rp_external_bins/metal/",
#		  "hvloc","./",


elsif ($clusters{co_ipsych}){
    %variables = (
		  "p2loc","/data/tools/plink2_sept2015",
		  "shloc","/data/tools/shapeit_sept_2015/bin",
		  "i2loc","/data/tools/impute-2.3.2/",
		  "liloc","/data/user_tools/rp_external_bins/liftover",
		  "eloc","/data/tools/eigensoft-6.0.1/bin",
		  "admloc","",
		  "reaploc","",
		  "priloc","",
		  "rloc","/data/tools/R-3.2.1/bin",
		  "rscloc","/data/tools/R-3.2.1/bin",
		  "rplloc","/data/tools/plink-1.07",
		  "hmloc","/data/user_tools/imputation_references",
		  "perlpack","/data/user_tools/rp_external_bins/perl_packages",
	);
}

#		  "ploc", "/data/tools/plink-1.07/",
#		  "ldloc","/data/user_tools/rp_external_bins/ldsc/",
#		  "rloc","/data/tools/R-3.2.1/bin/",
#		  "rpac","/data/user_tools/rp_external_bins/Rpackages/",
#		  "meloc","/data/tools/metal-20110325/",
#		  "hvloc","./",


elsif ($clusters{genomedk}){
    %variables = (
		  "p2loc","/project/ricopili/plink_1.9_jul4",
		  "shloc","/project/ricopili/3rd_bins/shapeit/bin",
		  "i2loc","/project/ricopili/3rd_bins/impute2",
		  "liloc","/project/ricopili/3rd_bins/liftover",
		  "eloc","/project/ricopili/3rd_bins/eigenstrat/bin",
		  "admloc","",
		  "reaploc","",
		  "priloc","",
		  "rloc","/com/extra/R/3.1.0/bin",
		  "rscloc","/com/extra/R/3.1.0/bin",
		  "rplloc","/project/ricopili/plink_src",
		  "hmloc","/project/ricopili/reference_dir",
		  "perlpack","/project/ricopili/perl_packages/",
	);
}

#		  "ploc", "/project/ricopili/plink_src/",
#		  "rloc","/com/extra/R/3.1.0/bin",
#		  "rpac","NA",
#		  "meloc","/project/ricopili/3rd_bins/metal/",
#		  "hvloc","./",


elsif ($clusters{mssm}){
    %variables = (
		  "p2loc","/hpc/users/xripkes01/ricopili/3rd_binaries/plink-1.09-src-aug4",
		  "shloc","/hpc/users/xripkes01/ricopili/3rd_binaries/shapeit/bin",
		  "i2loc","/hpc/users/xripkes01/ricopili/3rd_binaries/impute2",
		  "liloc","/hpc/users/xripkes01/ricopili/3rd_binaries/liftover",
		  "eloc","/hpc/packages/minerva-common/eigensoft/5.0.1/bin",
		  "admloc","",
		  "reaploc","",
		  "priloc","",
		  "rloc","/hpc/packages/minerva-common/R/2.15.3/lib64/R/bin",
		  "rscloc","/hpc/packages/minerva-common/R/2.15.3/lib64/R/bin",
		  "rplloc","/hpc/users/xripkes01/ricopili/3rd_binaries/plink-1.07-src-sripke",
		  "hmloc","/hpc/users/xripkes01/ricopili/reference_dir",
		  "perlpack","/hpc/users/xripkes01/perl_modules",
	);
}

#		  "ploc", "/hpc/users/xripkes01/ricopili/3rd_binaries/plink-1.07-src-sripke/",
#		  "rloc","/hpc/packages/minerva-common/R/2.15.3/lib64/R/bin/",
#		  "rpac","NA",
#		  "meloc","/hpc/users/xripkes01/ricopili/3rd_binaries/metal/",
#		  "hvloc","./",


foreach (keys %variables){

	if ($variables{$_} ne '' && (-d $variables{$_})){
	    print "Default location for $longvar{$_} is: \n\t$variables{$_}\n\n";
	    if ($defall == 0) {
			print "Do you want to use this location (y or n)?\n";
			while (1) {
		    	my $answer = lc <>;
				chomp $answer;
				if ($answer eq "y") {
					print "Using $variables{$_} for $longvar{$_}.\n\n";
					last;
				}
				elsif ($answer eq "n") {print "Please enter a new location to use for $longvar{$_}:\n";
					my $input = <>;
					chomp $input;
					$input =~ s/^~/$ENV{HOME}/g;
					$input =~ s/^\./$cd/g;
					unless ( -d $input ){print "Not a valid directory. Please try again.\n";next;
				}
				print "\n";
				last;
			}
		    else {
				print "Please answer with y or n.\n";
			}
		}
	    }
	}
	else {
	    while (1){
			unless($clusters{other} == 1){
				print "No default value available for $longvar{$_}\n";
			}
			print "Please enter a location:\n";
			my $input = "";
			$input = <>;
			chomp $input;
			$input =~ s/^~/$ENV{HOME}/g;
			$input =~ s/^\./$cd/g;
			unless ( -d $input ){print "Not a valid directory. Please try again.\n";next;}
			$variables{$_} = $input;
			print "\n";
			last;
	    }
	}
}

foreach (keys %variables){
    push (@text, "$_ $variables{$_}");
}

push (@text, "sloc $sloc");
push (@text, "init $initials");
push (@text, "email $email");

### define queue depending on cluster

if ($clusters{broad}){push (@text, "queue broad_uger")}
if ($clusters{lisa}){push (@text, "queue lisa")}
if ($clusters{computerome}){push (@text, "queue computerome")}
if ($clusters{co_ipsych}){push (@text, "queue computerome_ipsych")}
if ($clusters{genomedk}){push (@text, "queue genomedk")}
if ($clusters{mssm}){push (@text, "queue mssm")}


} # end if block for getting conf file info


unless ( -e $conf_file && $ans_ow eq "n") {
    die $! unless open FILE, "> $conf_file";
    foreach (@text) {print FILE "$_\n"};
    close FILE;
}



#############################
# read picopili.conf file with default parameters
#############################
my %conf = (); ## hash with config parameters

### Read config file
die $!."($conf_file)" unless open FILE, "< $conf_file";
while (my $line = <FILE>){
    my @cells = split /\s+/, $line;
    $conf{$cells[0]} = $cells[1];
}
close FILE;

print "\n";

############################
# check whether all binary directories exist
############################
my @fail_path = ();
my %locs = ( 
		  "p2loc", "",
		  "shloc","",
		  "i2loc","",
		  "liloc","",
		  "eloc","",
		  "admloc","",
		  "reaploc","",
		  "priloc","",
		  "rloc","",
		  "rscloc","",
		  "rplloc","",
		  "hmloc","",
		  "perlpack","",
);

die $!."($conf_file)" unless open FILE, "< $conf_file";
while (my $line = <FILE>){
    my @cells = split /\s+/, $line;
    my $path = $cells[1];
    my $variable = $cells[0];
    unless (-d $path) {
        if (exists $locs{$variable}) {push(@fail_path,$variable)};
    }
}
close FILE;

#############################
# print finish statement
#############################

my $email_on = 0;
if ($conf{'email'} =~ m/\@/) {
	$email_on = 1;
}

if ($#fail_path != -1) { 

	print "You will need to install the binaries as described here (https://sites.google.com/a/broadinstitute.org/ricopili/resources) and use a text editor (emacs,vim,etc.) to edit the file paths listed in $home/picopili.conf for the following variables:\n";

	foreach (@fail_path) {
		print "\t$_\n";
	}
}
else {
	print "Configuration completed successfully! Settings are stored in $conf_file\n";
	
	if ($email_on){
		print "If you do not receive an email with the subject picopili_config, please check your address is entered correctly at $conf_file\n";
		&mysystem("echo \"Configuration for picopili was successful.\" | mail -s picopili_config  $conf{'email'}");
	}
}



###################################################
###  Optional: Add bin to default search path
###################################################

system("bin_check_pico"); # dummy script that doesn't do anything
my $status_bin = ($? >> 8);


if ($clusters{lisa} == 1) {
    unless (-e "$home/.bash_profile") {
	die $! unless open FILE, "> $home/.bash_profile";
	     print FILE 'if [ -f ~/.bashrc ]; then '."\n";
	     print FILE '    . ~/.bashrc'."\n";
	     print FILE 'fi'."\n";
	close FILE;
    }
    unless (-e "$home/.bashrc") {
	system "touch ~/.bashrc\n";
    }
}


if ($status_bin != 0) {
    my $bash = "PATH=\$PATH:$cdir";    
    my $csh = "set path=(\$path $cdir)";
    
    print "\n----------------------------------------------------\n";     
    print "## You will probably want to add picopili to the default search path.\n";
	
	
	# Determine the shell
	my $shell = '';
	if (exists $ENV{SHELL}){$shell = basename($ENV{SHELL});}
	if ($shell eq "bash-login-check"){$shell = "bash";}
	if ($shell ne "bash" && $shell ne "tcsh") {
    	print "Warning! Shell not recognized: $shell\n";
		print "Please send email to rwalters\@broadinstitute.org\n";
	}
	print "Detected you are using the following shell: $shell\n\n";
	
	# provide commands, where possible
	# perm tracks if a command is generated for .bashrc or equivalent
	my $perm = 0;
	if ($shell eq "bash"){
		print "To do this in bash, run the following command:\n";
		print "$bash\n";
		if ($clusters{broad}){
	        if (-e "$home/.my.bashrc") {
	            print "echo \"$bash\" >> ~/.my.bashrc\n";
				$perm = 1;
	        }
		}
		else{
			if ($clusters{lisa} == 1){
				unless (-e "$home/.bashrc") {
					print "touch ~/.bashrc\n";
				}
				unless (-e "$home/.bash_profile") {
					print "echo \"if [ -f ~/.bashrc ]; then \" > $home/.bash_profile";
					print "echo \"    ~/.bashrc\" >> $home/.bash_profile";
					print "echo \"fi\" >> $home/.bash_profile";
				}
			}
	        if (-e "$home/.bashrc") {
	            print "echo \"$bash\" >> ~/.bashrc\n";
				$perm = 1;
	        }
		}
	}
	elsif ($shell eq "tcsh"){
		print "To do this in tcsh, run the following command:\n";
		print "$csh\n";
		if ($clusters{broad}){
	        if (-e "$home/.my.cshrc") {
	            print "echo \"$csh\" >> ~/.my.cshrc\n";
				$perm = 1;
	        }
		}
		else{
	        if (-e "$home/.cshrc") {
	            print "echo \"$csh\" >> ~/.cshrc\n";
				$perm = 1;
	        }
		}
	}
	# else if shell not determined
    else {
        print "You'll want to add the following path:\n";
        print "\t$cdir\n";
	}
	# additional instructions of not .bashrc equivalent provided
	if ($perm == 0){
        print "If possible, add these paths permanently. Otherwise, you will need to do this everytime you start a new session.\n";		
        print "For example instructions, see http://www.cyberciti.biz/faq/unix-linux-adding-path/\n";		
	}
}
else{
	print "Successfully found picopili directory in search path!\n\n";	
}

print "### Finished ###\n\n";
    
exit;
########## Done ##########
