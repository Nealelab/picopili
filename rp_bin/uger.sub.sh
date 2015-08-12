#!/usr/bin/env bash
#$ -j y
#$ -cwd
#$ -V

# wrapper script for job submission on Broad UGER cluster
#
# first parameter should be duration for 'sleep' before
# execution
# remainder of command line should be the job to be 
# submitted (including all agruments)
#
# The -V flag above will provoke a warning that 
# LD_LIBRARY_PATH won't be used for security reasons;
# this warning can be safely ignored

source /broad/software/scripts/useuse
reuse -q Anaconda
sleep $1
shift
"$@"
# eof
