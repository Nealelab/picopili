#!/usr/bin/env bash

# wrapper script for job submission on Broad UGER cluster
#
# The -V below above will provoke a warning that
# LD_LIBRARY_PATH won't be used for security reasons;
# this warning can be safely ignored

#$ -j y
#$ -cwd
#$ -V
#$ -N {job_name}
#$ -o {log_name}
#$ -q {queue_name}
#$ -l m_mem_free={mem_in_gb}g,h_vmem={mem_in_gb}g
::PICO_ARRAY_ONLY::#$ -t 1-{array_jobs}
::PICO_ARRAY_ONLY::#$ -tc {array_max}
::PICO_THREAD_ONLY::#$ -pe smp {threads}

# sleep option (for preventing race conditions on network file systems)
sleep {sleep_time}

# setup resources
source /broad/software/scripts/useuse
reuse -q Anaconda
reuse -q .curl-7.47.1
reuse -q .cairo-1.14.2

# main command line
{cmd_string}

# eof
