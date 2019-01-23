#! /bin/bash
#SBATCH --time={wall_hours}:00:00
#SBATCH --nodes=1
#SBATCH --job-name {job_name}
#SBATCH -o {log_name}
#SBATCH -a 1-{array_jobs}

# sleep option (for preventing race conditions on network file systems)
sleep {sleep_time}

# setup resources
cd {workdir}
module load R

# main command line
{cmd_string}

# eof
