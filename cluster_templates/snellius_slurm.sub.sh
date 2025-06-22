#! /bin/bash
#SBATCH --time={wall_hours}:00:00
#SBATCH --mem={mem_in_gb}G
#SBATCH --nodes=1
#SBATCH --job-name {job_name}
#SBATCH -o {log_name}
#SBATCH --parsable
#SBATCH -a 1-{array_jobs}

# sleep option (for preventing race conditions on network file systems)
sleep {sleep_time}

# setup resources
cd {workdir}
source ~/.bashrc
echo $LD_LIBRARY_PATH
module load 2022
module load R/4.2.1-foss-2022a
module load Python/2.7.18-GCCcore-11.3.0-bare
module unload SciPy-bundle/2022.05-foss-2022a

# main command line
{cmd_string}

# eof
