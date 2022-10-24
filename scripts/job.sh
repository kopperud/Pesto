#!/usr/bin/env sh
#SBATCH --job-name=bears_SCM
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem=2GB
#SBATCH --output=logs/bears_BDS-%a.log
#SBATCH --error=logs/bears_BDS-%a.err
#SBATCH --qos=normal_prio
#SBATCH --ntasks=2

module load gnu openmpi

#source "scripts/env.sh"

TIMESLICES=$(cat scripts/arg_list.txt | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

mpirun -np 2 rb-mpi --args ${TIMESLICES} --file scripts/mcmc_BDS_bears.Rev
