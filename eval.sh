#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --ntasks-per-node 24      # number of MPI tasks per node
#SBATCH --time 0:30:00            # max time (HH:MM:SS)
#SBATCH --mail-type FAIL,END
#
# File: evalab.sh
#
# Description: Runs evaluate_clustering.pl on Abacus cluster.
#  The cutoff value can be specified by means of $SLURM_JOB_NAME,
#  by setting -J option when submitting. The format is {name}-cutoff
#  This took ??min for n=5.
#
# Author: Yuki Koyanagi
# History:
#

echo Running on "$(hostname)"
echo Available nodes: "$SLURM_NODELIST"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"

# Load the modules previously used when compiling the application
module purge
module load python-intel

# Start in total [nodes]*24 MPI ranks on all available CPU cores
srun python ./eval.py

echo Done: "$(date)"
