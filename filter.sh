#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --ntasks-per-node 24      # number of MPI tasks per node
#SBATCH --time 1:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL
#
# File: filter.sh
#
# Author: Yuki Koyanagi
#

echo Running on "$(hostname)"
echo Available nodes: "$SLURM_NODELIST"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"
start=$(date +%s)

# Load the modules previously used when compiling the application
module purge
module load python-intel

# Start in total [nodes]*24 MPI ranks on all available CPU cores
srun python ./filter.py

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"
echo Done.
