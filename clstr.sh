#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim
#SBATCH --nodes 32                 # number of nodes
#SBATCH --time 6:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL
#
# Time-stamp: <2017-02-24 10:26:27 au447708>
#
# File: clstr.sh
#
# Description: This is the job script submitted by sbatch for
#  clustering analysis. We estimate 16 nodes take ca. 10min to run
#  clustering for a single step. Control the step numbers for which
#  we run the clustering by setting -J option in sbatch. Argument
#  in the form {name}.stepstart.stepend, or {name}-group#-#ofgroups,
#  where group# is 0-based.
#
# Author: Yuki Koyanagi
# History:
#

echo Running on "$(hostname)"
echo Running job: "$SLURM_JOB_NAME"
echo Available nodes: "$SLURM_NODELIST"
echo Number of nodes: "$SLURM_NNODES"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"
start=$(date +%s)

echo Writing hostnames
scontrol show hostnames $SLURM_NODELIST > /tmp/nodelist

echo Enable modules
module purge
module add python-intel pp R

echo Starting servers
srun ppserver.py -p 2048 -k 36000 -t 30 &
sleep 1 # sleep a bit to ensure that the servers have started

echo Starting Python program
python clstr.py

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"

sleep 31 # this should kill all remote ppservers

echo Done.
