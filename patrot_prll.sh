#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim
#SBATCH --nodes 32                 # number of nodes
#SBATCH --time 6:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type=FAIL,END
#
# File: patrot_prll.sh
#
# Description: This is the job script submitted by sbatch for
#  loacalpattern analysis. Control steps to analyse by setting -J option;
#  {name}-startstep-endstep. If -J option is not set, it uses the max_step
#  set in the config.yaml.file.
#
# Author: Yuki Koyanagi
#

echo Running on "$(hostname)"
echo Running job: "$SLURM_JOB_NAME"
echo Available nodes: "$SLURM_NODELIST"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"
start=$(date +%s)

# echo Clearing scratch folder
# rm -f ${SCRATCH}/*

echo Writing hostnames
scontrol show hostnames $SLURM_NODELIST > /tmp/nodelist

echo Enable modules
module purge
module add python-intel pp

echo Starting servers
srun ppserver.py -p 2048 -k 36000 &
sleep 1 # sleep a bit to ensure that the servers have started

echo Starting Python program
python patrot_prll.py

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"

sleep 31 # this should kill all remote ppservers

echo Done.
