#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim
#SBATCH --nodes 3                 # number of nodes
#SBATCH --time 2:00:00            # max time (HH:MM:SS)
#
# File: pred_prll.sh
#
# Description: This is the job script submitted by sbatch for
#  prediction. See pred_prll.py for arguments expected by it.
#
# Author: Yuki Koyanagi
#

echo Running on "$(hostname)"
echo Running job: "$SLURM_JOB_NAME"
echo Available nodes: "$SLURM_NODELIST"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"
start=$(date +%s)

echo Writing hostnames
scontrol show hostnames $SLURM_NODELIST > /tmp/nodelist

echo Enable modules
module purge
module add python-intel pp

echo Starting servers
srun ppserver.py -p 2048 -t 30 &
sleep 1 # sleep a bit to ensure that the servers have started

echo Starting Python program
python predict.py /work/austmathjea/cdp/n1assess \
       /work/austmathjea/cdp/test \
       100 /work/austmathjea/cdp/prediction.txt

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"

sleep 31 # this should kill all remote ppservers

echo Done.
