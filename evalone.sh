#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 1:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type FAIL,END
#
# File: evalone.sh
#
# Description: Runs evaluate_clustering.pl on a single step.
#
# Author: Yuki Koyanagi
# History:
#

echo Running on "$(hostname)"
echo Available nodes: "$SLURM_NODELIST"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"
start=$(date +%s)

# Load the modules previously used when compiling the application
module purge
module load python-intel

workdir=/work/austmathjea/cdp

echo Starting Python program
python ./evalone.py $workdir/step2000/n1/summary.pkl \
       2000 $SLURM_SUBMIT_DIR/eval.pl $workdir/opts

rm -f $SLURM_SUBMIT_DIR/step2000_summary2.txt
mv $SLURM_SUBMIT_DIR/step2000_assess $workdir/step2000/n1

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"
