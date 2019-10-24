#!/usr/bin/env bash

#SBATCH --account austmathjea_slim
#SBATCH --time=01:00:00
#SBATCH --output=logs/parallel.log
#SBATCH --ntasks=48
#SBATCH --exclusive
#SBATCH --mail-type=END,Fail

module load parallel

# the --exclusive to srun makes srun use distinct CPUs for each job step
# -N1 -n1 allocates a single core to each task
srun="srun --exclusive -N1 -n1"

# --delay .2 prevents overloading the controlling node
# -j is the number of tasks parallel runs so we set it to $SLURM_NTASKS
# --joblog makes parallel create a log of tasks that it has already run
# --resume makes parallel use the joblog to resume from where it has left off
# the combination of --joblog and --resume allow jobs to be resubmitted if
# necessary and continue from where they left off
parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog logs/runtask.log --resume"

# this runs the parallel command we want
# in this case, we are running a script named runtask
# parallel uses ::: to separate options. Here {0..99} is a shell expansion
# so parallel will run the command passing the numbers 0 through 99
# via argument {1}
wdir=/work/austmathjea/cdp
ls ${wdir}/test/*.txt \
    | ${parallel} ${srun} python pred_clst.py {} ${wdir}/assess \
		  ${wdir}/opts > log_file
#$parallel "$srun python pk.py {1} 100 > logs/parallel_{1}.log" ::: {0..99}
