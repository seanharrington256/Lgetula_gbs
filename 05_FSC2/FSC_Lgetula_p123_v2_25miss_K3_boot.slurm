#!/bin/bash

#SBATCH --job-name FSC_k5_Lget
#SBATCH -A inbreh
#SBATCH -p teton
#SBATCH -t 1-00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=seanharrington256@gmail.com
#SBATCH -e err_FSC_k5_Lget_%A_%a.err
#SBATCH -o std_FSC_k5_Lget_%A_%a.out
#SBATCH --array=1-5000

source ~/.bash_profile

cd /project/inbreh/fsc/boot_input/Reps   # set the correct working directory - make sure that the array number is correct, too

for x in */*; do
  repdirs=(${repdirs[@]} "${x}")
done

cd ${repdirs[($SLURM_ARRAY_TASK_ID-1)]}   ## here, I subtract 1 from the $SLURM_ARRAY_TASK_ID because bash indexing starts at zero - I think it's less confusing to subtract 1 here than to remember to do it when specifying the number of jobs for the array

est=$(echo *.est)
tpl=$(echo *.tpl)

fsc26 -t $tpl -e $est -n 100000 -L 40 -m -M -u -x -q
