#!/bin/bash

#SBATCH --job-name=inla
#SBATCH -p normal,cbirdq,gpu
#SBATCH --nodes=1
#SBATCH -t 1-00:00:00

#### Full Fitting inlabru model from param file ####

param_file=$1
echo $param_file

run_name=$2.${SLURM_ARRAY_TASK_ID}
echo $run_name

module load R/gcc/64/3.5.1

Rscript /home/jselwyn/Habitat/inlabru_model_v7.0.R ${param_file} ${run_name}
