#!/bin/bash

#SBATCH --job-name=inla
#SBATCH -p normal,cbirdq
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00
#SBATCH --output=/work/hobi/jselwyn/Habitat/SLURM_Out/run_inla_%j.out

#### Full Fitting inlabru model from param file ####

param_file=$1
echo $param_file

run_name=$2.${SLURM_ARRAY_TASK_ID}
echo $run_name

module load R/gcc/64/3.5.1

Rscript /home/jselwyn/Habitat/inlabru_model_v7.0.R ${param_file} ${run_name}
