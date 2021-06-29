#!/bin/bash

#SBATCH --job-name=setup_INLAbru
#SBATCH -p normal,gpu
#SBATCH --nodes=1
#SBATCH -t 1-00:00:00

#### Script to fit inlabru models for each site using inlabru_model_vX.0.R

model=$1
echo $model

module load R/gcc/64/3.5.1

Rscript /home/jselwyn/Habitat/setup_inlabru.R ${model}
