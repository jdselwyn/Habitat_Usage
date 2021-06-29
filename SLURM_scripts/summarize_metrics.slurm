#!/bin/bash

#SBATCH --job-name=metricSummary
#SBATCH -p normal,gpu
#SBATCH --nodes=1
#SBATCH -t 1-00:00:00

#### Create all site mean/sd for z-score transformation

model=$1
echo $model

module load R/gcc/64/3.5.1

Rscript /home/jselwyn/Habitat/summarize_metrics.R ${model}
