#!/bin/bash

#SBATCH --job-name=ModelChoice
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00

#### Script to choose model from models already fit using "classify_habitat.slurm"

model=$1
echo $model

readarray -t SITES < <(cut -d, -f2 /home/jselwyn/Habitat/site_list.txt) #Read in that individual indexing file
site=${SITES[${SLURM_ARRAY_TASK_ID}]} #Take the array IDth individual
echo $site

module load R/gcc/64/3.5.1

Rscript /home/jselwyn/Habitat/classification_selection.R ${model} ${site}
