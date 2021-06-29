#!/bin/bash

#SBATCH --job-name=reefTranslate
#SBATCH -p normal,gpu
#SBATCH --nodes=1
#SBATCH -t 1-00:00:00

#### Join all sites into a site mosaic

model=$1
echo $model

module load R/gcc/64/3.5.1

Rscript /home/jselwyn/Habitat/reefTranslate.R ${model}
