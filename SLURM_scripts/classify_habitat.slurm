#!/bin/bash

#SBATCH --job-name=HabClass
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00
#SBATCH --output=/work/hobi/jselwyn/Habitat/SLURM_Out/habitat_classification_%j.out

#### Script to run the initial model fitting of habitat classification model.
#### Run classification_selection.R after selecting top model

module load R/gcc/64/3.5.1

Rscript /home/jselwyn/Habitat/classify_habitat.R
