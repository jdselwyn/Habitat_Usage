#!/bin/bash

#SBATCH --job-name=ModelChoiceSummary
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00

#### Script to choose make all site plots and summaries from classification_selection.slurm

model=$1
echo $model

module load R/gcc/64/3.5.1

Rscript /home/jselwyn/Habitat/summarize_habitat_classification.R ${model}
