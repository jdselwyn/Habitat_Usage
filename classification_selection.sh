#!/bin/bash

#Run array across all sites listed in site_list.txt to predict sand or reef habitat type

model=$1
readarray -t SITES < <(cut -d, -f2 /home/jselwyn/Habitat/site_list.txt)

JOBID=$(sbatch --array=0-$((${#SITES[@]}-1))%6 --output=/work/hobi/jselwyn/Habitat/SLURM_Out/classification_selection_%A_%a.out /home/jselwyn/Habitat/classification_selection.slurm ${model})

NUMBER=$(echo ${JOBID} | sed 's/[^0-9]*//g')

#add in summary function call

sbatch --dependency=afterany:${NUMBER} --output=/work/hobi/jselwyn/Habitat/SLURM_Out/classification_summary_%j.out /home/jselwyn/Habitat/summarize_habitat_classification.slurm ${model}
