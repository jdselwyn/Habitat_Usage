#!/bin/bash

#Run array across all sites listed in site_list.txt to calculate topographic measures

model=$1
readarray -t SITES < <(cut -d, -f2 /home/jselwyn/Habitat/site_list.txt)

#Calculate topography
JOBID=$(sbatch --array=0-$((${#SITES[@]}-1))%6 --output=/work/hobi/jselwyn/Habitat/SLURM_Out/topography_%A_%a.out /home/jselwyn/Habitat/topography.slurm ${model})
NUMBER1=$(echo ${JOBID} | sed 's/[^0-9]*//g')

#Summarize Topography
JOBID=$(sbatch --dependency=afterany:${NUMBER1} --output=/work/hobi/jselwyn/Habitat/SLURM_Out/metricSummary_%j.out /home/jselwyn/Habitat/summarize_metrics.slurm ${model})
NUMBER2=$(echo ${JOBID} | sed 's/[^0-9]*//g')

#Join everything together for topography used in inlabru
JOBID=$(sbatch --dependency=afterany:${NUMBER1} --output=/work/hobi/jselwyn/Habitat/SLURM_Out/reefTranslate_%j.out /home/jselwyn/Habitat/reefTranslate.slurm ${model})
NUMBER3=$(echo ${JOBID} | sed 's/[^0-9]*//g')
#sbatch --output=/work/hobi/jselwyn/Habitat/SLURM_Out/reefTranslate_%j.out /home/jselwyn/Habitat/reefTranslate.slurm c5

#Set up everything for inlabru except mesh & fish positions
JOBID=$(sbatch --dependency=afterany:${NUMBER3} --output=/work/hobi/jselwyn/Habitat/SLURM_Out/setup_inlabru_%j.out /home/jselwyn/Habitat/setup_inlabru.slurm ${model})
NUMBER4=$(echo ${JOBID} | sed 's/[^0-9]*//g')
