#!/bin/bash

##Run many inla models sbatch array

#bash Habitat/multi_inla.sh /home/jselwyn/Habitat/params/inlabru_params.R 9.7.21 1000 0 20

param_file=$1
echo $param_file

run_name=$2
echo $run_name

number_fits=$3
echo $number_fits

start_number=$4
echo $start_number

nodes=$5
echo $nodes

#Spawn array of jobs to fit many inlabru models
JOBID=$(sbatch --array=${start_number}-$((${number_fits}-1))%${nodes} --output=/work/hobi/jselwyn/Habitat/SLURM_Out/run_inla_%A_%a.out /home/jselwyn/Habitat/run_inla.slurm ${param_file} ${run_name})
NUMBER1=$(echo ${JOBID} | sed 's/[^0-9]*//g')


#post-process inlabru models
JOBID=$(sbatch --dependency=afterany:${NUMBER1} --output=/work/hobi/jselwyn/Habitat/SLURM_Out/post_process_inlabru_%j.out /home/jselwyn/Habitat/post_process_inlabru.slurm ${run_name})
NUMBER2=$(echo ${JOBID} | sed 's/[^0-9]*//g')

#Make post-processed tifs with predictions of means etc.
readarray -t SITES < <(cut -d, -f2 /home/jselwyn/Habitat/site_list.txt)
JOBID=$(sbatch --dependency=afterany:${NUMBER1} --array=0-$((${#SITES[@]}-1))%6 --output=/work/hobi/jselwyn/Habitat/SLURM_Out/site_postProcess_%A_%a.out /home/jselwyn/Habitat/site_postProcess.slurm ${run_name})
