#!/bin/bash

##Run many inla models sbatch array

param_file=$1
echo $param_file

run_name=$2
echo $run_name

number_fits=$3
echo $number_fits

start_number=$4
echo $start_number

JOBID=$(sbatch --array=${start_number}-$((${number_fits}-1))%10 --output=/work/hobi/jselwyn/Habitat/SLURM_Out/run_inla_%A_%a.out /home/jselwyn/Habitat/run_inla.slurm ${param_file} ${run_name})
NUMBER1=$(echo ${JOBID} | sed 's/[^0-9]*//g')
