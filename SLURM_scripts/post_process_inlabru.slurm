#!/bin/bash

#SBATCH --job-name=postINLA
#SBATCH -p normal,cbirdq
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00
#SBATCH --output=/work/hobi/jselwyn/Habitat/SLURM_Out/post_process_inlabru_%j.out

#### Post-Process many inlabru models ####

run_name=$1
echo $run_name

module load R/gcc/64/3.5.1

Rscript /home/jselwyn/Habitat/post_process_inlabru.R ${run_name}
