#!/bin/bash

#SBATCH --job-name=sitePred
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00

#### post-process & plot predictions from inlabru onto individual sites

run_name=$1
echo $run_name


readarray -t SITES < <(cut -d, -f2 /home/jselwyn/Habitat/site_list.txt) #Read in that individual indexing file
site=${SITES[${SLURM_ARRAY_TASK_ID}]} #Take the array IDth individual
echo $site

module load R/gcc/64/3.5.1

#Old version to run buffered version
# echo "#! /bin/bash" > /work/hobi/jselwyn/Habitat/buffer.sh
# echo "Rscript /home/jselwyn/Habitat/buffered_topography.R ${model}" >> /work/hobi/jselwyn/Habitat/buffer.sh
#
# echo "#! /bin/bash" > /work/hobi/jselwyn/Habitat/site.sh
# echo "Rscript /home/jselwyn/Habitat/site_topography.R ${model}" >> /work/hobi/jselwyn/Habitat/site.sh
#
# chmod 777 /work/hobi/jselwyn/Habitat/buffer.sh
# chmod 777 /work/hobi/jselwyn/Habitat/site.sh
#
# srun --nodes=1 --ntasks=1 /work/hobi/jselwyn/Habitat/buffer.sh &
# srun --nodes=1 --ntasks=1 /work/hobi/jselwyn/Habitat/site.sh #&
# wait

Rscript /home/jselwyn/Habitat/site_postProcess.R ${run_name} ${site}
