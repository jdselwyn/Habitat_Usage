#!/bin/bash

#SBATCH --job-name=compress
#SBATCH -p normal,cbirdq
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00
#SBATCH --output=/work/hobi/jselwyn/Habitat/SLURM_Out/compress_%j.out

#### compress files ####

cd /work/hobi/jselwyn/Habitat/Results/INLA_models

#sbatch clean.slurm models_13.3.21
#sbatch clean.slurm models_7.7.21
#sbatch clean.slurm models_14.3.21
#sbatch clean.slurm models_16.3.21

#modelPrefix=models_8.7.21
modelPrefix=$1

dirList=$(ls -d ${modelPrefix}*)
tar -cjvf ${modelPrefix}.tar.bz ${dirList}
rm -rf ${dirList}
