#!/bin/bash

#SBATCH --job-name=compress
#SBATCH -p normal,cbirdq
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00
#SBATCH --output=/work/hobi/jselwyn/Habitat/SLURM_Out/compress_%j.out

#### compress files ####

cd /work/hobi/jselwyn/Habitat/Results/INLA_models

modelPrefix=models_17.5.21

dirList=$(ls -d ${modelPrefix}*)
tar -cjvf ${modelPrefix}.tar.bz ${dirList}
rm -r ${dirList}
