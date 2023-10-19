#!/bin/bash
#SBATCH --job-name=glv_type1
#SBATCH --account=pi-salesina
#SBATCH --output=glv_type1.out
#SBATCH --error=glv_type1.err
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

module load R/4.2.0
Rscript glv_type1.r

mv glv_type1.sbatch runfiles
mv glv_type1.* outputfiles
