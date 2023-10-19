#!/bin/bash
#SBATCH --job-name=glv_randomzoo
#SBATCH --account=pi-salesina
#SBATCH --output=glv_randomzoo.out
#SBATCH --error=glv_randomzoo.err
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

module load R/4.2.0
Rscript glv_randomzoo.r

mv glv_randomzoo.sbatch runfiles
mv glv_randomzoo.* outputfiles
