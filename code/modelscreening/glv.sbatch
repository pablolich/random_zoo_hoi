#!/bin/bash
#SBATCH --job-name=glv
#SBATCH --account=pi-salesina
#SBATCH --output=glv.out
#SBATCH --error=glv.err
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

module load R/4.2.0
Rscript glv.r

mv glv.sbatch runfiles
mv glv.* outputfiles
