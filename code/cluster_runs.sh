#!/bin/bash
#SBATCH --job-name=trialrun
#SBATCH --account=pi-salesina
#SBATCH --output=trialrun.out
#SBATCH --error=trialrun.out
#SBATCH --time=24:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000


module load julia/1.9

julia cluster_simulations.jl 1 1
