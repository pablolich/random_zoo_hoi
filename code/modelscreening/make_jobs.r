#file to create jobs for cluster based on names of files to run
args = commandArgs(trailingOnly = TRUE)
njobs = length(args)
for (i in 1:njobs){
  #name the job and output files
  jobname = paste0(args[[i]])
  output = paste0(jobname, ".out")
  error = paste0(jobname, ".err")
  #create a text file for ith job
  sink(paste0(jobname, ".sbatch"))
  cat(paste0("#!/bin/bash",
         "\n#SBATCH --job-name=", jobname,
         "\n#SBATCH --account=pi-salesina",
         "\n#SBATCH --output=", output,
         "\n#SBATCH --error=", error,
         "\n#SBATCH --time=10:00:00",
         "\n#SBATCH --partition=caslake",
         "\n#SBATCH --nodes=1",
         "\n#SBATCH --ntasks-per-node=1",
         "\n#SBATCH --mem-per-cpu=2000\n",
         "\nmodule load R/4.2.0",
         "\nRscript ", jobname,".r\n",
         "\nmv ", jobname, ".sbatch runfiles",
         "\nmv ", jobname, ".* outputfiles",
         "\nmv outputfiles/", jobname, ".r .\n"))
  sink()
}