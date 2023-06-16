#!/usr/bin/env python3

import sys
import os
import time

n_sweeps = int(sys.argv[1])
n_seeds = int(sys.argv[2])

for i in range(n_seeds):
    #write name of file including number of sweeps and used seed
    fname = "sweeps_"+str(n_sweeps)+"_seed_"+str(i+1)+".sh"
    #create bash file
    os.system("touch "+fname)
    #write common text to file
    os.system("cat job_common_text.txt >> "+fname)
    #add to last line the unique seed
    line = "\"\njulia cluster_simulations.jl "+str(n_sweeps)+" "+str(i+1)+"\""
    os.system("echo "+line+">>"+fname)
    #call file from terminal
    os.system("bash "+fname)
    #wait 2 seconds before sending next job
    time.sleep(2)
