#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@uchicago.edu)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from scipy.stats import qmc
import time

def consumer_resource(x, d, mu, kappa, C):
    """
    x (nx1 array): species abundance
    d (nx1 array): species death rates
    mu (mx1 array): resource dilution rate
    kappa (mx1 array): resource supply rate
    C (mxn array): resource preference matrix
    """
    Cx = np.dot(C, x)
    w = kappa/(mu + Cx)
    return np.dot(C.T, w) - d

def store_if_unique(stored, tostore):
    nrow = stored.shape[0]
    comparison_vec = np.zeros(nrow)
    #compare tostore with each row of stored
    comparison_vec = [all(stored[i,:] == tostore) for i in range(nrow)]
    #not add if already there
    if (nrow == 0 or not any(comparison_vec)):
        return np.vstack((stored, tostore))
    else:
        return stored

#minimize all components of fun_vec over different x, first withouth inputing
#the gradient

def main(argv):
    '''main function'''
    #set random seed
    np.random.seed(2)
    #max diversity
    nsppmax = 10
    nresmax = 10
    #how many solutions to find
    nsols_found = int(1e4)
    #how many simulations to carry for each diversity
    nsim = 100
    #how many trials for each simulation
    ntrialsmax=2*nsols_found
    #preallocate space for number of equlibria for each diversity
    neq_results = np.zeros((0, 4))
    for i in range(nresmax):
        nres = i+1
        for j in range(i):
            nspp = j+1
            for s in range(nsim):
                #initialize storing matrix for actual unique equilibria
                eq_results = np.zeros((0, nsppmax+3)) #res + spp + sim = 3
                print(" ")
                print("resources: ", nres, "\n")
                print("diversity: ", nspp, "\n")
                print("simulation: ", s, "\n")
                print(" ")
                nsols = 0
                #sample parameters of system
                d = np.random.rand(nspp)
                mu = np.random.rand(nres)
                kappa = np.random.rand(nres)
                C = np.random.rand(nres, nspp) 
                ntrials = 0
                #keep finding roots until nsols_found are found
                while nsols < nsols_found:
                    x0 = np.random.normal(scale = 1e6, size = nspp)
                    res = fsolve(consumer_resource, x0, \
                            args = (d, mu, kappa, C))
                    res_rounded = np.round(res, decimals = 2)
                    #store if its an actual solution
                    if sum(abs(consumer_resource(res, d, mu, kappa, C))) < 1e-6:
                        #complete res with trailing zeros to store
                        res_expand = np.hstack((s, nres, nspp, res_rounded, \
                            np.zeros(nsppmax-nspp)))
                        eq_results = store_if_unique(eq_results, res_expand)
                        nsols += 1
                    #break out of the while loop after enough trials 
                    ntrials += 1
                    print("Number of trials ", ntrials, end = "\r")
                    if (ntrials == nsols_found and nsols == 0) or \
                            (ntrials > ntrialsmax):
                        break
                nsol = eq_results.shape[0]
                npos = sum([all(eq_results[i,3:nspp]>0) for i in range(nsol)])
                add_row = np.array([nres, nspp, nsol, npos])
                neq_results = np.vstack((neq_results, add_row))
                #create dataframe and save
                df = pd.DataFrame(neq_results, 
                    columns = ["m", "n", "nsol", "npos"])
                df.to_csv("../data/neq_cr_.csv", index = False)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
