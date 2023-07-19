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

## CONSTANTS ##


## FUNCTIONS ##

def glv_functional_response(x, r, h, A):
    nspp = len(x)
    func_vec = np.zeros(nspp)
    for i in range(nspp):
        inter_term = np.zeros(nspp)
        for j in range(nspp):
            if j==i:
                inter_term[j] = 0
            else:
                inter_term[j] = A[i, j] * x[j]/(1+h[j]*x[j])
        func_vec[i] = r[i] - A[i,i] * x[i] - sum(inter_term)
    return func_vec

def jac_glv_functional_response_non_vec(x, r, h, A):
    nspp = len(x)
    J = np.zeros((nspp, nspp))
    for i in range(nspp):
        for j in range(nspp):
            if i==j:
                J[i,i] = -A[i, i]
            else:
                J[i,j] = -A[i,j]/(1+h[j]*x[j])**2
    return J

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
    nmax = 15
    #how many solutions to find
    nsols_found = int(1e4)
    #how many simulations to carry for each diversity
    nsim = 100
    #how many trials for each simulation
    ntrialsmax=2*nsols_found
    #preallocate space for number of equlibria for each diversity
    neq_results = np.zeros((0, 3))
    for i in range(nmax):
        for s in range(nsim):
            #initialize storing matrix for actual unique equilibria
            eq_results = np.zeros((0, nmax+2))
            nspp = i+1
            print("simulation: ", s, "\n")
            print("diversity: ", nspp, "\n")
            nsols = 0
            #sample parameters of system
            r = np.random.rand(nspp)
            h = np.random.rand(nspp)
            A = np.random.rand(nspp, nspp) 
            ntrials = 0
            #keep finding roots until nsols_found are found
            while nsols < nsols_found:
                x0 = np.random.normal(scale = 1e6, size = nspp)
                res = fsolve(glv_functional_response, x0, args = (r, h, A), \
                        fprime=jac_glv_functional_response_non_vec)
                res_rounded = np.round(res, decimals = 2)
                #store if its an actual solution
                if sum(abs(glv_functional_response(res, r, h, A))) < 1e-6:
                    #complete res with trailing zeros to store
                    res_expand = np.hstack((s, nspp, res_rounded, 
                        np.zeros(nmax-nspp)))
                    eq_results = store_if_unique(eq_results, res_expand)
                    nsols += 1
                #break out of the while loop after enough trials 
                ntrials += 1
                print("Number of trials ", ntrials, end = "\r")
                if (ntrials == nsols_found and nsols == 0) or \
                        (ntrials > ntrialsmax):
                    break
            nsol = eq_results.shape[0]
            npos = sum([all(eq_results[i,2:nspp]>0) for i in range(nsol)])
            add_row = np.array([nspp, nsol, npos])
            neq_results = np.vstack((neq_results, add_row))
            #create dataframe and save
            df = pd.DataFrame(neq_results, columns = ["n", "nsol", "npos"])
            df.to_csv("../data/neq_classic_models_rand.csv", index = False)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
