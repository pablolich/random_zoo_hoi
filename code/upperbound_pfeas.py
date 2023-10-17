#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@uchicago.edu)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
from scipy.optimize import minimize, Bounds, LinearConstraint, linprog
import pdb, traceback

## FUNCTIONS ##
def lp_n_d(n, d):
    nroots = d**n+1 #there can be zero positive roots too
    obj = [1]+(nroots-1)*[0] #objective function to minimize
    lhs_eq = [nroots*[1], #equality constraint (left hand side)
              list(range(nroots))]
    rhs_eq = [1, np.sqrt(nroots)/2**n] #equality constraint (right hand side)
    bnd = nroots*[(0,1)] #bounds for p
    opt = linprog(c = obj, A_eq = lhs_eq, b_eq = rhs_eq, bounds = bnd, 
                  method = "highs")
    return opt

def main(argv):
    '''Main function'''
    #choose arguments
    n_vec = np.arange(8, dtype = int) + 1
    d_vec = np.arange(6, dtype = int) + 1
    n_n = len(n_vec)
    n_d = len(d_vec)
    #preallocate
    data = np.zeros((n_n*n_d, 3))
    data[:,0] = np.repeat(n_vec, n_d)
    data[:,1] = np.tile(d_vec, n_n)
    #iterator
    it = 0
    for n in n_vec:
        print("Diversity ", n)
        for d in d_vec:
            print("Interaction order ", d)
            result = lp_n_d(n, d)
            #store
            data[it, 2] = 1 - result.fun #get probability of feasibility bound
            it+=1
    np.savetxt("../data/p_feas_upper_bound.csv", data)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

