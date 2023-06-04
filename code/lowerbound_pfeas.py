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

## CONSTANTS ##


## FUNCTIONS ##

def get_c_linprog(n, d):
    k = np.arange(d**n+1)
    return (1-1/2**n)**k

def min_pfeas_linprog(n, d):
    '''
    Minimize probability of feasibility using linear programming
    '''
    nsol_max = int(d**n)
    #coefficients of linear problem
    c = get_c_linprog(n, d)
    #equality constraints
    A = np.vstack([np.ones(nsol_max),    #probability constraint
                   np.arange(nsol_max), #expected value constraint
                   np.arange(1)]) #add variance
    b = np.array([1, np.sqrt(nsol_max)])   #equality constraint values
    #bounds
    upper_bound = parity_restricted_bounds(d, nsol_max)
    lower_bound = np.zeros(nsol_max)
    bounds = [tuple(i) for i in zip(upper_bound, lower_bound)]
    res = linprog(c, A_eq = A, b_eq = b, bounds=bounds)
    return res

def fun(p, n, d):
    #vector of posible number of real solutions
    ind_sum = np.arange(d**n+1)
    fun_value = sum(p*(1-1/2**n)**ind_sum)
    return fun_value

def parity_restricted_bounds(d, nsol_max):
    #vector of bounds
    upper_bounds = np.ones(nsol_max+1)
    if (d%2)==0:
        odd_indices = np.arange(1, nsol_max, 2)
        #fix to 0 odd probabilities
        upper_bounds[odd_indices] = 0
    else:
        even_indices = np.arange(0, nsol_max, 2)
        #fix to 0 even probabilities
        upper_bounds[even_indices] = 0
    return upper_bounds

def minimize_pfeas(fun, p0, nsol_max, n, d, var):
    #deal with the d=2 case separately
    if d == 1:
        bounds = Bounds(lb = 0, ub = 1)
    else:
        #bounds for p
        upper_bound = list(parity_restricted_bounds(d, nsol_max))
        lower_bound = list(np.zeros(nsol_max+1))
        bounds = Bounds(lower_bound, upper_bound)
    #get matrix of constraints
    A = np.vstack([np.ones(nsol_max+1),    #probability constraint
                   np.arange(nsol_max+1), #expected value constraint
                   (np.arange(nsol_max+1)-np.sqrt(nsol_max))**2]) #variance
    moments = np.array([1, np.sqrt(nsol_max), var])
    #create linear constraints
    constraints = LinearConstraint(A, 
                                   lb = moments,
                                   ub = moments)
    #minimize
    res = minimize(fun, p0, args = (n, d), 
                   method = 'trust-constr',
                   bounds = bounds,
                   constraints = [constraints]) 
    return res


def main(argv):
    '''Main function'''
    #choose arguments
    n_vec = np.arange(8, dtype = float) + 1
    d_vec = np.arange(5, dtype = float) + 2
    n_n = len(n_vec)
    n_d = len(d_vec)
    #load variance for each n, d
    df = pd.read_csv("../data/variance_data.csv", delimiter = " ")
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
            nsol_max = int(d**n) #count 0
            #get variance
            df_sub = df[(df["d"] == d+1) & (df["n"] == n)]
            var = float(df_sub["var_real"])
            p0 = np.random.dirichlet(np.ones(nsol_max+1))
            res = minimize_pfeas(fun, p0, nsol_max, n, d, var)
            #store
            data[it, 2] = res.fun
            data[it, 1] +=1
            it+=1
        #add one to d
        #save
        np.savetxt("../data/p_feas_min_d.csv", data)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
