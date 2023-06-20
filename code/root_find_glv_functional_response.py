#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@uchicago.edu)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import qmc

## CONSTANTS ##


## FUNCTIONS ##

def polynomial_form(x, r, h, A):
    nspp = len(x)
    fun_vec = np.zeros(nspp)
    for i in range(nspp):
        #create vector of functional response
        vec_func = 1+h*x
        #eliminate ith component from vector and take the product
        vec_func_jnoti = np.copy(vec_func)
        vec_func_jnoti[i] = 1
        prod_jnoti = np.prod(vec_func_jnoti)
        #preallocate interspecific competition term
        inter_comp = np.zeros(nspp)
        for j in range(nspp):
            if j == i:
                inter_comp[j] = 0
            else:
                #eliminate jth component from vector and take the product
                vec_func_knotj = np.copy(vec_func_jnoti)
                vec_func_knotj[j] = 1
                prod_knotj = np.prod(vec_func_knotj)
                #build component j
                inter_comp[j] = A[i,j]*x[j]*prod_knotj
        tot_inter_comp_i = sum(inter_comp)
        fun_vec[i] = r[i]*prod_jnoti - A[i,i]*x[i]*prod_jnoti - tot_inter_comp_i
    return fun_vec

##############################################################################
#test function: it works for at most 3 spp
#nspp = 2
#r = np.array([1/2, 2, 1/3])
#h = np.array([1,2,1])
#A = np.array([[2, 1, 1], [1,2, 1], [1, 1, 2]])
##equilibrium of this model is easy to compute in mathematica:
#x_star_analytically = np.array([0.110015,0.975757, -0.0481869])
#fun_value_test = fun(x_star_analytically, r, h, A)
#print(fun_value_test) #it should be close to 0
##############################################################################

def classic_form(x, r, h, A):
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

##############################################################################
#test function: it works for at most 3 spp
#r = np.array([1/2, 2, 1/3])
#h = np.array([1,2,1])
#A = np.array([[2, 1, 1], [1,2, 1], [1, 1, 2]])
##equilibrium of this model is easy to compute in mathematica:
#x_star_analytically = np.array([-1.61555, -1.83414, -1.4893])
#fun_value_test = classic_form(x_star_analytically, r, h, A)
#print(fun_value_test) #it should be close to 0
#import ipdb; ipdb.set_trace(context = 20)
##############################################################################

def fun(x, r, h, A):
    return sum(abs(classic_form(x, r, h, A)))


#minimize all components of fun_vec over different x, first withouth inputing
#the gradient

#sample x0

def main(argv):
    '''Main function'''
    #set random seed
    np.random.seed(2)
    #max diversity
    nmax = 5
    #how many solutions to find
    nsols_found = 50
    #preallocate results
    neq_results = np.zeros((nmax, 2))
    for i in range(nmax):
        nspp = i+1
        print("Diversity: ", nspp)
        nsols = 0
        #sample parameters of system
        r = np.random.rand(nspp)
        h = np.random.rand(nspp)
        A = np.random.rand(nspp, nspp)
        #sample many initial conditions
        sampler = qmc.LatinHypercube(d = nspp)
        results = np.zeros((nsols_found, nspp))
        #sampling bounds
        lb = nspp*[-1e6]
        ub = nspp*[1e6]
        #keep minimizing until 1000 solutions are found
        while nsols < nsols_found:
            sample = sampler.random(n=1)
            #scale to bounds
            x0 = qmc.scale(sample, lb, ub)[0]
            res = minimize(fun, x0, args = (r, h, A), 
                    method="Nelder-Mead", tol = 1e-20)
            #only save if its a solution
            if res.fun == 0:
                results[nsols,:] = res.x
                nsols += 1
                print("Number of equilibria found ", nsols)
        #get unique equilibria
        unique_eq = np.unique(results, axis = 0)
        #count how many
        neq = unique_eq.shape[0]
        #add to neq_results
        neq_results[i, 0] = nspp
        neq_results[i, 1] = neq
    import ipdb; ipdb.set_trace(context = 20)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

