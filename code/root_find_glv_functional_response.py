#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@uchicago.edu)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd

## CONSTANTS ##


## FUNCTIONS ##

def fun(x, r, h, A):
    nspp = len(x)
    fun_vec = np.zeros(nspp)
    for i in range(nspp):
        #create vector of functional response
        vec_func_i = 1+h[i]*x
        vec_func_i_not_j = np.delete(vec_func_i, i)
        prod_inotj = np.prod(vec_func_i_not_j)
        inter_comp = np.zeros(nspp)
        for j in range(nspp):
            if j==i:
                inter_comp[j] = 0
            else:
                vec_func_j_not_k = np.delete(vec_func_i, j)
                prod_jnotk = np.prod(vec_func_j_not_k)
                inter_comp[j] = A[i,j]*x[j]*prod_jnotk
        tot_inter_comp_i = sum(inter_comp)
        fun_vec[i] = r[i]*prod_inotj - A[i,i]*x[i]*prod_inotj - tot_inter_comp_i
    return fun_vec

#test function
nspp = 2
r = np.ones(nspp)
h = np.ones(nspp)
A = np.ones((nspp, nspp))
#equilibrium of this model is easy to compute:
x_star_analytically = np.array([0.6180339887498949, 0.618033988749894])
fun_value_test = fun(x_star_analytically, r, h, A)

def main(argv):
    '''Main function'''

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

