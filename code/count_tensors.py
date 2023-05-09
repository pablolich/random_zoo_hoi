#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@uchicago.edu)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
import itertools as it

## CONSTANTS ##


## FUNCTIONS ##

def get_all_sign_flips(dim):
    return it.product(range(2), repeat = dim)

def binary_to_signed(vec):
    return(-2*vec + 1)

def get_new_B(B_old, D, dim):
    #make copy of B_old
    B_new = np.copy(B_old)
    for i in range(dim):
        B_new[i,:,:] = D@B_new[i,:,:]@D
    return B_new

def is_new(list_elements, new_element):
    n_elements = len(list_elements)
    for el in range(n_elements):
        equal = (list_elements[el] == new_element).all()
        if equal:
            return False
    return True
def count_Bs(n):
    #get total number of sign flips
    n_flips = 2**n
    #initialize storage
    list_Bs_unique = list()
    list_Bs = list()
    list_D = list()
    #sample random B
    B = abs(np.random.rand(n, n, n))
    #get all sign flip vectors (1 is sign flip, 0 is no flip)
    flips_iter = it.product(range(2), repeat = n)
    #loop through all sign flips
    for flip in flips_iter:
        #get vector of 1s and -1s
        vec = binary_to_signed(np.array(flip))
        #embed in the diagonal of a matrix 
        D_vec = np.zeros((n,n))
        np.fill_diagonal(D_vec, vec)
        #construct new B
        B_flip = get_new_B(B, D_vec, n)
        #store if it hasn't been seen before
        list_Bs.append(B_flip)
        list_D.append(D_vec)
        if is_new(list_Bs_unique, B_flip):
            list_Bs_unique.append(B_flip)
    return (list_Bs_unique)

def main(argv):
    '''Main function'''
    #number of species and total sign flips
    n_vec = np.array([3, 4, 5, 6]) 
    n_eq = list()
    for i in n_vec:
        print("order of hois: ", i)
        list_Bs = count_Bs(i)
        n_eq.append(len(list_Bs))
    print('Number of different Bs up to a sign change: ', n_eq)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
