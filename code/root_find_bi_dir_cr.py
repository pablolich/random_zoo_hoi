import sys
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from scipy.stats import qmc
import time

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

def jac_glv_functional_response(x, r, h, A):
    nspp = len(x)
    J = np.zeros((nspp, nspp))
    for i in range(nspp):
        for j in range(nspp):
            if i==j:
                J[i,i] = -A[i, i]
            else:
                J[i,j] = -A[i,j]/(1+h[j]*x[j])**2
    return J
