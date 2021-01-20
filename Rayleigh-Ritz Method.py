import numpy as np
import matplotlib.pyplot as plt

'''I use the Rayleigh-Ritz method to solve a second order  ordinary differential equation with a nonzero RHS
The d.e is of the form :  - d/dx(  p(x)dy/dx) + q(x)y = f(x).
The integration is performed on the interval 0<=x<=1'''

def rhs(x):
    '''This computes the rhs of the d.e.
    INPUT: the variable x'''
    return np.cos(x)
N =100
x = np.linspace(0,1, N)

def basis_function(x):
    '''Basis functions of the Rayleigh-Ritz method'''
    h = []                                  #intialize space step
    phi = np.zeros([N, N], dtype='float')                  #initialize basis function
    for i in range(N-1):
        h.append(x[i+1] - x[i])
    for j in range(1, N-1):
        for idx in range(1, N-1):
            if (0<=x[idx]<=x[j-1]):
                phi[j, idx] = 0
            elif (x[j-1]<=x[idx]<=x[j]):
                phi[j, idx] = (x[idx] - x[idx-1])/h[idx-1]
            elif (x[j]<=x[idx]<=x[j+1]):
                phi[j, idx] = ((x[idx+1] - x[idx])/h[idx])
            elif (x[j+1]<x[idx]<1):
                phi[j, idx] = 0
    return phi
#basis_function(x)

#coefficient functions 
def p(x):
    '''approximating the coefficient p in the d.e.'''
    pass

def q(x):
    '''approximating the coefficient q in the d.e.'''
    pass
#initialize vectors for computing integrals
Q1 = np.zeros(N, dtype='float')
Q2 = np.zeros(N, dtype='float')
Q3 = np.zeros(N, dtype='float')
Q4 = np.zeros(N, dtype='float')
Q5 = np.zeros(N, dtype='float')
Q6 = np.zeros(N, dtype='float')
#approximating the 6 integrals
for i in range(N-1):
    pass

