import numpy as np
import matplotlib.pyplot as plt

'''I use the Rayleigh-Ritz method to solve a second order  ordinary differential equation with a nonzero RHS
The d.e is of the form :  - d/dx(  p(x)dy/dx) + q(x)y = f(x).
The integration is performed on the interval 0<=x<=1'''

def rhs(x):
    '''This computes the rhs of the d.e.
    INPUT: the variable x'''
    return np.sin(np.pi*x)
N =5
x = 2*np.pi*np.linspace(0,1, N)


h = []          #intialize space step
for i in range(N-1):
        h.append(x[i+1] - x[i])

def basis_function(x):
    '''Basis functions of the Rayleigh-Ritz method'''
                                    
    phi = np.zeros([N, N], dtype='float')                               #initialize basis function
    for j in range(1, N-1):             
        for idx in range(1, N-1):
            if (0<=x[idx]<=x[j-1]):
                phi[j, idx] = 0
            elif (x[j-1]<x[idx]<=x[j]):
                phi[j, idx] = (x[idx] - x[j-1])/h[j-1]
            elif (x[j]<x[idx]<=x[j+1]):
                phi[j, idx] = ((x[j+1] - x[idx])/h[j])
            elif (x[j+1]<x[idx]<=1):
                phi[j, idx] = 0
    return phi
#basis_function(x)

#coefficient functions 
def p(x):
    '''Defining the coefficient p in the d.e.'''
    return 1

def q(x):
    '''Defining the coefficient q in the d.e.'''
    return np.pi**2
#initialize vectors for computing integrals
Q1 = []
Q2 = []
Q3 = []
Q4 = []
Q5 = []
Q6 = []
#approximating the 6 integrals
for i in range(N-1):
    Q1.append((h[i]/12)*(p(x[i]) + q(x[i+1])))
    Q2.append((h[i-1]/12)*(3*q(x[i]) + q(x[i-1])))
    Q3.append((h[i])*(3*q(x[i]) + q(x[i+1])))
    Q4.append((h[i-1]/2)*(p(x[i]) + q(x[i-1])))
    Q5.append((h[i-1]/6)*(2*rhs(x[i]) + rhs(x[i-1])))
    Q6.append((h[i]/6)*(2*rhs(x[i]) + rhs(x[i+1])))
    #now compute  Q4,n 
    Q4.append((h[N-2]/2)*(p(x[N-1]) + q(x[N-2])))

alpha = np.zeros(N, dtype='float')
beta = np.zeros(N, dtype='float')
b = np.zeros(N, dtype='float')
a = np.zeros(N, dtype='float')
for i in range(N-2):
    alpha[i] = Q4[i] + Q4[i+1] + Q2[i] + Q3[i]
    beta[i] = Q1[i] - Q4[i+1]
    b[i] = Q5[i] + Q6[i]

print(len(Q4))
alpha[N-1] = Q4[N-1] + Q4[N] + Q2[N-1] + Q3[N-1]



