"""
Implementation of the 1D blood flow equations
using Richtmyer's two step Lax-Wendroff method
Author : Syver Døving Agdestein
--------------------------------------

dU/dt + dF/dx = S

             | A(x,t) |                |       q            |                  |       0        |
U = U(x,t) = |        |,    F = F(U) = |                    |,      S = S(U) = |                |
             | q(x,t) |                | q²/A + sqrt(A_0*A) |                  | -2πR/ð_bRe q/A |

A = πR²

r_0 = constant

--------------------------------------
"""

####################################################################

# Import useful libraries

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

####################################################################

# Data for our problem, considered global constants that can be used within functions

data_q = np.genfromtxt('../data/example_inlet.csv', delimiter = ',')
#print(data_q)
#plt.plot(data_q[:,0], data_q[:,1])
#plt.savefig('data.png')

t = data_q[:,0]
q_inlet = data_q[:,1]

N = len(t)
M = 100

L = 20.8
Tf = t[-1]

dx = L/M # Constant
# dt is not constant.


r0 = 0.37
E = 1.0
H = 0.01
Re = 1.0
nu = 1.0
T = 1.0
db = np.sqrt(nu*T/2/np.pi)
A0 = pi*r0**2

f = 4*E*H/3/r0

#####################################################################

# We define F and S as functions of U, where U[0] = A and U[1] = q

def F(U):
	return np.array([U[1], U[1]**2 + np.sqrt(A0*U[0])])

def S(U):
	return np.array([0, -2*np.sqrt(np.pi)/db/Re*U[1]/np.sqrt(U[0])])
	
#####################################################################

# We want to store the discrete solution U_m^n, 0<=m<=M, 0<=n<=N not as a matrix, but as a vector containing all the columns of the matrix. This allows us to write the passage from time step n to n+1 as an affine transformation.

# Returns the index k of the vector associated with the matrix (U_m^n)mn
def indk(m,n):
	return m + (M+1)*n

# Returns the index mn of the matrix associated with the vector (U_k)k
def indmn(k):
	return k - k%(M+1, k%(M+1)

# Returns the vector form of the matrix (U_m^n)mn
def vect(U):
	vectU = np.zeros((M+1)*(N+1))
	for n in range(N+1):
		for m in range(M+1):
			vectU[indk(m,n)] = U[m,n]
	return vectU

# Returns the matrix for of the vector (U_k)k
def mat(U):
	matU = np.zeros([(M+1),(N+1)])
	for n in range(N+1):
		for m in range(M+1):
			matU[m,n] = U[indk(m,n)]
	return matU

####################################################################

# We write 




