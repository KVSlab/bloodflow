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

q(0,t) = q_inlet(t)
q(L,t) = q(0,t), as the fluid is incompressible (?).
q(x,0) = q_inlet(0) since r_0 = constant, which makes the artery a perfect cylinder.

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

N = len(t)-1
M = 100

L = 20.8
Tf = t[-1]

dx = L/M            # constant
dt = t[1:] - t[:-1] # dt is not constant


r0 = 0.37
E = 1.0
H = 0.01
Re = 1.0
nu = 1.0
T = 1.0
db = np.sqrt(nu*T/2/np.pi)
A0 = np.pi*r0**2

f = 4*E*H/3/r0

#####################################################################

# We define F and S as functions of U, where U[0] = A and U[1] = q

def F(U):
	return np.array([U[1], U[1]**2 + np.sqrt(A0*U[0])])

def S(U):
	return np.array([0, -2*np.sqrt(np.pi)/db/Re*U[1]/np.sqrt(U[0])])
	
#####################################################################

# Initialisation
U = np.zeros([2,M+1,N+1])

U[0,:,0] = A0*np.ones(M+1)
U[0,0,:] = A0*np.ones(N+1)
U[0,-1,:] = A0*np.ones(N+1)

U[1,0,:] = q_inlet
U[1,-1,:] = q_inlet
U[1,:,0] = q_inlet[0]*np.ones(M+1)

# Time loop
for n in range(N):
	# Definitions for compact formulae
	dtn = dt[n]
	Un = U[:,:,n]
	Fn = F(Un)
	Sn = S(Un)
	# Half step
	U_half = (Un[:,1:]+Un[:,:-1])/2 - dtn/dx/2*(Fn[:,1:]-Fn[:,:1]) + dtn/4*(Sn[:,1:]+Sn[:,:1])
	F_half = F(U_half)
	S_half = S(U_half)
	#Complete step
	Un1 = Un[:,1:-1] - dtn/dx*(F_half[:,1:]-F_half[:,:-1]) + dtn/2*(S_half[:,1:]+S_half[:,:-1])
	U[:,1:-1,n+1] = Un1

plt.imshow(U)
plt.colorbar()
plt.savefig('flow.png')
