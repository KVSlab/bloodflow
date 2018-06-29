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
M = 5

L = 20.8
Tf = t[-1]

dx = L/M            # constant
dt = t[1:] - t[:-1] # dt is not constant

x = np.linspace(0,L,M+1)

r0 = 0.37
E = 1.0e+6
H = 0.01
Re = 1.0
nu = 0.04
T = 1.0
db = np.sqrt(nu*T/2/np.pi)
A0 = np.pi*r0**2

f = 4*E*H/3/r0

#####################################################################

# We define F and S as functions of U, where U[0] = A and U[1] = q

def F(A,q):
	return np.array([q, q**2 + f*np.sqrt(A0*A)])

def S(A,q):
	return np.array([0*A, -2*np.sqrt(np.pi)/db/Re*q/np.sqrt(A)])
	
#####################################################################

# Initialisation
A = np.zeros([M+1,N+1])
q = np.zeros([M+1,N+1])

A[:,0] = A0*np.ones(M+1)
A[0,:] = A0*np.ones(N+1)
A[-1,:] = A0*np.ones(N+1)

q[0,:] = q_inlet
q[-1,:] = q_inlet
q[:,0] = q_inlet[0]*np.ones(M+1)

# Time loop
for n in range(N):

	# Definitions for compact formulae
	dtn = dt[n]
	An = A[:,n]
	qn = q[:,n]
	Fn = F(An,qn)
	Sn = S(An,qn)
	
	# Half step
	A_half = (An[1:]+An[:-1])/2 - dtn/dx/2*(Fn[0,1:]-Fn[0,:1]) + dtn/4*(Sn[0,1:]+Sn[0,:1])
	q_half = (qn[1:]+qn[:-1])/2 - dtn/dx/2*(Fn[1,1:]-Fn[1,:1]) + dtn/4*(Sn[1,1:]+Sn[1,:1])
	F_half = F(A_half, q_half)
	S_half = S(A_half, q_half)
	
	#Complete step
	An1 = An[1:-1] - dtn/dx*(F_half[0,1:]-F_half[0,:-1]) + dtn/2*(S_half[0,1:]+S_half[0,:-1])
	qn1 = qn[1:-1] - dtn/dx*(F_half[1,1:]-F_half[1,:-1]) + dtn/2*(S_half[1,1:]+S_half[1,:-1])
	A[1:-1,n+1] = An1
	q[1:-1,n+1] = qn1
	
fig1 = plt.figure()
plt.imshow(A)
plt.colorbar()
plt.savefig('area.png')

fig2 = plt.figure()
plt.imshow(q)
plt.colorbar()
plt.savefig('flow.png')
