"""
Implementation of the 1D blood flow equations.
Author : Syver Døving Agdestein
"""

from fenics import *
import numpy as np
import matplotlib.pyplot as plt

data_q = np.genfromtxt('../data/example_inlet.csv', delimiter = ',')
#print(data_q)
#plt.plot(data_q[:,0], data_q[:,1])
#plt.savefig('data.png')

tt = data_q[:,0]
qq = data_q[:,1]

L, Tfinal = 20.8, data_q[-1,0]
Nx, Nt = 100, len(data_q[:,0])

r0 = 0.37
E = 1.0
H = 0.01
Re = 1.0
nu = 1.0
T = 1.0
deltab = np.sqrt(nu*T/2/np.pi)
A0 = pi*pow(r0,2)

f = 4*E*H/3/r0

mesh = RectangleMesh(Point(0,0), Point(L, Tfinal), Nx, Nt)
mesht = IntervalMesh(Nt,0,Tfinal)
V = FunctionSpace(mesh, 'P', 1)
W = FunctionSpace(mesht, 'P', 1)

R = Function(V)
q = Function(V)
A = pi*pow(R,2)

q_inlet = project(tt,qq, W)


U = as_vector([A,q])

F = as_vector([q, pow(q,2)/A + sqrt(A0*A)])

S = as_vector([0, 2*pi*R/deltab/Re*q/A])



bc_q = DirichletBC(V, ql, on_boundary and near(x[0], 0))
