"""
Implementation of the 1D blood flow equations.
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

from fenics import *
import numpy as np
import matplotlib.pyplot as plt

data_q = np.genfromtxt('../data/example_inlet.csv', delimiter = ',')
#print(data_q)
#plt.plot(data_q[:,0], data_q[:,1])
#plt.savefig('data.png')

tt = data_q[:,0]
qq = data_q[:,1]

L, T = 20.8, data_q[-1,0]
Nx, Nt = 10, 10 # len(data_q[:,0])

r0 = 0.37
E = 1.0e+6
H = 0.01
Re = 1.0
nu = 0.046
db = np.sqrt(nu*T/2/pi)

f = 4*E*H/3/r0

mesh = IntervalMesh(0, L, Nx)

V = FiniteElement("CG", mesh.ufl_cell(), 1)
V2 = FunctionSpace(mesh, V*V)

q_inlet = Constant(5.0)
q0 = q_inlet

A0 = Constant(pi*pow(r0,2))


tol = 1.e-14

def inlet_bdry(x, on_boundary):
	return on_boundary and near(x[0],0,tol)
	
def outlet_bdry(x, on_boundary):
	return on_boundary and near(x[0],L,tol)
	
def init_bdry(x, on_boundary):
	return on_boundary and near(x[1],0,tol)

bc_inlet = DirichletBC(V2.sub(1), q_inlet, inlet_bdry)
bc_outlet = DirichletBC(V2.sub(1), q_outlet, outlet_bdry)

#bcs = [bc_inlet, bc_outlet]
bcs = [bc_inlet]


U = Function(V2)
A, q = split(U)

U_n = Constant((A0,q0))

v1, v2 = TestFunctions(V2)

for n in range(Nt):

	FF = A*v1*dx\
	   + q*v2*dx\
	   + dt*grad(q)*v1*dx\
	   + dt*(pow(q,2)/A+f*sqrt(A0*A)*v2*ds\
	   + dt*grad(pow(q,2)/A+f*sqrt(A0*A))*grad(v2)*dx\
	   + dt*2*sqrt(pi)/db/Re*q/sqrt(A)*v2*dx\
	   - dot(U_n,(v1,v2))*dx
#+ dt*grad(pow(q,2)/A+f*sqrt(A0*A))*v2*dx\
	solve(FF == 0, U, bcs)
	
	U_n.assign(U)







