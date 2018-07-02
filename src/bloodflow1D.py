"""
Implementation of the 1D blood flow equations.
Author : Syver Døving Agdestein
--------------------------------------

dU/dt + dF/dx = S

             | A(x,t) |                |       q            |                  |       0        |
U = U(x,t) = |        |,    F = F(U) = |                    |,      S = S(U) = |                |
             | q(x,t) |                | q²/A + sqrt(A_0*A) |                  | -2πR/dbRe q/A  |

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

dt = T/Nt

r0 = 0.37
E = 1.0e+6
H = 0.01
Re = 1.0
nu = 0.046
db = np.sqrt(nu*T/2/pi)

f = 4*E*H/3/r0

mesh = IntervalMesh(Nx, 0, L)

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


bc_inlet_A = DirichletBC(V2.sub(0), A0, inlet_bdry)
bc_outlet_A = DirichletBC(V2.sub(0), A0, outlet_bdry)
bc_inlet_q = DirichletBC(V2.sub(1), q0, inlet_bdry)
bc_outlet_q = DirichletBC(V2.sub(1), q0, outlet_bdry)

#bcs = [bc_inlet_A, bc_outlet_A, bc_inlet_q, bc_outlet_q]
#bcs = [bc_inlet, bc_outlet]
bcs = [bc_inlet_q]


U = Function(V2)
A, q = split(U)

v1, v2 = TestFunctions(V2)

U_n = Function(V2)
U_n.assign(Constant((A0,q0)))

xdmffile_U = XDMFFile('bloodflow1D.xdmf')
#xdmffile_A = XDMFFile('bloodflow1D_A.xdmf')
#xdmffile_q = XDMFFile('bloodflow1D_q.xdmf')

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

	xdmffile_U.write(U, n*T/Nt)
	#xdmffile_A.write(A, n*T/Nt)
	#xdmffile_q.write(q, n*T/Nt)
	



