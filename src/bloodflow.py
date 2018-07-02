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

L, Tf = 20.8, data_q[-1,0]
Nx, Nt = 10, len(data_q[:,0])

r0 = 0.37
E = 1.0e+6
H = 0.01
Re = 1.0
nu = 0.04
T = 1.0
db = np.sqrt(nu*T/2/pi)

f = 4*E*H/3/r0

mesh = RectangleMesh(Point(0,0), Point(L, Tf), Nx, Nt)

V = FiniteElement("CG", mesh.ufl_cell(), 1)
V2 = FunctionSpace(mesh, V*V)


U0 = Function(V2)
#A00, q0 = split(U0)
#q0 = U0.sub(1)

q0 = Constant(5.0)
q_inlet = q0
q_outlet = q0

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
bc_init = DirichletBC(V2.sub(0), Constant(A0), init_bdry)

#bcs = [bc_inlet, bc_outlet, bc_init]
bcs = [bc_inlet, bc_init]


U = Function(V2)
A, q = split(U)

#UR = Function(V2)
#R, q = split(UR)

v1, v2 = TestFunctions(V2)

e1 = Constant((1,0))
e2 = Constant((0,1))


FF = div(A*e2)*v1*dx\
  + div(q*(v1*e1+v2*e2))*dx\
  + div((pow(q,2)/A+f*sqrt(A0*A))*e1)*v2*dx\
  + 2*sqrt(pi)/db/Re*q/sqrt(A)*v2*dx


"""
FFR = 2*pi*R*dot(grad(R),e2)*v1*dx\
    + dot(grad(q),e1)*v1*dx\
    + dot(grad(q),e2)*v2*dx\
    + (2/pi*q/pow(R,2)*dot(grad(q),e1)-4/pi*pow(q,2)/pow(R,2)+f*sqrt(pi*A0)*dot(grad(R),e1))*v2*dx\
    + 2/db/Re*q/R*v2*dx
"""
solve(FF == 0, U, bcs)





