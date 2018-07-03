"""
Implementation of the 1D blood flow equations.
Author : Syver Døving Agdestein
--------------------------------------

dU/dt + dF/dx = S

             | A(x,t) |                |            q            |                  |       0        |
U = U(x,t) = |        |,    F = F(U) = |                         |,      S = S(U) = |                |
             | q(x,t) |                | q²/A + f(r0)sqrt(A_0*A) |                  | -2πR/dbRe q/A  |

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
import scipy.interpolate as ip
from mpl_toolkits.mplot3d import Axes3D

data_q = np.genfromtxt('../data/example_inlet.csv', delimiter = ',')
#print(data_q)
#plt.plot(data_q[:,0], data_q[:,1])
#plt.savefig('data.png')

ttt = data_q[:,0]
qqq = data_q[:,1]

L, T = 20.8, data_q[-1,0]
Nx, Nt = 100, 100 # len(data_q[:,0])

xx = np.linspace(0,L,Nx)

qt = ip.interp1d(ttt, qqq)
tt = np.linspace(0,T,Nt)
qq = qt(tt)
#qq = np.linspace(0,25,Nt)
#qq = 5.0*np.ones(Nt)

#plt.plot(tt,qq)
#plt.savefig('ttqq.png')

#dt = ttt[1:]-ttt[:-1]
#dt = min(ttt[1:]-ttt[:-1])
dt = T/Nt

r0 = 0.37
E = 1.0e+6
H = 0.01
nu = 0.046
Re = 10.0/nu/1.0
db = np.sqrt(nu*T/2/pi)

f = 4*E*H/3/r0

mesh = IntervalMesh(Nx, 0, L)

elV = FiniteElement("CG", mesh.ufl_cell(), 2)
V = FunctionSpace(mesh, elV)
V2 = FunctionSpace(mesh, elV*elV)

#q_inlet = Constant(5.0)
#q0 = Constant(qq[0])
q0 = Function(V)
q0.assign(Constant(qq[0]))
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
#bcs = [bc_inlet_q, bc_outlet_q]
#bcs = [bc_inlet_A, bc_inlet_q]
bcs = [bc_inlet_q]

U = Function(V2)
A, q = split(U)

v1, v2 = TestFunctions(V2)

U_n = Function(V2)
U_n.assign(Constant((A0,q0(0))))

xdmffile_U = XDMFFile('bloodflow1D.xdmf')
#xdmffile_A = XDMFFile('bloodflow1D_A.xdmf')
#xdmffile_q = XDMFFile('bloodflow1D_q.xdmf')

FF = A*v1*dx\
   + q*v2*dx\
   + dt*grad(q)[0]*v1*dx\
   + dt*grad(pow(q,2)/(A+1.e-14)+f*sqrt(A0*(A+1.e-14)))[0]*v2*dx\
   + dt*2*sqrt(pi)/db/Re*q/sqrt(A+1.e-14)*v2*dx\
   - U_n[0]*v1*dx\
   - U_n[1]*v2*dx

qmat = np.zeros([Nx, Nt])

qmat[:,0] = qq[0]*np.ones(Nx)

for i in range(3):
	
	t = 0

	for n in range(Nt-1):
		
		print('Iteration '+str(n))
		
		t += dt

		solve(FF == 0, U, bcs)
		#if n % (int(Nt/100)) == 0:
		#	plot(A)
			
		U_n.assign(U)
		
		q0.assign(Constant(qq[n]))
		
		xdmffile_U.write(U, dt)
		#xdmffile_A.write(A, dt)
		#xdmffile_q.write(q, dt)
		
		
		qmat[:,n+1] = np.array([q([xx[i]]) for i in range(Nx)])
		
		

plt.imshow(qmat)

#plt.scatter(0,qq[0])
#plt.ylim(-5, 30)

X, Y = np.meshgrid(tt, xx)

fig = plt.figure(figsize=(12,8))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, qmat, cmap='viridis', linewidth=0, antialiased=False)

plt.savefig('q.png')



