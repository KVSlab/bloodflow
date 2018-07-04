"""
------------------------------------------------------
--- Implementation of the 1D blood flow equations. ---
--- Author : Syver Døving Agdestein                ---
------------------------------------------------------


------------------------------
- Definition of the problem: -
------------------------------

Variables: x and t:

	x in [0,L]
	t in [0,T]
	
	where L is the length of the artery and T is the duration of one cardiac cycle (period).


    t
    ^
    |
    |
  T +--------------------------------+
    |                                .
    |                                .
    |                                .
    |                                .
    |                                .
    |                                .
    |                                .
    |                                .
    |                                .
    |                                .
  0 +––––––––––––––––––––––––––––––––+–––> x
    0                                L


dU/dt + dF/dx = S

             | A(x,t) |                |            q            |                  |       0        |
U = U(x,t) = |        |,    F = F(U) = |                         |,      S = S(U) = |                |
             | q(x,t) |                | q²/A + f(r0)sqrt(A_0*A) |                  | -2πR/dbRe q/A  |

A = πR²

p(x,t) - p0 = f(r_0)(1 - sqrt(A_0/A))

r_0(x) = r_u*(r_d/r_u)^(x/L)


------------------------
- Boundary Conditions: -
------------------------


Left (inlet) border:

	q(0,t) = q_inlet(t), a given function


Bottom border:

	q(x,0) = q_inlet(0) since r_0 = constant, which makes the artery a perfect cylinder.
	R(x,0) = r_0(x) <=> A(x,0) = A_0(x)

Right (outlet) border:

	p = p_0  <=> A = A_0

     t
     ^
     |
     |
   T +--------------------------------+
     |                                '
     |                                '
     |                                '
     |                                '
q=q_0|                                ' p=p_0 (A = A_0)
     |                                '
     |                                '
     |                                '
     |                                '
     |                                '
   0 +––––––––––––––––––––––––––––––––+–––> x
     0          q=q_0, R=r_0          L


-----------------------------------------------------------
-----------------------------------------------------------
"""

from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as ip
from mpl_toolkits.mplot3d import Axes3D

# Pressure unit converting functions ('unit' = g cm-1 s-2)
def unit_to_mmHg(p):
	return 76/101325*p
	
def mmHg_to_unit(p):
	return 101325/76*p


# Import the inlet flow data
data_q = np.genfromtxt('../data/example_inlet.csv', delimiter = ',')
#plt.plot(data_q[:,0], data_q[:,1])
#plt.savefig('../output/r0/data.png')

ttt = data_q[:,0]
qqq = data_q[:,1]

L, T = 20.8, data_q[-1,0]
#Nx, Nt = 100, len(data_q[:,0])
Nx, Nt = 100, 100

xx = np.linspace(0,L,Nx)

qt = ip.interp1d(ttt, qqq)
tt = np.linspace(0,T,Nt)
qq = qt(tt)

#plt.plot(tt,qq)
#plt.savefig('../output/r0/interpolated_data.png')

#dt = ttt[1:]-ttt[:-1]
#dt = min(ttt[1:]-ttt[:-1])
dt = T/Nt

nu = 0.046
Re = 10.0/nu/1.0
db = np.sqrt(nu*T/2/pi)
p0 = mmHg_to_unit(90) # Unit: g cm-1 s-2

ru = 0.37
rd = 0.37
k = ln(rd/ru)/L

k1 = 2.0e7
k2 = -22.53
k3 = 8.65e5

Eh = ru*(k1*exp(k2*ru)+k3)

# Definition of mesh and function spaces
mesh = IntervalMesh(Nx, 0, L)

elV = FiniteElement("CG", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, elV)
V2 = FunctionSpace(mesh, elV*elV)

# Definition of trial function
U = Function(V2)
A, q = split(U)

# Definition of test functions
v1, v2 = TestFunctions(V2)

# Initial vessel-radius and deduced quantities, all functions of the spatial variable
r0 = Expression('ru*pow(rd/ru, x[0]/L)', degree = 2, ru = ru, rd = rd, L = L)
A0 = Expression('pi*pow(ru,2)*pow(rd/ru,2*x[0]/L)', degree = 2, ru = ru, rd = rd, L = L)
f = Expression('4/3*Eh/ru*pow(ru/rd,x[0]/L)', degree = 2, ru = ru, rd = rd, L = L, Eh = Eh)
dfdr = Expression('4/3*k1*k2*exp(k2*ru*pow(rd/ru,x[0]/L))', degree = 2, ru = ru, rd = rd, L = L, Eh = Eh, k1 = k1, k2 = k2)
drdx = Expression('log(rd/ru)/L*ru*pow(rd/ru,x[0]/L)', degree = 2, ru = ru, rd = rd, L = L)

# Inlet flow, defined at one single given time t_n (starting at t_0)
q_in = Function(V)
q_in.assign(Constant(qq[0]))

# Outlet area
#A_out = Constant(A0(L)/pow(1+p0/f(L),2))
A_out = Constant(A0(L))

# The initial value of the trial function is deduced from the bottom boundary conditions
U_n = Function(V2)
U_n.assign(Expression(('pi*pow(ru,2)*pow(rd/ru,2*x[0]/L)', 'q00'), degree = 2, ru = ru, rd = rd, L = L, q00 = qq[0]))


# Spatial boundary conditions
tol = 1.e-14

def inlet_bdry(x, on_boundary):
	return on_boundary and near(x[0],0,tol)
	
def outlet_bdry(x, on_boundary):
	return on_boundary and near(x[0],L,tol)

bc_outlet = DirichletBC(V2.sub(0), A_out, outlet_bdry)
bc_inlet = DirichletBC(V2.sub(1), q_in, inlet_bdry)

bcs = [bc_inlet, bc_outlet]

"""
FF = A*v1*dx\
   + q*v2*dx\
   + dt*grad(q)[0]*v1*dx\
   + dt*grad(pow(q,2)/(A+1.e-16)+4/3*Eh/ru*pow(ru/rd,x[0]/L)*sqrt(A0*(A+1.e-16)))[0]*v2*dx\
   + dt*2*sqrt(pi)/db/Re*q/sqrt(A+1.e-16)*v2*dx\
   - (2*sqrt(A)*(sqrt(pi)*4/3*Eh/ru*pow(ru/rd,x[0]/L)\
                +sqrt(A0)*4/3*k1*k2*exp(k2*ru*pow(rd/ru,x[0]/L)))\
     -A*4/3*k1*k2*exp(k2*ru*pow(rd/ru,x[0]/L))\
     )*ln(rd/ru)/L*ru*pow(rd/ru,x[0]/L)*v2*dx\
   - U_n[0]*v1*dx\
   - U_n[1]*v2*dx

# Variational form: FF == 0
FF = A*v1*dx\
   + q*v2*dx\
   + dt*grad(q)[0]*v1*dx\
   + dt*grad(pow(q,2)/(A+1.e-16)+f*sqrt(A0*(A+1.e-16)))[0]*v2*dx\
   + dt*2*sqrt(pi)/db/Re*q/sqrt(A+1.e-16)*v2*dx\
   - dt*(2*sqrt(A+1.e-16)*(sqrt(pi)*f+sqrt(A0)*dfdr)-(A+1.e-16)*dfdr)*drdx*v2*dx\
   - U_n[0]*v1*dx\
   - U_n[1]*v2*dx
"""

# Variational form: FF == 0
FF = A*v1*dx\
   + q*v2*dx\
   + dt*grad(q)[0]*v1*dx\
   + dt*(pow(q,2)/(A+1.e-16)+f*sqrt(A0*(A+1.e-16)))*v2*ds\
   - dt*(pow(q,2)/(A+1.e-16)+f*sqrt(A0*(A+1.e-16)))*grad(v2)[0]*dx\
   + dt*2*sqrt(pi)/db/Re*q/sqrt(A+1.e-16)*v2*dx\
   - dt*(2*sqrt(A+1.e-16)*(sqrt(pi)*f+sqrt(A0)*dfdr)-(A+1.e-16)*dfdr)*drdx*v2*dx\
   - U_n[0]*v1*dx\
   - U_n[1]*v2*dx


# Matrices for storing the solution
qmat = np.zeros([Nx, Nt])
Amat = np.zeros([Nx, Nt])
pmat = np.zeros([Nx, Nt])

qmat[:,0] = qq[0]*np.ones(Nx)
Amat[:,0] = [A0([x]) for x in xx]

#xdmffile_U = XDMFFile('../output/r0/bloodflow1Dr.xdmf')

# Progress bar
progress = Progress('Time-stepping')
set_log_level(PROGRESS)

t = 0

# Time-stepping
for n in range(Nt-1):
	
	print('Iteration '+str(n))
	
	t += dt
	
	# U_n+1 is solution of FF == 0
	solve(FF == 0, U, bcs)
	
	# Update previous solution
	U_n.assign(U)
	
	# Update inlet boundary condition
	q_in.assign(Constant(qq[n]))
	
	#xdmffile_U.write(U, dt)
	
	# Store solution at time t_n+1
	qmat[:,n+1] = [q([x]) for x in xx]
	Amat[:,n+1] = [A([x]) for x in xx]
	
	# Update progress bar
	progress.update((t+dt)/T)
	

X, Y = np.meshgrid(tt, xx)

# Assembly of the pressure matrix
for n in range(Nt):
	for i in range(Nx):
		pmat[i,n] = unit_to_mmHg(p0 + f(xx[i])*(1-np.sqrt(A0([xx[i]])/Amat[i,n])))


# Area plot
fig = plt.figure(figsize=(8,6))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Amat, rstride=1, cstride=1,  cmap='viridis', linewidth=0, antialiased=False)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('A')
ax.set_ylim(min(xx), max(xx))
ax.set_xlim(min(tt), max(tt))
plt.savefig('../output/r0/arear.png')


# Flow plot
fig = plt.figure(figsize=(8,6))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, qmat, rstride=1, cstride=1,  cmap='viridis', linewidth=0, antialiased=False)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('q')
ax.set_ylim(min(xx), max(xx))
ax.set_xlim(min(tt), max(tt))
#ax.set_zlim(-15,0.0)
plt.savefig('../output/r0/flowr.png')


# Pressure plot
fig = plt.figure(figsize=(8,6))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, pmat, rstride=1, cstride=1,  cmap='viridis', linewidth=0, antialiased=False)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('p')
ax.set_ylim(min(xx), max(xx))
ax.set_xlim(min(tt), max(tt))
plt.savefig('../output/r0/pressurer.png')


