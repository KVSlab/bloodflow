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



dU/dt + dF/dx = S

             | A(x,t) |                |            q            |                  |       0        |
U = U(x,t) = |        |,    F = F(U) = |                         |,      S = S(U) = |                |
             | q(x,t) |                | q²/A + f(r0)sqrt(A_0*A) |                  | -2πR/dbRe q/A  |

A = πR²

p(x,t) - p0 = f(r_0)(1 - sqrt(A_0/A))

r_0 = constant


------------------------
- Boundary Conditions: -
------------------------


Left (inlet) border:

	q(0,t) = q_inlet(t), a given function


Bottom border:

	q(x,0) = q_inlet(0) since r_0 = constant, which makes the artery a perfect cylinder.
	R(x,0) = r_0 <=> A(x,0) = A_0

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

# Libraries
import numpy as np

from fenics import *

import matplotlib.pyplot as plt
import scipy.interpolate as ip
from mpl_toolkits.mplot3d import Axes3D


# Pressure unit converting functions ('unit' = g cm-1 s-2)
def unit_to_mmHg(p):
	return 76/101325*p
	
def mmHg_to_unit(p):
	return 101325/76*p


# Import the inlet flow data
data_q = np.genfromtxt('../data/example_inlet.csv', delimiter=',')
#plt.plot(data_q[:,0], data_q[:,1])
#plt.savefig('../output/constant_r0/data.png')

ttt = data_q[:, 0]
qqq = data_q[:, 1]

L, T = 20.8, data_q[-1, 0]
#Nx, Nt = 100, len(data_q[:,0])

# Nx/Nt ratio should be 1/3 (based on experience)
Nx, Nt = 100, 300

xx = np.linspace(0, L, Nx)

qt = ip.interp1d(ttt, qqq, kind='cubic')
tt = np.linspace(0, T, Nt)
qq = qt(tt)

#plt.plot(tt,qq)
#plt.savefig('../output/constant_r0/interpolated_data.png')

#dt = ttt[1:]-ttt[:-1]
#dt = min(ttt[1:]-ttt[:-1])
dt = T/Nt

nu = 0.046
Re = 10.0/nu/1.0
db = np.sqrt(nu*T/2/pi)
p0 = mmHg_to_unit(90)  # Unit: g cm-1 s-2

r0 = 0.37

k1 = 2.0e7
k2 = -22.53
k3 = 8.65e5
f = 4/3*(k1*exp(k2*r0)+k3)

mesh = IntervalMesh(Nx, 0, L)

elV = FiniteElement("CG", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, elV)
V2 = FunctionSpace(mesh, elV*elV)

# Initial area
A0 = Constant(pi*pow(r0, 2))



R1 = 25300
R2 = 13900
CT = 1.3384e-6

def F_from_equation(U):
	return np.array([U[1], U[1]**2 + f*np.sqrt(A0(L)*U[0])])

def S_from_equation(U):
	return np.array([0, -2*np.sqrt(np.pi)/db/Re*U[1]/np.sqrt(U[0])])

# Computes the outlet pressure at time t_n+1 from the values of the solution at the three end points m-2, m-1 and m (m=Nx-1) at time t_n.
# q_(m-1)^(n+1) is computed using Richtmyer's two step Lax-Wendroff method.
# q_m^(n+1) is computed using the Windkessel model, based on an initial estimate of p_m^(n+1) (starting at p_m^n).
def outlet_area(U_n, k_max=100, tol=1.0e-7):
	
	# Spatial step, many times larger than the one used in the finite elements scheme (to ensure convergence).
	deltax = 10*L/Nx
	
	Um2, Um1, Um0 = U_n(L-2*deltax), U_n(L-deltax), U_n(L)
	
	# Values at time step n
	Fm2, Sm2 = F_from_equation(Um2), S_from_equation(Um2)
	Fm1, Sm1 = F_from_equation(Um1), S_from_equation(Um1)
	Fm0, Sm0 = F_from_equation(Um0), S_from_equation(Um0)
	
	# Values at time step n+1/2
	U_half_21 = (Um1+Um2)/2 - dt/deltax*(Fm1-Fm2) + dt/4*(Sm1+Sm2)
	U_half_10 = (Um0+Um1)/2 - dt/deltax*(Fm0-Fm1) + dt/4*(Sm0+Sm1)
	F_half_21, S_half_21 = F_from_equation(U_half_21), S_from_equation(U_half_21)
	F_half_10, S_half_10 = F_from_equation(U_half_10), S_from_equation(U_half_10)
	
	# Value at time step n+1
	qm1 = Um1[1] - dt/deltax*(F_half_10[1]-F_half_21[1]) + dt/2*(S_half_10[1]+S_half_21[1])
	
	# Fixed point iteration
	pn = p0 + f*(1-np.sqrt(A0(L)/Um0[0]))
	p = pn
	for k in range(k_max):
		p_old = p
		qm0 = Um0[1] + (p-pn)/R1 + dt/R1/R2/CT*pn - dt*(R1+R2)/R1/R2/CT*Um0[1]
		Am0 = Um0[0] - dt/deltax*(qm0-qm1)
		p = p0 + f*(1-np.sqrt(A0(L)/Am0))
		if abs(p-p_old) < tol:
			break
	
	return Am0


# Definition of trial functions
U = Function(V2)
A, q = split(U)

# Definition of test functions
v1, v2 = TestFunctions(V2)

# Inlet flow at a given time t_n (initially t_0)
q_in = Function(V)
q_in.assign(Constant(qq[0]))

# Outlet area
A_out = Function(V)
A_out.assign(Constant(A0(L)))

# The initial value of the trial function is deduced from the bottom boundary conditions
U_n = Function(V2)
U_n.assign(Constant((A0, q_in(0))))

# Spatial boundary conditions
tol = 1.e-14

def inlet_bdry(x, on_boundary):
	return on_boundary and near(x[0], 0, tol)
	
def outlet_bdry(x, on_boundary):
	return on_boundary and near(x[0], L, tol)

bc_outlet = DirichletBC(V2.sub(0), A_out, outlet_bdry)
bc_inlet = DirichletBC(V2.sub(1), q_in, inlet_bdry)

bcs = [bc_inlet, bc_outlet]


# Variational form: FF == 0
# In order for FEniCS to handle the expression correctly, a small perturbation of A seems necessary. Since A > 0.4, this should not numerically change the solution.
FF = A*v1*dx\
   + q*v2*dx\
   + dt*q*v1*ds\
   + dt*(pow(q, 2)/(A+DOLFIN_EPS)+f*sqrt(A0*(A+DOLFIN_EPS)))*v2*ds\
   - dt*q*grad(v1)[0]*dx\
   - dt*(pow(q, 2)/(A+DOLFIN_EPS)+f*sqrt(A0*(A+DOLFIN_EPS)))*grad(v2)[0]*dx\
   + dt*2*sqrt(pi)/db/Re*q/sqrt(A+DOLFIN_EPS)*v2*dx\
   - U_n[0]*v1*dx\
   - U_n[1]*v2*dx



# Number of cardiac cycles
N_cycles = 1

# File for storing the solution
xdmffile_A = XDMFFile(mpi_comm_world(), '../output/constant_r0/area.xdmf')
xdmffile_q = XDMFFile(mpi_comm_world(), '../output/constant_r0/flow.xdmf')

# Initial values (t_0)
xdmffile_A.write_checkpoint(U_n.split()[0], 'area', 0)
xdmffile_q.write_checkpoint(U_n.split()[1], 'flow', 0)

# Progress bar
progress = Progress('Time-stepping')
set_log_level(PROGRESS)

t = 0

# Time-stepping
for n_cycle in range(N_cycles):

	for n in range(0, Nt-1):
		
		print('Iteration '+str(n))
		
		t += dt
		
		# Update inlet boundary condition
		q_in.assign(Constant(qq[n+1]))
		
		# Update outlet boundary condition
		A_out_value = outlet_area(U_n)
		A_out.assign(Constant(A_out_value))

		# U is solution of FF == 0
		solve(FF == 0, U, bcs)
		
		# Update previous solution
		U_n.assign(U)
		
		# Update XDMF-files
		xdmffile_A.write_checkpoint(U.split()[0], 'area', t)
		xdmffile_q.write_checkpoint(U.split()[1], 'flow', t)
		
		# Update progress bar
		progress.update((n_cycle*T+t+dt)/N_cycles/T)
		
	U_n.assign(U)

xdmffile_A.close()
xdmffile_q.close()

########################## assert


