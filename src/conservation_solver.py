import sys
import numpy as np
from fenics import *

sys.path.insert(0, '../src')
import artery_network

class Conservation_Solver(object):

	def F_from_equation(a, U, x):
		"""Compute the flux term.
		:param U: Value of the solution
		:param x: Point of evaluation
		:return: F(x)
		U should be the value of the solution at the point x.
		"""
		return np.array([U[1], U[1]**2 + a.f(x)*np.sqrt(a.A0(x)*U[0])])

	def S_from_equation(a, U, x):
		"""Compute the source term.
		:param U: Value of the solution
		:param x: Point
		:return: S(x)
		U should be the value of the solution at the point x.
		"""
		S1 = 0
		S2 = -2*np.sqrt(np.pi)/a.db/a.Re*U[1]/np.sqrt(U[0])\
			+(2*np.sqrt(U[0])*(np.sqrt(np.pi)*a.f(x)\
							  +np.sqrt(a.A0(x))*a.dfdr(x))\
				-U[0]*a.dfdr(x))*a.drdx(x)
		return np.array([S1, S2])

	def compute_A_out(a, U_n, k_max=100, tol=1.0e-7):
		"""Compute the outlet boundary condition.
		:param U_n: Solution (function) at time step t_n
		:param k_max: Maximum number of iterations in Piccards scheme
		:param tol: Tolerance for Piccards fixed point iteration scheme
		:return: Outlet boundary value of A at time step t_(n+1).
		"""
		# Spatial step, scaled to satisfy the CFL condition
		deltax = 10*a.L/a.Nx
		x2, x1, x0 = a.L-2*deltax, a.L-deltax, a.L
		x21, x10 = a.L-1.5*deltax, a.L-0.5*deltax

		Um2, Um1, Um0 = U_n(x2), U_n(x1), U_n(x0)

		# Values at time step n
		Fm2, Sm2 = F_from_equation(a, Um2, x2), S_from_equation(a, Um2, x2)
		Fm1, Sm1 = F_from_equation(a, Um1, x1), S_from_equation(a, Um1, x1)
		Fm0, Sm0 = F_from_equation(a, Um0, x0), S_from_equation(a, Um0, x0)

		# Values at time step n+1/2
		U_half_21 = (Um1+Um2)/2 - a.dt/deltax*(Fm1-Fm2) + a.dt/4*(Sm1+Sm2)
		U_half_10 = (Um0+Um1)/2 - a.dt/deltax*(Fm0-Fm1) + a.dt/4*(Sm0+Sm1)
		F_half_21 = F_from_equation(a, U_half_21, x21)
		S_half_21 = S_from_equation(a, U_half_21, x21)
		F_half_10 = F_from_equation(a, U_half_10, x10)
		S_half_10 = S_from_equation(a, U_half_10, x10)

		# Value at time step n+1
		qm1 = Um1[1]\
			- a.dt/deltax*(F_half_10[1]-F_half_21[1])\
			+ a.dt/2*(S_half_10[1]+S_half_21[1])

		# Fixed point iteration
		pn = a.outlet_pressure(Um0[0])
		p = pn
		for k in range(k_max):
			p_old = p
			qm0 = Um0[1]\
				+ (p-pn)/a.R1\
				+ a.dt/a.R1/a.R2/a.CT*pn\
				- a.dt*(a.R1+a.R2)/a.R1/a.R2/a.CT*Um0[1]
			Am0 = Um0[0] - a.dt/deltax*(qm0-qm1)
			p = a.outlet_pressure(Am0)
			if abs(p-p_old) < tol:
				break

		return Am0


	def solve(a, q_ins):
		"""Compute and store the solution to dU/dt + dF/dx = S.
		:param artery a: Artery on which the solution is to be computed
		:param q_ins: Vector containing inlet flow
		"""
		# Define trial function
		U = Function(a.V2)
		A, q = split(U)

		# Define test functions
		v1, v2 = TestFunctions(a.V2)

		# Inlet flow
		q_in = Function(a.V)
		q_in.assign(Constant(q_ins[0]))

		# Outlet area
		A_out = Function(a.V)
		A_out.assign(Constant(a.A0(a.L)))

		# Initial value deduced from bottom boundary conditions
		U_n = Function(a.V2)
		U_n.assign(Expression(('pi*pow(Ru, 2)*pow(Rd/Ru, 2*x[0]/L)', 'q00'),
			degree=2, Ru=a.Ru, Rd=a.Rd, L=a.L, q00=q_ins[0]))

		# Spatial boundary conditions
		tol = 1.e-14
		def inlet_bdry(x, on_boundary):
			return on_boundary and near(x[0], 0, tol)
		def outlet_bdry(x, on_boundary):
			return on_boundary and near(x[0], a.L, tol)
		bc_outlet = DirichletBC(a.V2.sub(0), A_out, outlet_bdry)
		bc_inlet = DirichletBC(a.V2.sub(1), q_in, inlet_bdry)
		bcs = [bc_inlet, bc_outlet]

		# Variational form
		FF = A*v1*dx\
		   + q*v2*dx\
		   + a.dt*grad(q)[0]*v1*dx\
		   + a.dt*(pow(q, 2)/(A+DOLFIN_EPS)\
				  +a.f*sqrt(a.A0*(A+DOLFIN_EPS)))*v2*ds\
		   - a.dt*(pow(q, 2)/(A+DOLFIN_EPS)\
				  +a.f*sqrt(a.A0*(A+DOLFIN_EPS)))*grad(v2)[0]*dx\
		   + a.dt*2*sqrt(pi)/a.db/a.Re*q/sqrt(A+DOLFIN_EPS)*v2*dx\
		   - a.dt*(2*sqrt(A+DOLFIN_EPS)*(sqrt(pi)*a.f+sqrt(a.A0)*a.dfdr)\
				  -(A+DOLFIN_EPS)*a.dfdr)*a.drdx*v2*dx\
		   - U_n[0]*v1*dx\
		   - U_n[1]*v2*dx

		# Progress bar
		progress = Progress('Time-stepping')
		set_log_level(PROGRESS)

		# Initialise time
		t = 0

		# Cardiac cycle iteration
		for n_cycle in range(N_cycles):
			
			# Store solution at multiples of time t_Nt (beginning of cycle)
			a.solution[0].assign(U_n)
			
			# Time-stepping for one period
			for n in range(Nt-1):

				print('Iteration '+str(n))

				t += a.dt

				# U_n+1 is solution of FF == 0
				solve(FF == 0, U, bcs)

				# Update previous solution
				U_n.assign(U)

				# Update inlet boundary condition
				q_in.assign(Constant(q_ins[n]))

				# Update outlet boundary condition
				A_out_value = compute_A_out(U_n)
				A_out.assign(Constant(A_out_value))
				
				# Store solution at time t_(n+1)
				a.solution[n+1].assign(U)

				# Update progress bar
				progress.update((t+dt)/N_cycles/T)
