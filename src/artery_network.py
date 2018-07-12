__author__ = 'Syver Døving Agdestein'


import numpy as np
import numpy.linalg as npl
from fenics import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from artery import Artery


class Artery_Network(object):
	"""Describe a network of arteries.
	:param order: Number of arterial levels
	:param Ru: Upstream radii
	:param Rd: Downstream radii
	:param L: Vessel lengths
	:param k: Vectors containing k1, k2 and k3 from the relation Eh/R0
	:param Re: Reynolds' number
	:param p0: Diastolic pressure
	"""
	def __init__(self, order, Ru, Rd, L, k1, k2, k3, Re, p0):
		self.order = order
		self.arteries = [0] * (2**order-1)
		for i in range(len(arteries)):
			arteries[i] = Artery(Ru[i], Rd[i], L[i], k1[i],
								 k2[i], k3[i], Re, p0)


	def daughter_vessels(self, i):
		"""Find daughter vessels.
		:param i: Index of parent vessel
		:return: Indices of daughter vessels
		"""
		if i < 0 or i >= 2**(self.order-1):
			raise Exception('Vessel index out of range')
			
		return 2*i+1, 2*i+2


	def parent_vessel(self, i):
		"""Find parent vessel.
		:param i: Index of daughter vessel
		:return: Index of parent vessel
		"""
		if i <= 0 or i >= 2**self.order:
			raise Exception('Vessel index out of range')
			
		# d1 is odd, d2=d1+1 is pair
		return int((i-1)/2)
		
	def sister_vessel(self, i):
		"""Find sister vessel.
		:param i: Index of vessel
		:return: Index of sister vessel
		"""
		if i%2 == 0:
			return i-1
		else:
			return i+1


	def flux(self, a, U, x):
		"""Compute the flux term.
		:param a: Artery on which the flux term is computed
		:param U: Value of solution
		:param x: Point of evaluation
		:return: F(x)
		U should be the value of the solution at the point x.
		"""
		return np.array([U[1], U[1]**2 + a.f(x)*np.sqrt(a.A0(x)*U[0])])


	def source(self, a, U, x):
		"""Compute the source term.
		:param a: Artery on which the source term is computed
		:param U: Value of solution
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
	
	def compute_U_half(self, a, U0, U1, x0, x1):
		"""Compute a half step.
		:param a: Current artery
		:param U0: Left point
		:param U1: Right point
		:param x0: x-value at which U0 is the solution
		:param x1: x-value at which U1 is the solution
		:return: Middle point, half a time step later
		"""
		# Value of terms at time t_n
		F0, S0 = self.flux(U0, x0), self.source(U0, x0)
		F1, S1 = self.flux(U1, x1), self.source(U1, x1)
		
		return (Um1+Um0)/2 - a.dt/(x1-x0)*(F1-F0) + a.dt/4*(Sm1+Sm0) 
		

	def compute_A_out(self, a, k_max=100, tol=1.0e-7):
		"""Compute the outlet boundary condition.
		:param k_max: Maximum number of iterations in Piccards scheme
		:param tol: Tolerance for Piccards fixed point iteration scheme
		:return: Outlet boundary value of A at time step t_(n+1).
		"""
		# Spatial step, scaled to satisfy the CFL condition
		deltax = 10*a.dex
		x2, x1, x0 = a.L-2*deltax, a.L-deltax, a.L
		x21, x10 = a.L-1.5*deltax, a.L-0.5*deltax
		Um2, Um1, Um0 = U_n(x2), U_n(x1), U_n(x0)

		# Values at time step n
		Fm2, Sm2 = self.flux(Um2, x2), self.source(Um2, x2)
		Fm1, Sm1 = self.flux(Um1, x1), self.source(Um1, x1)
		Fm0, Sm0 = self.flux(Um0, x0), self.source(Um0, x0)

		# Values at time step n+1/2
		U_half_21 = self.compute_U_half(a, Um2, Um1, x2, x1)
		U_half_10 = self.compute_U_half(a, Um1, Um0, x1, x0)
		F_half_21 = self.flux(U_half_21, x21)
		S_half_21 = self.source(U_half_21, x21)
		F_half_10 = self.flux(U_half_10, x10)
		S_half_10 = self.source(U_half_10, x10)

		# Value at time step n+1
		qm1 = Um1[1]\
			- self.dt/deltax*(F_half_10[1]-F_half_21[1])\
			+ self.dt/2*(S_half_10[1]+S_half_21[1])

		# Fixed point iteration
		pn = self.outlet_pressure(Um0[0])
		p = pn
		for k in range(k_max):
			p_old = p
			qm0 = Um0[1]\
				+ (p-pn)/self.R1\
				+ self.dt/self.R1/self.R2/self.CT*pn\
				- self.dt*(self.R1+self.R2)/self.R1/self.R2/self.CT*Um0[1]
			Am0 = Um0[0] - self.dt/deltax*(qm0-qm1)
			p = self.outlet_pressure(Am0)
			if abs(p-p_old) < tol:
				break

		return Am0


	def problem_function(self, p, d1, d2, x):
		"""Compute the function representing the system of equations <system>.
		If x is the solution to the problem, then function(x) = 0.
		:param p: Parent artery
		:param d1: First daughter vessel
		:param d2: Second daughter vessel
		:param x: Current point, an 18-dimensional vector
		:return: Value of function in x
		"""
		# Abbreviations
		A0p, A01, A02 = p.A0(p.L), d1.A0(0), d2.A0(0)
		fp, f1, f2 = p.f(p.l),  d1.f(0), d2.f(0)
		dbp, db1, db2 = p.db, d1.db, d2.db
		Rep, Re1, Re2 = p.Re, d1.Re, d2.Re
		dfdrp, dfdr1, dfdr2 = p.dfdr(p.L), d1.dfdr(0), d2.dfdr(0)
		drdxp, drdx1, drdx2 = p.drdx(p.L), d1.drdx(0), d2.drdx(0)
		rpi = sqrt(np.pi)
		
		# Ghost half terms
		Fp = self.flux(p, np.array([x[11], x[2]), p.L+p.dex/2)
		F1 = self.flux(d1, np.array([x[14], x[5]), -d1.dex/2)
		F2 = self.flux(d2, np.array([x[17], x[8]), -d2.dex/2)
		Sp = self.source(p, np.array([x[11], x[2]), p.L+p.dex/2)
		S1 = self.source(d1, np.array([x[14], x[5]), -d1.dex/2)
		S2 = self.source(d2, np.array([x[17], x[8]), -d2.dex/2)
		
		# Compute half-time-step-values in M-1/2 for p and 1/2 for d1 and d2
		Um1p, Um0p = p.U_n(p.L-p.dex), U_n(p.L)
		U0d1, U1d1 = d1.U_n(0), d1.U_n(d1.dex)
		U0d2, U1d2 = d2.U_n(0), d2.U_n(d2.dex)

		U_half_p = compute_U_half(Um1p, Um0p, p.L-p.dex, p.L)
		U_half_1 = compute_U_half(U0d1, U1d1, 0, d1.dex)
		U_half_2 = compute_U_half(U0d2, U1d2, 0, d2.dex)
		
		F_half_p = self.flux(U_half_p, L-p.dex/2)
		S_half_p = self.source(U_half_p, L-p.dex/2)
		F_half_1 = self.flux(U_half_1, d1.dex/2)
		S_half_1 = self.source(U_half_1, d1.dex/2)
		F_half_2 = self.flux(U_half_1, d1.dex/2)
		S_half_2 = self.source(U_half_1, d1.dex/2)
		
		# Function value
		y = np.zeros(18)
		
		# Entries from equation (20)
		y[0] = 2*x[1] - U_half_p[1] - x[2]
		y[1] = 2*x[4] - x[5] - U_half_1[1]
		y[2] = 2*x[7] - x[8] - U_half_2[1]
		
		# Entries from equation (21)
		y[3] = 2*x[10] - U_half_p[0] - x[11]
		y[4] = 2*x[13] - x[14] - U_half_1[0]
		y[5] = 2*x[16] - x[17] - U_half_2[0]
		
		# Entries from equation (22)
		y[6] = x[1] - x[4] - x[7]
		y[7] = x[2] - x[5] - x[8]
		
		# Entries from equation (23)
		y[8] = fp*(1-np.sqrt(A0p/x[10])) - f1*(1-np.sqrt(A01/x[13]))
		y[9] = fp*(1-np.sqrt(A0p/x[10])) - f2*(1-np.sqrt(A02/x[16]))
		y[10] = fp*(1-np.sqrt(A0p/x[9])) - f1*(1-np.sqrt(A01/x[12]))
		y[11] = fp*(1-np.sqrt(A0p/x[9])) - f2*(1-np.sqrt(A02/x[15]))
		
		# Entries from equation (26) 
		y[12] = x[0] - Um0p[1] + p.dt/p.dex(Fp[1] - F_half_p[1])\
			  - p.dt/2*(Sp[1]+S_half_p[1])
		y[13] = x[3] - U0d1[1] + d1.dt/d1.dex(F1[1] - F_half_1[1])\
			  - d1.dt/2*(S2[1]+S_half_2[1])
		y[14] = x[6] - U0d2[1] + d2.dt/d2.dex(F2[1] - F_half_2[1])\
			  - d2.dt/2*(S2[1]+S_half_2[1])
		
		# Entries from equation (27)
		y[15] = x[9] - Um0p[0] + p.dt/p.dex(Fp[0] - F_half_p[0])
		y[16] = x[12] - U0d1[0] + d1.dt/d1.dex(F1[0] - F_half_1[0])
		y[17] = x[15] - U0d2[0] + d2.dt/d2.dex(F2[0] - F_half_2[0])
		
		return y


	def jacobian(self, p, d1, d2, x):
		"""Compute the jacobian matrix for the system of equations <system>.
		:param p: Parent artery
		:param d1: First daughter vessel
		:param d2: Second daughter vessel
		:param x: Current point, an 18-dimensional vector
		:return: Jacobian matrix
		"""
		# Abbreviations
		A0p, A01, A02 = p.A0(p.L), d1.A0(0), d2.A0(0)
		fp, f1, f2 = p.f(p.l),  d1.f(0), d2.f(0)
		dbp, db1, db2 = p.db, d1.db, d2.db
		Rep, Re1, Re2 = p.Re, d1.Re, d2.Re
		dfdrp, dfdr1, dfdr2 = p.dfdr(p.L), d1.dfdr(0), d2.dfdr(0)
		drdxp, drdx1, drdx2 = p.drdx(p.L), d1.drdx(0), d2.drdx(0)
		rpi = sqrt(np.pi)
		
		# Jacobian matrix
		J = np.zeros([18, 18])
		
		# Entries from equation (20)
		J[1, 0] = 2
		J[2, 0] = -1
		J[4, 1] = 2
		J[5, 1] = -1
		J[7, 2] = 2
		J[8, 2] = -1
		
		# Entries from equation (21)
		J[10, 3] = 2
		J[11, 3] = -1
		J[13, 4] = 2
		J[14, 4] = -1
		J[16, 5] = 2
		J[17, 5] = -1
		
		# Entries from equation (22)
		J[1, 6] = 1
		J[4, 6] = -1
		J[7, 6] = -1
		J[2, 7] = 1
		J[5, 7] = -1
		J[8, 7] = -1
		
		# Entries from equation (23)
		J[10, 8] = fp*np.sqrt(A0p)/2/pow(x[10], 3/2)
		J[13, 8] = -f1*np.sqrt(A01)/2/pow(x[13], 3/2)
		J[10, 9] = fp*np.sqrt(A0p)/2/pow(x[10], 3/2)
		J[16, 9] = -f2*np.sqrt(AO2)/2/pow(x[16], 3/2)
		J[9, 10] = fp*np.sqrt(A0p)/2/pow(x[9], 3/2)
		J[12, 10] = -f2*np.sqrt(A01)/2/pow(x[12], 3/2)
		J[9, 11] = fp*np.sqrt(A0p)/2/pow(x[9], 3/2)
		J[15, 11] = -f2*np.sqrt(A02)/2/pow(x[15], 3/2)
		
		# Entries from equation (26)
		J[0, 12] = 1
		J[2, 12] = p.dt/p.dex*2*x[2]/x[11] + p.dt*rpi/dbp/Rep/np.sqrt(x[11])
		J[11, 12] = p.dt/p.dex*(-pow(x[2]/x[11], 2) + fp/2*np.sqrt(A0p/x[11]))\
				  - p.dt/2*(rpi/dbp/Rep*x[2]/pow(x[11], 3/2)\
						   + (1/np.sqrt(x[11])*(rpi*fp+np.sqrt(A0p)*dfdrp)\
											   - dfdrp)*drdxp)
		J[3, 13] = 1
		J[5, 13] = -d1.dt/d1.dex*2*x[5]/x[14] + d1.dt*rpi/db1/Re1/np.sqrt(x[14])
		J[14, 13] = d1.dt/d1.dex*(pow(x[5]/x[14], 2) - f1/2*np.sqrt(A01/x[14]))\
				  - d1.dt/2*(rpi/db1/Re1*x[5]/pow(x[14], 3/2)\
							+ (1/np.sqrt(x[14])*(rpi*f1+np.sqrt(A01)*dfdr1)\
												- dfdr1)*drdx1)
		J[6, 14] = 1
		J[8, 14] = -d2.dt/d2.dex*2*x[8]/x[17] + d2.dt*rpi/db2/Re2/np.sqrt(x[17])
		J[17, 14] = d2.dt/d2.dex*(pow(x[8]/x[17], 2) - f2/2*np.sqrt(A02/x[17]))\
				  - d2.dt/2*(rpi/db2/Re2*x[8]/pow(x[17], 3/2)\
							+ (1/np.sqrt(x[17])*(rpi*f2+np.sqrt(A02)*dfdr2)\
												- dfdr2)*drdx2)
		
		# Entries from equation (27)
		J[2, 15] = p.dt/p.dex
		J[11, 15] = 1
		J[5, 16] = d1.dt/d1.dex
		J[14, 16] = 1
		J[8, 17] = d1.dt/d1.dex
		J[17, 17] = 1
		
		return J
		

	def newton(self, p, d1, d2, x=np.zeros(18), k_max=1000, tol=1.e-10):
		"""Compute solution to the system of equations.
		:param p: Parent artery
		:param d1: First daughter vessel
		:param d2: Second daughter vessel
		:param x: Current solution, 18-dimensional vector
		:param k_max: Maximum number of iterations
		:param tol: Tolerance for difference between two steps
		:return: Solution to the system of equations
		"""
		for k in range(k_max):
			x_old = np.copy(x)
			x -= npl.solve(jacobian(x), problem_function(x))
			if npl.norm(x-x_old) < tol:
				break
		return x


	def set_inner_bc(self, p, d1, d2):
		""" Compute the inter-arterial boundary conditions for one bifurcation.
		:param p: Parent artery
		:param d1: First daughter vessel
		:param d2: Second daughter vessel
		:return: Area and flow for the three vessels
		"""
		self.x = newton(p, d1, d2, self.x)
		p.U_out = Constant((self.x[9], self.x[0]))
		d1.U_in = Constant((self.x[12], self.x[3]))
		d2.U_in = Constant((self.x[15], self.x[6]))


	def define_geometry(self, Nx, Nt, T, N_cycles)
		"""Calls define_geometry on each artery.
		:param Nx: Number of spatial steps
		:param Nt: Number of temporal steps
		:param T: Period of one cardiac cycle
		:param N_cycles: Number of cardiac cycles
		"""
		for i in range(len(arteries)):
			arteries[i].define_geometry(Nx, Nt, T, N_cycles)


	def define_solution(self, q0):
		"""Computes q0 on each artery, befor calling define_solution.
		The daughter vessel gets a flow proportional to its share of the area.
		:param q0: Initial flow of the first vessel
		"""
		arteries[0].define_solution(q0)
		for i in range(1, len(arteries)):
			p = self.parent_vessel(i)
			s = self.sister_vessel(i)
			q0 = arteries[i].A0(0)/(arteries[i].A0(0)+arteries[s].A0(0))
			   * arteries[p].q0
			arteries[i].define_solution(q0)
		self.x = np.zeros([2**(self.order-1), 18]
		for i in range(2**(self.order-1)):
			#x[i] = <initial x>
			None


	def solve(self, q_ins):
		"""Solve the equation on the entire arterial network.
		:param q_ins: Vector containing inlet flow for the first artery.
		"""
		
		# Progress bar
		progress = Progress('Time-stepping')
		set_log_level(PROGRESS)

		# Initialise time
		t = 0

		# Cardiac cycle iteration
		for n_cycle in range(self.N_cycles):
			
			# Time-stepping for one period
			for n in range(self.Nt-1):

				print('Iteration '+str(n))


				# Update current solution
				for artery in arteries:
					artery.U_n.assign(artery.U)

				# Update inlet boundary conditions
				arteries[0].q_in.assign(Constant(q_ins[n+1]))
				
				# Update bifurcation boundary conditions
				for ip in range(2**(self.order-1):
					i1, i2 = daughter_vessels(ip)
					p = self.arteries[ip]
					d1, d2 = self.arteries[d1], self.arteries[d2]
					set_inner_bc(p, d1, d2)
				
				# Update outlet boundary condition
				A_out_value = self.compute_A_out(U_n)
				A_out.assign(Constant(A_out_value))
				
				# U_n+1 is solution of FF == 0
				solve(variational_form == 0, U, bcs)
				
				# Store solution at time t_(n+1)
				None
				
				t += self.dt
				
				# Update progress bar
				progress.update((t+self.dt)/self.N_cycles/self.T)
