__author__ = 'Syver Døving Agdestein'

import sys
import numpy as np
import numpy.linalg as npl

import configparser
from fenics import *
from mshr import *

from artery import Artery


set_log_level(30)


class Artery_Network(object):
	"""
	:param order: Number of arterial levels
	:param rc: Characteristic radius (length)
	:param qc: Characteristic flow
	:param Ru: Upstream radii
	:param Rd: Downstream radii
	:param L: Vessel lengths
	:param k1: First constant from the relation Eh/r0
	:param k2: Second constant from the relation Eh/r0
	:param k3: Third constant from the relation Eh/R0
	:param rho: Density of blood
	:param Re: Reynolds' number
	:param nu: Viscosity of blood
	:param p0: Diastolic pressure
	:param R1: First resistance from Windkessel model
	:param R2: Second resistance from Windkessel model
	:param CT: Compliance from Windkessel model
	"""
	def __init__(self, order, rc, qc, Ru, Rd, L, k1, k2, k3,
				 rho, Re, nu, p0, R1, R2, CT):
		"""Artery network constructor.
		"""
		self.order = order
		self.arteries = [0] * (2**self.order-1)
		self.range_arteries = range(2**self.order-1)
		self.range_parent_arteries = range(2**(self.order-1)-1)
		self.range_daughter_arteries = range(1, 2**self.order-1)
		self.range_end_arteries = range(2**(self.order-1)-1, 2**self.order-1)
		self.rc, self.qc, self.rho = rc, qc, rho
		self.R1, self.R2, self.CT = R1, R2, CT
		for i in self.range_arteries:
			root_vessel = (i==0)
			end_vessel = (i in self.range_end_arteries)
			self.arteries[i] = Artery(root_vessel, end_vessel, rc, qc, Ru[i],
									  Rd[i], L[i], k1, k2, k3, rho, Re, nu, p0)


	def daughter_arteries(self, i):
		"""Find daughter vessels.
		:param i: Parent vessel index
		:return: Daughter vessel indices
		"""
		#if i < 0 or i >= 2**(self.order-1):
		#	raise Exception('Vessel index out of range')

		return 2*i+1, 2*i+2


	def parent_artery(self, i):
		"""Find parent vessel.
		:param i: Daughter vessel index
		:return: Parent vessel index
		"""
		#if i <= 0 or i >= 2**self.order:
		#	raise Exception('Vessel index out of range')

		return (i-1)//2  # d1 is odd, d2=d1+1 is pair

	def sister_artery(self, i):
		"""Find sister vessel.
		:param i: Vessel index
		:return: Sister vessel index
		"""
		if i%2 == 0:
			return i-1
		else:
			return i+1


	def define_geometry(self, Nx, Nt, T, N_cycles):
		"""Define FEniCS geometry on the entire arterial network.
		:param Nx: Spatial discretisation number
		:param Nt: Number of time-steps
		:param T: Duration of one cardiac cycle
		:param N_cycles: Number of cardiac cycles
		"""
		self.Nx = Nx
		self.Nt = Nt
		self.T = T
		self.N_cycles = N_cycles
		self.dt = T/Nt
		for i in self.range_arteries:
			self.arteries[i].define_geometry(Nx, Nt, T, N_cycles)


	def define_solution(self, output_location, q0, theta=0.5):
		"""Computes q0 on each artery, before calling define_solution.
		The daughter vessel gets a flow proportional to its share of the area.
		:param string output_location: Output location
		:param q0: Initial flow of the root vessel
		:param theta: Crank-Nicolson weight parameter, in the interval [0, 1]
		"""
		self.output_location = output_location
		self.theta = theta
		self.arteries[0].define_solution(q0, theta)
		for i in self.range_daughter_arteries:
			p = self.parent_artery(i)
			s = self.sister_artery(i)
			q0 = self.arteries[i].A0(0)/(self.arteries[i].A0(0)\
										+self.arteries[s].A0(0))\
			   * self.arteries[p].q0
			self.arteries[i].define_solution(q0, theta)

		self.x = np.ones([len(self.range_parent_arteries), 18])


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
		   + (2*np.sqrt(U[0])*(np.sqrt(np.pi)*a.f(x)\
							  +np.sqrt(a.A0(x))*a.dfdr(x))\
			 -U[0]*a.dfdr(x))*a.drdx(x)
		return np.array([S1, S2])


	def compute_U_half(self, a, x0, x1, U0, U1):
		"""Compute solution at half step.
		:param a: Current artery
		:param x0: Left point
		:param x1: Right point
		:param U0: Solution at right point
		:param U1: Solution at left point
		:return: Middle point, half a time step later
		"""
		# Value of terms at time t_n
		F0, S0 = self.flux(a, U0, x0), self.source(a, U0, x0)
		F1, S1 = self.flux(a, U1, x1), self.source(a, U1, x1)

		return (U0+U1)/2 - a.dt/(x1-x0)*(F1-F0) + a.dt/4*(S0+S1)


	def compute_A_out(self, a, k_max=100, tol=1.0e-12):
		"""Compute the outlet boundary condition.
		:param a: Artery on which the outlet area is to be computed
		:param k_max: Maximum number of iterations in Piccards scheme
		:param tol: Tolerance for Piccards fixed point iteration scheme
		:return: Outlet boundary value of A at time step t_(n+1)
		"""
		a.adjust_dex(a.L, a.Un(a.L)[0], a.Un(a.L)[1])

		# Spatial step, scaled to satisfy the CFL condition
		x2, x1, x0 = a.L-2*a.dex, a.L-a.dex, a.L
		x21, x10 = a.L-1.5*a.dex, a.L-0.5*a.dex
		Um2, Um1, Um0 = a.Un(x2), a.Un(x1), a.Un(x0)

		# Values at time step n
		Fm2, Sm2 = self.flux(a, Um2, x2), self.source(a, Um2, x2)
		Fm1, Sm1 = self.flux(a, Um1, x1), self.source(a, Um1, x1)
		Fm0, Sm0 = self.flux(a, Um0, x0), self.source(a, Um0, x0)

		# Values at time step n+1/2
		U_half_21 = self.compute_U_half(a, x2, x1, Um2, Um1)
		U_half_10 = self.compute_U_half(a, x1, x0, Um1, Um0)
		F_half_21 = self.flux(a, U_half_21, x21)
		S_half_21 = self.source(a, U_half_21, x21)
		F_half_10 = self.flux(a, U_half_10, x10)
		S_half_10 = self.source(a, U_half_10, x10)

		# Value at time step n+1
		qm1 = Um1[1]\
			- a.dt/a.dex*(F_half_10[1]-F_half_21[1])\
			+ a.dt/2*(S_half_10[1]+S_half_21[1])

		# Fixed point iteration
		pn = a.compute_outlet_pressure(Um0[0])
		p = pn
		for k in range(k_max):
			p_old = p
			qm0 = Um0[1]\
				+ (p-pn)/self.R1\
				+ self.dt/self.R1/self.R2/self.CT*pn\
				- self.dt*(self.R1+self.R2)/self.R1/self.R2/self.CT*Um0[1]
			Am0 = Um0[0] - self.dt/a.dex*(qm0-qm1)
			p = a.compute_outlet_pressure(Am0)
			if abs(p-p_old) < tol:
				break

		return Am0


	def initial_x(self, p, d1, d2):
		""" Make an initial guess for x at a bifurcation point.
		Set same value at time t_(n+1) and t_(n+1/2) as time t_n.
		At point M+-1/2, set same value as in point M.
		:param p: Parent artery
		:param d1: First daughter artery
		:param d2: Second daughter artery
		:return: x, an 18-dimensional vector containing guess values
		"""
		x = np.zeros(18)
		x[:3] = p.q0*np.ones(3)
		x[3:6] = d1.q0*np.ones(3)
		x[6:9] = d2.q0*np.ones(3)
		x[9:12] = p.A0(p.L)*np.ones(3)
		x[12:15] = d1.A0(0)*np.ones(3)
		x[15:] = d2.A0(0)*np.ones(3)
		return x


	def define_x(self):
		"""Make an initial guess for the solutions to the systems of equations.
		"""
		for ip in self.range_parent_arteries:
			i1, i2 = self.daughter_arteries(ip)
			p, d1, d2 = self.arteries[ip], self.arteries[i1], self.arteries[i2]
			self.x[ip] = self.initial_x(p, d1, d2)


	def problem_function(self, p, d1, d2, x):
		"""Function representing the system of equations.
		If x is the solution to the problem, then function(x) = 0.
		:param p: Parent artery
		:param d1: First daughter vessel
		:param d2: Second daughter vessel
		:param x: Current solution, an 18-dimensional vector
		:return: Value of function in x
		"""
		# Abbreviations
		A0p, A01, A02 = p.A0(p.L), d1.A0(0), d2.A0(0)
		fp, f1, f2 = p.f(p.L),  d1.f(0), d2.f(0)
		dbp, db1, db2 = p.db, d1.db, d2.db
		Rep, Re1, Re2 = p.Re, d1.Re, d2.Re
		dfdrp, dfdr1, dfdr2 = p.dfdr(p.L), d1.dfdr(0), d2.dfdr(0)
		drdxp, drdx1, drdx2 = p.drdx(p.L), d1.drdx(0), d2.drdx(0)
		rpi = sqrt(np.pi)

		# Ghost half terms
		Fp = self.flux(p, np.array([x[11], x[2]]), p.L+p.dex/2)
		F1 = self.flux(d1, np.array([x[14], x[5]]), -d1.dex/2)
		F2 = self.flux(d2, np.array([x[17], x[8]]), -d2.dex/2)
		Sp = self.source(p, np.array([x[11], x[2]]), p.L+p.dex/2)
		S1 = self.source(d1, np.array([x[14], x[5]]), -d1.dex/2)
		S2 = self.source(d2, np.array([x[17], x[8]]), -d2.dex/2)

		# Compute half-time-step-values in M-1/2 for p and 1/2 for d1 and d2
		Um1p, Um0p = p.Un(p.L-p.dex), p.Un(p.L)
		U0d1, U1d1 = d1.Un(0), d1.Un(d1.dex)
		U0d2, U1d2 = d2.Un(0), d2.Un(d2.dex)

		U_half_p = self.compute_U_half(p, p.L-p.dex, p.L, Um1p, Um0p)
		U_half_1 = self.compute_U_half(d1, 0, d1.dex, U0d1, U1d1)
		U_half_2 = self.compute_U_half(d2, 0, d2.dex, U0d2, U1d2)

		F_half_p = self.flux(p, U_half_p, p.L-p.dex/2)
		S_half_p = self.source(p, U_half_p, p.L-p.dex/2)
		F_half_1 = self.flux(d1, U_half_1, d1.dex/2)
		S_half_1 = self.source(d1, U_half_1, d1.dex/2)
		F_half_2 = self.flux(d2, U_half_2, d2.dex/2)
		S_half_2 = self.source(d2, U_half_2, d2.dex/2)

		# Function return array
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
		y[6] = x[0] - x[3] - x[6]
		y[7] = x[1] - x[4] - x[7]

		# Entries from equation (23)
		y[8] = fp*(1-np.sqrt(A0p/x[10])) - f1*(1-np.sqrt(A01/x[13]))
		y[9] = fp*(1-np.sqrt(A0p/x[10])) - f2*(1-np.sqrt(A02/x[16]))
		y[10] = fp*(1-np.sqrt(A0p/x[9])) - f1*(1-np.sqrt(A01/x[12]))
		y[11] = fp*(1-np.sqrt(A0p/x[9])) - f2*(1-np.sqrt(A02/x[15]))

		# Entries from equation (26)
		y[12] = x[0] - Um0p[1] + p.dt/p.dex*(Fp[1] - F_half_p[1])\
			  - p.dt/2*(Sp[1] + S_half_p[1])
		y[13] = x[3] - U0d1[1] + d1.dt/d1.dex*(F_half_1[1] - F1[1])\
			  - d1.dt/2*(S_half_1[1] + S1[1])
		y[14] = x[6] - U0d2[1] + d2.dt/d2.dex*(F_half_2[1] - F2[1])\
			  - d2.dt/2*(S_half_2[1] + S2[1])

		# Entries from equation (27)
		y[15] = x[9] - Um0p[0] + p.dt/p.dex*(Fp[0] - F_half_p[0])
		y[16] = x[12] - U0d1[0] + d1.dt/d1.dex*(F_half_1[0] - F1[0])
		y[17] = x[15] - U0d2[0] + d2.dt/d2.dex*(F_half_2[0] - F2[0])

		return y


	def jacobian(self, p, d1, d2, x):
		"""Compute the analytical jacobian matrix for the system of equations.
		:param p: Parent artery
		:param d1: First daughter vessel
		:param d2: Second daughter vessel
		:param x: Current solution, an 18-dimensional vector
		:return: Jacobian matrix
		"""
		# Abbreviations
		A0p, A01, A02 = p.A0(p.L), d1.A0(0), d2.A0(0)
		A0hp, A0h1, A0h2 = p.A0(p.L+p.dex), d1.A0(-d1.dex), d2.A0(-d2.dex)
		fp, f1, f2 = p.f(p.L),  d1.f(0), d2.f(0)
		fhp, fh1, fh2 = p.f(p.L+p.dex), d1.f(-d1.dex), d2.f(-d2.dex)
		dbp, db1, db2 = p.db, d1.db, d2.db
		Rep, Re1, Re2 = p.Re, d1.Re, d2.Re
		dfdrp, dfdr1, dfdr2 = p.dfdr(p.L+p.dex),\
							  d1.dfdr(-d1.dex), d2.dfdr(-d2.dex)
		drdxp, drdx1, drdx2 = p.drdx(p.L+p.dex),\
							  d1.drdx(-d1.dex), d2.drdx(-d2.dex)
		rpi = np.sqrt(np.pi)

		# Jacobian matrix
		J = np.zeros([18, 18])

		# Entries from equation (20)
		J[0, 1] = 2
		J[0, 2] = -1
		J[1, 4] = 2
		J[1, 5] = -1
		J[2, 7] = 2
		J[2, 8] = -1

		# Entries from equation (21)
		J[3, 10] = 2
		J[3, 11] = -1
		J[4, 13] = 2
		J[4, 14] = -1
		J[5, 16] = 2
		J[5, 17] = -1

		# Entries from equation (22)
		J[6, 0] = 1
		J[6, 3] = -1
		J[6, 6] = -1
		J[7, 1] = 1
		J[7, 4] = -1
		J[7, 7] = -1

		# Entries from equation (23)
		J[8, 10] = fp*np.sqrt(A0p)/2/x[10]**(3/2)
		J[8, 13] = -f1*np.sqrt(A01)/2/x[13]**(3/2)
		J[9, 10] = fp*np.sqrt(A0p)/2/x[10]**(3/2)
		J[9, 16] = -f2*np.sqrt(A02)/2/x[16]**(3/2)
		J[10, 9] = fp*np.sqrt(A0p)/2/x[9]**(3/2)
		J[10, 12] = -f1*np.sqrt(A01)/2/x[12]**(3/2)
		J[11, 9] = fp*np.sqrt(A0p)/2/x[9]**(3/2)
		J[11, 15] = -f2*np.sqrt(A02)/2/x[15]**(3/2)

		# Entries from equation (26)
		J[12, 0] = 1
		J[12, 2] = p.dt/p.dex*2*x[2]/x[11] + p.dt*rpi/dbp/Rep/np.sqrt(x[11])
		J[12, 11] = p.dt/p.dex*(-(x[2]/x[11])**2 + fhp/2*np.sqrt(A0hp/x[11]))\
				  - p.dt/2*(rpi/dbp/Rep*x[2]/x[11]**(3/2)\
						   + (1/np.sqrt(x[11])*(rpi*fhp+np.sqrt(A0hp)*dfdrp)\
											   - dfdrp)*drdxp)
		J[13, 3] = 1
		J[13, 5] = -d1.dt/d1.dex*2*x[5]/x[14] + d1.dt*rpi/db1/Re1/np.sqrt(x[14])
		J[13, 14] = d1.dt/d1.dex*((x[5]/x[14])**2 - fh1/2*np.sqrt(A0h1/x[14]))\
				  - d1.dt/2*(rpi/db1/Re1*x[5]/x[14]**(3/2)\
							+ (1/np.sqrt(x[14])*(rpi*fh1+np.sqrt(A0h1)*dfdr1)\
												- dfdr1)*drdx1)
		J[14, 6] = 1
		J[14, 8] = -d2.dt/d2.dex*2*x[8]/x[17] + d2.dt*rpi/db2/Re2/np.sqrt(x[17])
		J[14, 17] = d2.dt/d2.dex*((x[8]/x[17])**2 - fh2/2*np.sqrt(A0h2/x[17]))\
				  - d2.dt/2*(rpi/db2/Re2*x[8]/x[17]**(3/2)\
							+ (1/np.sqrt(x[17])*(rpi*fh2+np.sqrt(A0h2)*dfdr2)\
												- dfdr2)*drdx2)

		# Entries from equation (27)
		J[15, 2] = p.dt/p.dex
		J[15, 9] = 1
		J[16, 5] = - d1.dt/d1.dex
		J[16, 12] = 1
		J[17, 8] = - d1.dt/d1.dex
		J[17, 15] = 1

		return J


	def newton(self, p, d1, d2, x=np.ones(18), k_max=30, tol=1.e-10):
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

			J = self.jacobian(p, d1, d2, x)
			func = self.problem_function(p, d1, d2, x)

			if npl.norm(func) < tol:
				break

			try:
				x -= npl.solve(J, func)
			except npl.LinAlgError:
				print('Singular')
				eps = 1.e-6  # Perturbation value
				J += eps*np.eye(18)
				func[0] += eps
				x -= npl.solve(J, func)

		return x


	def adjust_bifurcation_step(self, p, d1, d2, margin=0.05):
		"""Set dex to respect CFL-condition for all arteries in a bifurcation.
		:param p: Parent artery
		:param d1: First daughter artery
		:param d2: Second daughter artery
		:param margin: Number greater than or equal to one
		"""
		Mp = p.CFL_term(p.L, p.Un(p.L)[0], p.Un(p.L)[1])
		M1 = d1.CFL_term(0, d1.Un(0)[0], d1.Un(0)[1])
		M2 = d2.CFL_term(0, d2.Un(0)[0], d2.Un(0)[1])

		# dex is chosen to respect all three CFL-conditions
		p.dex = d1.dex = d2.dex = (1+margin)*self.dt/min([Mp, M1, M2])


	def set_inner_bc(self, ip, i1, i2):
		"""Compute the inter-arterial boundary conditions for one bifurcation.
		:param ip: Parent artery index
		:param i1: First daughter vessel index
		:param i2: Second daughter vessel index
		"""
		p, d1, d2 = self.arteries[ip], self.arteries[i1], self.arteries[i2]

		self.adjust_bifurcation_step(p, d1, d2)

		self.x[ip] = self.newton(p, d1, d2, self.x[ip])

		p.U_out = [self.x[ip, 9], self.x[ip, 0]]
		d1.U_in = [self.x[ip, 12], self.x[ip, 3]]
		d2.U_in = [self.x[ip, 15], self.x[ip, 6]]


	def set_bcs(self, q_in):
		""" Update boundary conditions for time t_(n+1) at all boundaries.
		:param q_in: Value of inflow for root artery at time t_(n+1)
		"""
		# Update inlet boundary conditions
		self.arteries[0].q_in = q_in

		# Update bifurcation boundary conditions
		for ip in self.range_parent_arteries:
			i1, i2 = self.daughter_arteries(ip)
			self.set_inner_bc(ip, i1, i2)

		# Update outlet boundary conditions
		for i in self.range_end_arteries:
			self.arteries[i].A_out = self.compute_A_out(self.arteries[i])


	def dump_metadata(self, Nt_store, N_cycles_store, store_area,
					  store_pressure):
		"""Save mesh.
		Dump metadata necessary for interpretation of xdmf-files.
		:param boolean store_area: Store area if True
		:param boolean store_pressure: Store pressure if True
		"""
		# Assemble strings
		mesh_locations = ''
		for i in self.range_arteries:
			mesh_location = self.output_location + '/mesh_%i.xml.gz' % (i)
			File(mesh_location) << self.arteries[i].mesh  # Save mesh
			if i > 0:
				mesh_locations += ','
			mesh_locations += mesh_location
		L = ''
		for artery in self.arteries[:-1]:
			L += str(artery.L)+','
		L += str(self.arteries[-1].L)
		names = ''
		locations = ''
		names += 'flow'
		locations += self.output_location + '/flow'
		if store_area:
			names += ',area'
			locations += ',' + self.output_location + '/area'
		if store_pressure:
			names += ',pressure'
			locations += ',' + self.output_location + '/pressure'

		# Save metadata
		config = configparser.RawConfigParser()
		config.add_section('data')
		config.set('data', 'order', str(self.order))
		config.set('data', 'Nx', str(self.Nx))
		config.set('data', 'Nt', str(Nt_store*N_cycles_store))
		config.set('data', 'T0', str(self.T*(self.N_cycles-N_cycles_store)))
		config.set('data', 'T', str(self.T*self.N_cycles))
		config.set('data', 'L', L)
		config.set('data', 'rc', str(self.rc))
		config.set('data', 'qc', str(self.qc))
		config.set('data', 'rho', str(self.rho))
		config.set('data', 'mesh_locations', mesh_locations)
		config.set('data', 'names', names)
		config.set('data', 'locations', locations)
		with open(self.output_location+'/data.cfg', 'w') as configfile:
			config.write(configfile)


	def solve(self, q_ins, Nt_store, N_cycles_store=1, store_area=False,
			  store_pressure=True):
		"""Solve the equation on the entire arterial network.
		:param q_ins: Vector containing inlet flow for the first artery
		:param Nt_store: Number of time-steps to be stored per cycle
		:param N_cycles_store: Number of cycles to be stored (starting at last)
		:param boolean store_area: Store area if true
		:param boolean store_pressure: Store pressure if true
		"""
		self.define_x()

		# Store parameters necessary for postprocessing
		self.dump_metadata(Nt_store, N_cycles_store, store_area, store_pressure)

		# Setup storage files
		xdmffile_flow = [0] * len(self.range_arteries)
		if store_area:
			xdmffile_area = [0] * len(self.range_arteries)
		if store_pressure:
			xdmffile_pressure = [0] * len(self.range_arteries)
		for i in self.range_arteries:
			xdmffile_flow[i] = XDMFFile('%s/flow/flow_%i.xdmf'\
										% (self.output_location, i))
			if store_area:
				xdmffile_area[i] = XDMFFile('%s/area/area_%i.xdmf'\
											% (self.output_location, i))
			if store_pressure:
				xdmffile_pressure[i] = XDMFFile('%s/pressure/pressure_%i.xdmf'\
												% (self.output_location, i))

		# Initialise time
		t = 0

		# Cardiac cycle iteration
		for n_cycle in range(self.N_cycles):

			# Time-stepping for one period
			for n in range(self.Nt):

				print('Current cycle: %i\nCycle iteration: %i\nTime-step t_%i'\
					  % (n_cycle, n, n_cycle*self.Nt+n))

				# Apply boundary conditions for time t_(n+1)
				self.set_bcs(q_ins[(n+1) % (self.Nt)])

				# Solve equation on each artery
				for i, artery in enumerate(self.arteries):

					# Store solution at time t_n
					cycle_store = (n_cycle >= self.N_cycles-N_cycles_store)

					if cycle_store and n % (self.Nt/Nt_store) == 0:

						# Split solution for storing, with deepcopy
						area, flow = artery.Un.split(True)

						xdmffile_flow[i].write_checkpoint(flow, 'flow', t)

						if store_area:
							xdmffile_area[i].write_checkpoint(area, 'area', t)

						if store_pressure:
							artery.update_pressure()
							xdmffile_pressure[i].write_checkpoint(artery.pn,
																  'pressure', t)

					# Solve problem on artery for time t_(n+1)
					artery.solve()

					# Update current solution on artery
					artery.update_solution()

				t += self.dt
