__author__ = 'Syver Døving Agdestein'

import sys
import numpy as np

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from fenics import *


class Artery(object):
	"""Represent an artery.
	:param boolean root_vessel: True if the artery is root-vessel (no parent)
	:param boolean end_vessel: True if the artery is end-vessel (no daughter)
	:param rc: Characteristic radius (length)
	:param qc: Characteristic flow
	:param Ru: Upstream radius
	:param Rd: Downstream radius
	:param L: Vessel length
	:param k1: First constant from the relation Eh/r0
	:param k2: Second constant from the relation Eh/r0
	:param k3: Third constant from the relation Eh/R0
	:param rho: Density of blood
	:param Re: Reynolds' number
	:param nu: Blood viscosity
	:param p0: Diastolic pressure
	"""
	def __init__(self, root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3,
				 rho, Re, nu, p0):
		""" Construct artery.
		Add its intrinsic characteristics.
		"""
		self.root_vessel = root_vessel
		self.end_vessel = end_vessel
		self.rc = rc
		self.qc = qc
		self.Ru = Ru
		self.Rd = Rd
		self.L = L
		self.k1 = k1
		self.k2 = k2
		self.k3 = k3
		self.rho = rho
		self.Re = Re
		self.nu = nu
		self.p0 = p0


	def define_geometry(self, Nx, Nt, T, N_cycles):
		"""Initialise geometry.
		Define FEniCS objects.
		:param Nx: Number of spatial steps
		:param Nt: Number of temporal steps
		:param T: Duration of one cardiac cycle
		:param N_cycles: Number of cardiac cycles
		"""
		self.Nx = Nx
		self.Nt = Nt
		self.T = T
		self.N_cycles = N_cycles
		
		self.dx = self.L/self.Nx
		self.dt = self.T/self.Nt
		
		# Step for boundary condition computations, starting at normal size
		self.dex = self.dx
		
		self.db = np.sqrt(self.nu*self.T/2/np.pi)
		
		self.mesh = IntervalMesh(self.Nx, 0, self.L)
		self.elV = FiniteElement('CG', self.mesh.ufl_cell(), 1)
		self.V = FunctionSpace(self.mesh, self.elV)
		self.V2 = FunctionSpace(self.mesh, self.elV*self.elV)

		# Initial vessel-radius and deduced quantities
		self.r0 = Expression('Ru*pow(Rd/Ru, x[0]/L)',
							degree=2, Ru=self.Ru, Rd=self.Rd, L=self.L)
		self.A0 = Expression('pi*pow(r0, 2)', degree=2, r0=self.r0)
		self.f = Expression('4.0/3.0*(k1*exp(k2*r0) + k3)', degree=2,
							k1=self.k1, k2=self.k2, k3=self.k3, r0=self.r0)
		self.dfdr = Expression('4.0/3.0*k1*k2*exp(k2*r0)', degree=2,
						 	   k1=self.k1, k2=self.k2, r0=self.r0)
		self.drdx = Expression('log(Rd/Ru)/L*Ru*pow(Rd/Ru, x[0]/L)', degree=2,
							   Ru=self.Ru, Rd=self.Rd, L=self.L)


	def define_solution(self, q0, theta=0.5, bc_tol=1.e-14):
		"""Set up FEniCS solution objects.
		Define boundary conditions.
		Define variational form.
		:param q0: Initial flow
		:param theta: Crank-Nicolson parameter
		:param bc_tol: Inlet and outlet boundary thickness (tolerance)
		"""
		# Initial flow value
		self.q0 = q0
		
		# Crank-Nicolson parameter
		self.theta = theta
		
		# Trial function
		self.U = Function(self.V2)
		A, q = split(self.U)

		# Test functions
		v1, v2 = TestFunctions(self.V2)
 
		# Initial solution
		self.Un = Function(self.V2)
		self.Un.assign(Expression(('A0', 'q0'), degree=2,
								  A0=self.A0, q0=self.q0))
		
		# Initial pressure
		self.pn = Function(self.V)
		self.pn.assign(Expression('f*(1-sqrt(A0/A))', degree=2,
								  f=self.f, A0=self.A0, A=self.Un.split()[0]))
		
		# Boundary conditions (spatial)
		def inlet_bdry(x, on_boundary):
			return on_boundary and near(x[0], 0, bc_tol)

		def outlet_bdry(x, on_boundary):
			return on_boundary and near(x[0], self.L, bc_tol)

		# Inlet boundary conditions
		if self.root_vessel:
			self._q_in = Expression('value', degree=0, value=self.q0)
			bc_inlet = DirichletBC(self.V2.sub(1), self._q_in, inlet_bdry)
		else:
			self._U_in = Expression(('A', 'q'), degree = 0,
									A=self.A0(0), q=self.q0)
			bc_inlet = DirichletBC(self.V2, self._U_in, inlet_bdry)
			
		# Outlet boundary conditions
		if self.end_vessel:
			self._A_out = Expression('value', degree=0, value=self.A0(self.L))
			bc_outlet = DirichletBC(self.V2.sub(0), self._A_out, outlet_bdry)
		else:
			self._U_out = Expression(('A', 'q'), degree=0,
									 A=self.A0(self.L), q=self.q0)
			bc_outlet = DirichletBC(self.V2, self._U_out, outlet_bdry)

		self.bcs = [bc_inlet, bc_outlet]
		
		# Terms for variational form
		U_v_dx = A*v1*dx + q*v2*dx

		Un_v_dx = self.Un[0]*v1*dx + self.Un[1]*v2*dx


		F2_v2_ds = (pow(q, 2)/(A+DOLFIN_EPS)\
				   +self.f*sqrt(self.A0*(A+DOLFIN_EPS)))*v2*ds
		
		F2_dv2_dx = (pow(q, 2)/(A+DOLFIN_EPS)\
					+self.f*sqrt(self.A0*(A+DOLFIN_EPS)))*grad(v2)[0]*dx

		dF_v_dx = grad(q)[0]*v1*dx + F2_v2_ds - F2_dv2_dx


		Fn_v_ds = (pow(self.Un[1], 2)/(self.Un[0])\
				  +self.f*sqrt(self.A0*(self.Un[0])))*v2*ds

		Fn_dv_dx = (pow(self.Un[1], 2)/(self.Un[0])\
				   +self.f*sqrt(self.A0*(self.Un[0])))*grad(v2)[0]*dx

		dFn_v_dx = grad(self.Un[1])[0]*v1*dx + Fn_v_ds - Fn_dv_dx


		S_v_dx = - 2*sqrt(pi)/self.db/self.Re*q/sqrt(A+DOLFIN_EPS)*v2*dx\
			   + (2*sqrt(A+DOLFIN_EPS)*(sqrt(pi)*self.f
									   +sqrt(self.A0)*self.dfdr)\
				 -(A+DOLFIN_EPS)*self.dfdr)*self.drdx*v2*dx

		Sn_v_dx = -2*sqrt(pi)/self.db/self.Re*self.Un[1]/sqrt(self.Un[0])*v2*dx\
				+ (2*sqrt(self.Un[0])*(sqrt(pi)*self.f+sqrt(self.A0)*self.dfdr)\
				  -(self.Un[0])*self.dfdr)*self.drdx*v2*dx

		# Variational form
		self.variational_form = U_v_dx\
		 	- Un_v_dx\
			+ self.dt*self.theta*dF_v_dx\
			+ self.dt*(1-self.theta)*dFn_v_dx\
			- self.dt*self.theta*S_v_dx\
			- self.dt*(1-self.theta)*Sn_v_dx


	def solve(self):
		"""Solve problem for one iteration.
		U_(n+1) is solution of variational_form(U_n) == 0.
		"""
		solve(self.variational_form == 0, self.U, self.bcs)


	def update_solution(self):
		"""Assign new values to Un.
		"""
		self.Un.assign(self.U)


	def update_pressure(self):
		"""Assign new values to pn.
		"""
		self.pn.assign(Expression('p0 + f*(1-sqrt(A0/A))', degree=2, p0=self.p0,
								  f=self.f, A0=self.A0, A=self.Un.split(True)[0]))


	def compute_pressure(self, f, A0, A):
		""" Compute the pressure at a given point x and time t.
		:param f: Value of f(r0) in x
		:param A0: Value of A0 in x
		:param A: Area in x at a given time t
		:return: Pressure in x at time t
		"""
		return self.p0 + f*(1-np.sqrt(A0/A))


	def compute_outlet_pressure(self, A):
		""" Compute the outlet pressure at a given time t.
		:param A: Area in L at a given time t
		:return: Pressure in L at time t
		"""
		return self.p0 + self.f(self.L)*(1-np.sqrt(self.A0(self.L)/A))


	def CFL_term(self, x, A, q):
		"""Compute the term for the CFL condition.
		:param x: Point at which the condition is to be checked
		:param A: Value of area at x
		:param q: Value of flow at x
		:return: CFL term
		"""
		return 1/np.abs(q/A+np.sqrt(self.f(x)/2/self.rho\
								   *np.sqrt(self.A0(x)/A)))


	def check_CFL(self, x, A, q):
		"""Check the CFL condition.
		:param x: Point at which the condition is to be checked
		:param A: Value of area at x
		:param q: Value of flow at x
		:return: True if condition is verified
		"""
		return self.dt/self.dex < self.CFL_term(x, A, q)


	def adjust_dex(self, x, A, q, margin=0.05):
		"""Adjust boundary step-size so that the CFL condition is verified.
		:param x: Point at which the condition is to be checked
		:param A: Value of area at x
		:param q: Value of flow at x
		:param margin: A number greater than or equal to one
		"""
		M = self.CFL_term(x, A, q)
		#if self.dt/self.dex > M:
		self.dex = (1+margin)*self.dt/M

		
	@property
	def q_in(self):
		"""Inlet flow (only for root-artery)
		:return: Inlet flow Expression-object
		"""
		return self._q_in.value
	
	@q_in.setter
	def q_in(self, value):
		self._q_in.value = value


	@property
	def U_in(self):
		"""Inlet boundary conditions (only for non-root arteries (daughters))
		:return: Inlet solution value
		"""
		return np.array([self._U_in.A, self._U_in.q])
	
	@U_in.setter
	def U_in(self, U):
		self._U_in.A = U[0]
		self._U_in.q = U[1]

	
	@property
	def A_out(self):
		"""Outlet area (only in use for end arteries)
		:return: Outlet flow value
		"""
		return self._A_out.value
	
	@A_out.setter
	def A_out(self, value):
		self._A_out.value = value


	@property
	def U_out(self):
		"""Outlet boundary conditions (only for parent arteries)
		:return: Outlet solution Expression-object
		"""
		return np.array([self._U_out.A, self._U_out.q])
	
	@U_out.setter
	def U_out(self, U):
		self._U_out.A = U[0]
		self._U_out.q = U[1]
