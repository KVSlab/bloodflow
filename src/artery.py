__author__ = 'Syver Døving Agdestein'

import sys
import numpy as np
from fenics import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class Artery(object):
	""" Represent an artery.
	:param Ru: Upstream radius
	:param Rd: Downstream radius
	:param L: Vessel length
	:param k: Vector containing k1, k2 and k3 from the relation Eh/R0
	:param Re: Reynolds number
	:param p0: Diastolic pressure
	"""
	def __init__(self, root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0):
		""" Construct artery.
		Add its intrinsic characteristics, not its numerical solution.
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


	def CFL_term(self, x, A, q):
		"""Compute the term of the CFL condition.
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
		return self.dt/self.dex < self.CFL_term


	def adjust_dex(self, x, A, q, margin=1.01):
		"""Adjust boundary step-size so that the CFL condition is verified.
		:param x: Point at which the condition is to be checked
		:param A: Value of area at x
		:param q: Value of flow at x
		:param margin: A number greater than or equal to one
		"""
		M = self.CFL_term(x, A, q)
		if self.dt/self.dex > M:
			self.dex = margin*self.dt/M
		

	def define_geometry(self, Nx, Nt, T, N_cycles):
		"""Define FEniCS parameters.
		:param Nx: Number of spatial steps
		:param Nt: Number of temporal steps
		:param T: Duration of one cardiac cycle
		:param N_cycles: Number of cardiac cycles
		:param q0: Initial flow value.
		"""
		self.Nx = Nx
		self.Nt = Nt
		self.T = T
		self.N_cycles = N_cycles
		
		self.dt = self.T/self.Nt
		self.dx = self.L/self.Nx
		
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
		self.A0 = Expression('pi*pow(Ru, 2)*pow(Rd/Ru, 2*x[0]/L)',
							degree=2, Ru=self.Ru, Rd=self.Rd, L=self.L)
		self.f = Expression('4/3*(k1*exp(k2*Ru*pow(Ru/Rd, x[0]/L)) + k3)',
							degree=2, Ru=self.Ru, Rd=self.Rd, L=self.L,
							k1=self.k1, k2=self.k2, k3=self.k3)
		self.dfdr = Expression('4/3*k1*k2*exp(k2*Ru*pow(Rd/Ru, x[0]/L))',
						 	   degree=2, Ru=self.Ru, Rd=self.Rd, L=self.L,
							   k1=self.k1, k2=self.k2)
		self.drdx = Expression('log(Rd/Ru)/L*Ru*pow(Rd/Ru, x[0]/L)',
							   degree=2, Ru=self.Ru, Rd=self.Rd, L=self.L)


	def define_solution(self, q0, theta=0.5):
		
		# Define Crank-Nicolson parameter
		self.theta = theta
		
		# Define initial flow
		self.q0 = q0
		
		# Define trial function
		self.U = Function(self.V2)
		A, q = split(self.U)

		# Define test functions
		v1, v2 = TestFunctions(self.V2)
 
		# Inlet flow
		self.q_in = Function(self.V)
		self.q_in.assign(Constant(self.q0))
		
		# Inlet bc
		self.U_in = Function(self.V2)
		self.U_in.assign(Constant((self.Ru, self.q0)))

		# Outlet area
		self.A_out = Function(self.V)
		self.A_out.assign(Constant(self.A0(self.L)))

		# Outlet bc
		self.U_out = Function(self.V2)
		self.U_out.assign(Constant((self.Rd, self.q0)))

		# Initial value deduced from bottom boundary conditions
		self.Un = Function(self.V2)
		self.Un.assign(Expression(('pi*pow(Ru, 2)*pow(Rd/Ru, 2*x[0]/L)', 'q0'),
			degree=2, Ru=self.Ru, Rd=self.Rd, L=self.L, q0=self.q0))
		
		# Spatial boundary conditions
		tol = 1.e-14
		def inlet_bdry(x, on_boundary):
			return on_boundary and near(x[0], 0, tol)
		def outlet_bdry(x, on_boundary):
			return on_boundary and near(x[0], self.L, tol)

		if 1:#self.root_vessel:
			bc_inlet = DirichletBC(self.V2.sub(1), self.q_in, inlet_bdry)
		else:
			bc_inlet = DirichletBC(self.V2, self.U_in, inlet_bdry)
			
		if self.end_vessel:
			bc_outlet = DirichletBC(self.V2.sub(0), self.A_out, outlet_bdry)
		else:
			bc_outlet = DirichletBC(self.V2, self.U_out, outlet_bdry)

		self.bcs = [bc_inlet, bc_outlet]
		
		# Terms for variational form
		U_v_dx = A*v1*dx + q*v2*dx

		Un_v_dx = self.Un[0]*v1*dx + self.Un[1]*v2*dx


		F_v_ds = (pow(q, 2)/(A+DOLFIN_EPS)\
				 +self.f*sqrt(self.A0*(A+DOLFIN_EPS)))*v2*ds
		F_dv_dx = -grad(q)[0]*v1*dx\
				+ (pow(q, 2)/(A+DOLFIN_EPS)\
				  +self.f*sqrt(self.A0*(A+DOLFIN_EPS)))*grad(v2)[0]*dx
		dF_v_dx = F_v_ds - F_dv_dx

		Fn_v_ds = (pow(self.Un[1], 2)/(self.Un[0])\
				  +self.f*sqrt(self.A0*(self.Un[0])))*v2*ds
		Fn_dv_dx = -grad(self.Un[1])[0]*v1*dx\
				 + (pow(self.Un[1], 2)/(self.Un[0])\
				   +self.f*sqrt(self.A0*(self.Un[0])))*grad(v2)[0]*dx
		dFn_v_dx = Fn_v_ds - Fn_dv_dx


		S_v_dx = - 2*sqrt(pi)/self.db/self.Re*q/sqrt(A+DOLFIN_EPS)*v2*dx\
			   + (2*sqrt(A+DOLFIN_EPS)*(sqrt(pi)*self.f
									   +sqrt(self.A0)*self.dfdr)\
				 -(A+DOLFIN_EPS)*self.dfdr)*self.drdx*v2*dx

		Sn_v_dx = -2*sqrt(pi)/self.db/self.Re*self.Un[1]/sqrt(self.Un[0])*v2*dx\
				+ (2*sqrt(self.Un[0])*(sqrt(pi)*self.f+sqrt(self.A0)*self.dfdr)\
				  -(self.Un[0])*self.dfdr)*self.drdx*v2*dx

		# Variational form
		self.variational_form =\
			  U_v_dx\
		 	- Un_v_dx\
			+ self.dt*self.theta*dF_v_dx\
			+ self.dt*(1-self.theta)*dFn_v_dx\
			- self.dt*self.theta*S_v_dx\
			- self.dt*(1-self.theta)*Sn_v_dx


	def solve(self):
		"""Solve problem for one iteration.
		U_(n+1) is solution of variational_form == 0
		"""
		solve(self.variational_form == 0, self.U, self.bcs)


	def update_solution(self):
		"""Assign new values to U_n.
		"""
		self.Un.assign(self.U)


	def pressure(self, f, A0, A):
		""" Compute the pressure at a given point x and time t.
		:param f: Value of f(r0) in x
		:param A0: Value of A0 in x
		:param A: Area in x at a given time t
		:return: Pressure in x at time t
		"""
		return self.p0 + f*(1-np.sqrt(A0/A))


	def outlet_pressure(self, A):
		""" Compute the outlet pressure at a given time t.
		:param A: Area in L at a given time t
		:return: Pressure in L at time t
		"""
		return self.p0 + self.f(self.L)*(1-np.sqrt(self.A0(self.L)/A))
