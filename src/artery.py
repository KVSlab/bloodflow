__author__ = 'Syver Døving Agdestein'


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
	def __init__(self, Ru, Rd, L, k1, k2, k3, nu, p0, R1, R2, CT):
		""" Construct artery.
		Add its intrinsic characteristics, not its numerical solution.
		"""
		self.Ru = Ru
		self.Rd = Rd
		self.L = L
		self.k1 = k1
		self.k2 = k2
		self.k3 = k3
		self.nu = nu
		self.Re = 10.0/nu/1.0
		self.p0 = p0
		self.R1 = R1
		self.R2 = R2
		self.CT = CT


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
		self.q0 = q0
		
		self.dt = self.T/self.Nt
		self.dex = self.L/self.Nx
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


	def define_solution(self, q0):
	
		# Define trial function
		self.U = Function(self.V2)
		self.A, self.q = split(self.U)

		# Define test functions
		self.v1, self.v2 = TestFunctions(self.V2)
 
		# Inlet flow
		self.q_in = Function(self.V)
		self.q_in.assign(Constant(self.q0))

		# Outlet area
		self.A_out = Function(self.V)
		self.A_out.assign(Constant(self.A0(self.L)))

		# Initial value deduced from bottom boundary conditions
		self.U_n = Function(self.V2)
		self.U_n.assign(Expression(('pi*pow(Ru, 2)*pow(Rd/Ru, 2*x[0]/L)', 'q0'),
			degree=2, Ru=self.Ru, Rd=self.Rd, L=self.L, q0=self.q0))

		# Spatial boundary conditions
		tol = 1.e-14
		def inlet_bdry(x, on_boundary):
			return on_boundary and near(x[0], 0, tol)
		def outlet_bdry(x, on_boundary):
			return on_boundary and near(x[0], self.L, tol)
		bc_outlet = DirichletBC(self.V2.sub(0), self.A_out, outlet_bdry)
		bc_inlet = DirichletBC(self.V2.sub(1), self.q_in, inlet_bdry)
		self.bcs = [bc_inlet, bc_outlet]

		# Variational form
		self.variatonal_form = A*self.v1*dx\
			+ q*self.v2*dx\
			+ self.dt*grad(q)[0]*self.v1*dx\
			+ self.dt*(pow(q, 2)/(A+DOLFIN_EPS)\
				  +self.f*sqrt(self.A0*(A+DOLFIN_EPS)))*self.v2*ds\
			- self.dt*(pow(q, 2)/(A+DOLFIN_EPS)\
				  +self.f*sqrt(self.A0*(A+DOLFIN_EPS)))*grad(self.v2)[0]*dx\
			+ self.dt*2*sqrt(pi)/self.db/self.Re*q/sqrt(A+DOLFIN_EPS)*self.v2*dx\
			- self.dt*(2*sqrt(A+DOLFIN_EPS)*(sqrt(pi)*self.f+sqrt(self.A0)*self.dfdr)\
				  -(A+DOLFIN_EPS)*self.dfdr)*self.drdx*self.v2*dx\
			- self.U_n[0]*self.v1*dx\
			- self.U_n[1]*self.v2*dx


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
