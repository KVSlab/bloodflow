import numpy as np
from fenics import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import conservation_solver as cs


class Artery(object):
	""" Represent an artery.
	:param Ru: Upstream radius
	:param Rd: Downstream radius
	:param L: Vessel length
	:param k: Vector containing k1, k2 and k3 from the relation Eh/R0
	:param Re: Reynolds number
	:param p0: Diastolic pressure
	"""
	def __init__(self, Ru, Rd, L, k, Re, p0):
		""" Construct artery.
		Add its intrinsic characteristics, not its numerical solution.
		"""
		self.Ru = Ru
		self.Rd = Rd
		self.L = L
		self.k = k
		self.Re = Re
		self.p0 = p0

	def define_problem(self, Nx, Nt, T):
		"""Define the FEniCS objects for the problem.
		:param Nx: Number of spatial steps
		:param Nt: Number of temporal steps
		:param T: Period of one cardiac cycle
		"""
		self.Nx = Nx
		self.Nt = Nt
		self.T = T
		self.mesh = IntervalMesh(Nx, 0, L)
		self.elV = FiniteElement('CG', mesh.ufl_cell(), 1)
		self.V = FunctionSpace(mesh, elV)
		self.V2 = FunctionSpace(mesh, elV*elV)
		
		self.U = TestFunction(V2)
		
		# Flux function
		self.F = Function(V2)
		
		# Source function
		self.S = Function(V2)

	def solve(self, q_in):
		"""Compute the solution U on the artery.
		:param q_in: Vector containing the inlet flow of the artery.
		"""
		cs.solve(self.mesh, self.V, self.V2, self.Nx, self.Nt, self.L, self.T,
			self.q_ins, self.A0, self.U, self.F, self.S, self.compute_A_out)


class Artery_Network(object):
	"""Describe a network of arteries.
	:param order: Number of arterial levels
	:param Ru: Upstream radii
	:param Rd: Downstream radii
	:param L: Vessel lengths
	:param k: Vectors containing k1, k2 and k3 from the relation Eh/R0
	:param Re: Reynolds number
	:param p0: Diastolic pressure
	"""
	def __init__(self, order, Ru, Rd, L, k, Re, p0):
		self.order = order
		self.arteries = [0] * (2**order-1)
		for i in range(len(arteries)):
			arteries[i] = Artery(Ru[i], Rd[i], L[i], k[i], Re, p0)

	def solve(self, Nx, Nt, T, q_in):
		"""
		:param Nx: Number of spatial steps
		:param Nt: Number of temporal steps
		:param T: Period of one cardiac cycle
		:param q_in: Vector of inlet flow for the first artery.
		"""
		for i in range(len(arteries)):
			arteries[i].define_problem(Nx, Nt, T)
			arteries[i].solve(q_in)

