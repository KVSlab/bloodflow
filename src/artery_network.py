import numpy as np
from fenics import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import conservation_solver as cs

class Artery(object):
	"""
	:param Ru: Upstream radius
	:param Rd: Downstream radius
	:param L: Vessel length
	:param k: Vector containing k1, k2 and k3 from the relation Eh/R0
	:param Re: Reynolds number
	:param p0: Diastolic pressure
	"""
	def __init__(self, Ru, Rd, L, k, Re, p0):
		""" Constructor for the artery.
		Add its fundamental characteristics.
		"""
		self.Ru = Ru
		self.Rd = Rd
		self.L = L
		self.k = k
		self.Re = Re
		self.p0 = p0
		
	def define_problem(self, Nx, Nt, T):
		"""
		"""
		self.mesh = IntervalMesh(Nx, 0, L)
		self.elV = FiniteElement('CG', mesh.ufl_cell(), 1)
		self.V = FunctionSpace(mesh, elV)
		self.V2 = FunctionSpace(mesh, elV*elV)
		
	
	def solve(self, Nx, Nt, T, q_in):
		"""
		"""
		None
		
class ArteryNetwork(object):
	"""
	"""
	def __init__(self, order, Ru, Rd, L, k, Re, p0):
		self.order = order
		self.arteries = [0] * (order**2-1)
		for i in range(len(arteries)):
			arteries[i] = Artery(Ru[i], Rd[i], L[i], k[i], Re, p0)
	
	def solve(self, Nx, Nt, T, q_in):
		for i in range(len(arteries)):
			arteries[i].define_problem(Nx, Nt, T)
			arteries[i].solve(q_in)
			
