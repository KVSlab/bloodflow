__author__ = 'Syver Døving Agdestein'


import numpy as np
from fenics import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from artery import Artery
import conservation_solver as cs


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
	def __init__(self, order, Ru, Rd, L, k1, k2, k3, Re, p0):
		self.order = order
		self.arteries = [0] * (2**order-1)
		for i in range(len(arteries)):
			arteries[i] = Artery(Ru[i], Rd[i], L[i], k1[i], k2[i], k3[i], Re, p0)
	
	def jacobian():
		"""Compute the jacobian matrix for the system of equations <system>.
		:param x:
		"""
		J = np.zeros([18, 18])
		
		# Equation (20)
		J[1, 0] = 2
		J[2, 0] = -1
		J[4, 1] = 2
		J[5, 1] = -1
		J[7, 2] = 2
		J[8, 2] = -1
		
		# Equation (21)
		J[10, 3] = 2
		J[11, 3] = -1
		J[13, 4] = 2
		J[14, 4] = -1
		J[16, 5] = 2
		J[17, 5] = -1
		
		# Equation (22)
		J[1, 6] = 1
		J[4, 6] = -1
		J[7, 6] = -1
		J[2, 7] = 1
		J[5, 7] = -1
		J[8, 7] = -1
		
		# Equation (23)
		J[10, 8] = 
		J[13, 8] = 
		J[10, 9] = 
		J[16, 9] = 
		J[9, 10] = 
		J[12, 10] = 
		J[9, 11] = 
		J[15, 11] = 
		
		# Equation (26)
		J[, 12] = 
		J[, 13] = 
		J[, 14] = 
		
		# Equation (27)
		J[, 15] = 
		J[, 16] = 
		J[, 17] = 
		return J

	def solve(self, Nx, Nt, T, N_cycles, q_in):
		"""Solve the equation on the entire arterial network.
		:param Nx: Number of spatial steps
		:param Nt: Number of temporal steps
		:param T: Period of one cardiac cycle
		:param N_cycles: Number of cardiac cycles
		:param q_in: Vector of inlet flow for the first artery.
		"""
		for i in range(len(arteries)):
			arteries[i].define_geometry(Nx, Nt, T, N_cycles)
			cs.solve(arteries[i], q_ins)
