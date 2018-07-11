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
