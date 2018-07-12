__author__ = 'Syver Døving Agdestein'


import numpy as np
import numpy.linalg as npl
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
			arteries[i] = Artery(Ru[i], Rd[i], L[i], k1[i],
								 k2[i], k3[i], Re, p0)
	
	
	def problem_function(p, d1, d2, x):
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
		
		y = np.zeros(18)
		
		# Entries from equation (20)
		y[0] = 2*x[1]
		y[1] = 
		y[2] = 
		
		# Entries from equation (21)
		y[3] = 
		y[4] = 
		y[5] = 
		
		# Entries from equation (22)
		y[6] = 
		y[7] = 
		
		# Entries from equation (23)
		y[8] = 
		y[9] = 
		y[10] = 
		y[11] =
		
		# Entries from equation (26) 
		y[12] = 
		y[13] = 
		y[14] = 
		
		# Entries from equation (27)
		y[15] = 
		y[16] = 
		y[17] = 
		
		return y


	def jacobian(p, d1, d2, x):
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
		J[11, 12] = p.dt/p.dex*(-pow(x[2]/x[11], 2) + fp/2*np.sqrt(A0p/x[11]))
				  - p.dt/2*(rpi/dbp/Rep*x[2]/pow(x[11], 3/2)
						   + (1/np.sqrt(x[11])*(rpi*fp+np.sqrt(A0p)*dfdrp)
											   - dfdrp)*drdxp)
		J[3, 13] = 1
		J[5, 13] = -d1.dt/d1.dex*2*x[5]/x[14] + d1.dt*rpi/db1/Re1/np.sqrt(x[14])
		J[14, 13] = d1.dt/d1.dex*(pow(x[5]/x[14], 2) - f1/2*np.sqrt(A01/x[14]))
				  - d1.dt/2*(rpi/db1/Re1*x[5]/pow(x[14], 3/2)
							+ (1/np.sqrt(x[14])*(rpi*f1+np.sqrt(A01)*dfdr1)
												- dfdr1)*drdx1)
		J[6, 14] = 1
		J[8, 14] = -d2.dt/d2.dex*2*x[8]/x[17] + d2.dt*rpi/db2/Re2/np.sqrt(x[17])
		J[17, 14] = d2.dt/d2.dex*(pow(x[8]/x[17], 2) - f2/2*np.sqrt(A02/x[17]))
				  - d2.dt/2*(rpi/db2/Re2*x[8]/pow(x[17], 3/2)
							+ (1/np.sqrt(x[17])*(rpi*f2+np.sqrt(A02)*dfdr2)
												- dfdr2)*drdx2)
		
		# Entries from equation (27)
		J[2, 15] = p.dt/p.dex
		J[11, 15] = 1
		J[5, 16] = d1.dt/d1.dex
		J[14, 16] = 1
		J[8, 17] = d1.dt/d1.dex
		J[17, 17] = 1
		
		return J
		

	def newton(p, d1, d2, k_max=1000, tol=1.e-14):
		""" Compute the boundary conditions.
		:param p: Parent artery
		:param d1: First daughter vessel
		:param d2: Second daughter vessel
		:param x: Current point, an 18-dimensional vector
		:return: Solution to the system of equations
		"""
		x = np.zeros(18)
		for k in range(k_max):
			x_old = x
			x -= npl.solve(jacobian(x), problem_function(x))
			if npl.norm(x-x_old) < tol:
				break
		return x
	
	def compute_inner_bc(p, d1, d2):
		""" Compute the inter-arterial boundary conditions for one bifurcation.
		:param p: Parent artery
		:param d1: First daughter vessel
		:param d2: Second daughter vessel
		:return: Area and flow for the three vessels
		"""
		x = newton(p, d1, d2)
		
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
