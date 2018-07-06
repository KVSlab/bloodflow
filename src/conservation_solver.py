import numpy as np
from fenics import *

def solve(mesh, V, V2, Nx, Nt, L, T, q_in, A0, F, S, A_out):
	""" Return the solution of dU/dt + dF/dx = S.
	:param mesh: 1D mesh of the interval [0, L] with Nx points
	:param V: Scalar function space on mesh
	:param V2: Vector function space on mesh
	:param Nx: Number of spatial discretisation
	:param 
	:param q_in: Inlet boundary condition for the second component of U
	:param A0: Initial value of the first component of U
	:param F: Flux function
	:param S: Source function
	:param A_out: Function that computes the next outlet boundary value of A
	"""
	
