import numpy as np

from fenics import *
from mshr import *
import matplotlib.pyplot as plt
from configparser import ConfigParser


# Pressure unit converting functions ('unit' = g cm-1 s-2)
def unit_to_mmHg(p):
	return 76/101325*p
	
def mmHg_to_unit(p):
	return 101325/76*p


def read_output(filename):
	config = ConfigParser()
	config.read(filename)
	
	order = config.getint('data', 'order')
	Nx = config.getint('data', 'Nx')
	Nt = config.getint('data', 'Nt')
	T = config.getfloat('data', 'T')
	L = [float(f) for f in config.get('data', 'L').split(',')]
	mesh_location = config.get('data', 'mesh_location')
	locations = config.get('data', 'locations').split(',')
	names = config.get('data', 'names').split(',')
	
	return order, Nx, Nt, T, L, mesh_location, locations, names


def XDMF_to_matrix(Nx, Nt, mesh_location, location, name):
	M = np.zeros([Nx+1, Nt])
	f = XDMFFile(location)
	mesh = Mesh(mesh_location)
	V = FunctionSpace(mesh, 'CG', 1)
	u = Function(V)
	for n in range(Nt):
		f.read_checkpoint(u, name, n)
		M[:, n] = u.vector().array()[::-1]
	return M


def plot_matrix(t, x, M, label, output):
	"""Create a plot of a matrix.
	Store it in a given location.
	:param t: Vector containing time values
	:param x: Vector containing space values
	:param M: Matrix representing function to be plotted, with dimension t*x
	:param string label: Name of the function
	:param string output: Location (with filename) of output file
	"""
	T, X = np.meshgrid(t, x)
	fig = plt.figure(figsize=(8, 6))
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(T, X, M, rstride=1, cstride=1,  cmap='viridis',
						   linewidth=0, antialiased=False)
	ax.set_xlabel('t')
	ax.set_ylabel('x')
	ax.set_zlabel(label)
	ax.set_ylim(min(x), max(x))
	ax.set_xlim(min(t), max(t))
	plt.savefig(output)
