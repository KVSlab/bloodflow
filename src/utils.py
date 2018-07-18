import numpy as np

from fenics import *
from mshr import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from configparser import ConfigParser


# Pressure unit converting functions ('unit' = g cm-1 s-2)
def unit_to_mmHg(p):
	return 76/101325*p
	
def mmHg_to_unit(p):
	return 101325/76*p


def adimensionalise(rc, qc, Ru, Rd, L, k1, k2, k3, rho, nu, p0, q_ins, T):
	"""Make quantities independent of units.
	:param rc: Characteristic length (radius)
	:param qc: Characteristic flow
	:param Ru: Radius
	:param Rd: Radius
	:param L: Length
	:param k1: Pressure
	:param k2: Inverse length
	:param k3: Pressure
	:param rho: Density
	:param nu: Viscosity
	:param p0: Pressure
	:param q_ins: Flow
	:param T: Time
	:return: Adimensionalised quantities, including Reynolds' number
	"""
	Ru = Ru/rc
	Rd = Rd/rc
	L = L/rc
	k1 = k1*rc**4/rho/qc**2
	k2 = k2*rc
	k3 = k3*rc**4/rho/qc**2
	Re = qc/nu/rc
	nu = nu*rc/qc
	p0 = p0*rc**4/rho/qc**2
	q_ins = q_ins/qc
	T = T*qc/rc**3
	return Ru, Rd, L, k1, k2, k3, Re, nu, p0, q_ins, T


def redimensionalise(rc, qc, rho, x, nature):
	"""Give a quantity its unit back
	:param rc: Characteristic radius
	:param qc: Characteristic flow
	:param rho: Density
	:param x: Quantity to be redimensionalised
	:param string nature: Nature of the quantity to be redimensionalised.
	"""
	if nature == 'time':
		x = x*rc**3/qc
	elif nature == 'area':
		x = x*rc**2
	elif nature == 'flow':
		x = x*qc
	elif nature == 'pressure':
		x = x*rho*qc**2/rc**4
	return x


def read_output(filename):
	config = ConfigParser()
	config.read(filename)
	
	order = config.getint('data', 'order')
	Nx = config.getint('data', 'Nx')
	Nt = config.getint('data', 'Nt')
	T = config.getfloat('data', 'T')
	L = [float(f) for f in config.get('data', 'L').split(',')]
	rc = config.getfloat('data', 'rc')
	qc = config.getfloat('data', 'qc')
	rho = config.getfloat('data', 'rho')
	mesh_location = config.get('data', 'mesh_location')
	locations = config.get('data', 'locations').split(',')
	names = config.get('data', 'names').split(',')
	
	return order, Nx, Nt, T, L, rc, qc, rho, mesh_location, locations, names


def XDMF_to_matrix(Nx, Nt, mesh_location, location, name):
	M = np.zeros([Nx, Nt])
	f = XDMFFile(location)
	mesh = Mesh(mesh_location)
	V = FunctionSpace(mesh, 'CG', 1)
	u = Function(V)
	for n in range(Nt):
		f.read_checkpoint(u, name, n)
		fig = plt.figure()
		plot(u)
		if name == 'area':
			plt.ylim(0, 0.7)
		else:
			plt.ylim(0, 27)
		plt.savefig('output/plot/%s_%i' % (name, n))
		plt.close(fig)
		M[:, n] = u.vector().get_local()[:-1][::-1]
	f.close()
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
