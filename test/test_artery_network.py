import sys

import numpy as np
from configparser import ConfigParser

sys.path.insert(0, 'src/')
sys.path.insert(0, 'test/')

import test_artery as ta
from test_artery import near
from artery_network import Artery_Network


def get_parameters(config_location):
	"""Read parameters for tests from file.
	:param config_location: Location of config file
	:return: The parameters needed for testing
	"""
	config = ConfigParser()
	config.read(config_location)

	# Constructor parameters
	order = config.getint('Parameters', 'order')
	rc = config.getfloat('Parameters', 'rc')
	qc = config.getfloat('Parameters', 'qc')
	Ru = np.array([float(f) for f in config.get('Parameters', 'Ru').split(',')])
	Rd = np.array([float(f) for f in config.get('Parameters', 'Rd').split(',')])
	L = np.array([float(f) for f in config.get('Parameters', 'L').split(',')])
	k1 = config.getfloat('Parameters', 'k1')
	k2 = config.getfloat('Parameters', 'k2')
	k3 = config.getfloat('Parameters', 'k3')
	rho = config.getfloat('Parameters', 'rho')
	Re = config.getfloat('Parameters', 'Re')
	nu = config.getfloat('Parameters', 'nu')
	p0 = config.getfloat('Parameters', 'p0')
	R1 = config.getfloat('Parameters', 'R1')
	R2 = config.getfloat('Parameters', 'R2')
	CT = config.getfloat('Parameters', 'CT')
	

	# Geometry parameters
	Nt = config.getint('Geometry', 'Nt')
	Nx = config.getint('Geometry', 'Nx')
	T = config.getfloat('Geometry', 'T')
	N_cycles = config.getint('Geometry', 'N_cycles')

	# Solution parameters
	output_location = config.get('Solution', 'output_location')
	theta = config.getfloat('Solution', 'theta')
	Nt_store = config.getint('Solution', 'Nt_store')
	N_cycles_store = config.getint('Solution', 'N_cycles_store')
	store_area = config.getint('Solution', 'store_area')
	store_pressure = config.getint('Solution', 'store_pressure')
	q0 = config.getfloat('Solution', 'q0')
	
	return order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT, Nt, Nx, T, N_cycles, output_location, theta, Nt_store, N_cycles_store, store_area, store_pressure, q0




def test_constructor(order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT):
	"""Construct artery network.
	Test correct assignment of parameters.
	Test correct structure of network.
	:param order: Number of arterial levels
	:param rc: Characteristic radius (length)
	:param qc: Characteristic flow
	:param Ru: Upstream radii
	:param Rd: Downstream radii
	:param L: Vessel lengths
	:param k1: First constant from the relation Eh/r0
	:param k2: Second constant from the relation Eh/r0
	:param k3: Third constant from the relation Eh/R0
	:param rho: Density of blood
	:param Re: Reynolds' number
	:param nu: Viscosity of blood
	:param p0: Diastolic pressure
	:param R1: First resistance from Windkessel model
	:param R2: Second resistance from Windkessel model
	:param CT: Compliance from Windkessel model
	:return: Constructed artery network
	"""
	an = Artery_Network(order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT)
	
	assert(an.order == order)
	assert(len(an.arteries) == 2**order-1)
	assert(an.range_arteries == range(2**order-1))
	assert(an.range_parent_arteries == range(2**(order-1)-1))
	assert(an.range_daughter_arteries == range(1, 2**order-1))
	assert(an.range_end_arteries == range(2**(order-1)-1, 2**order-1))
	assert(an.rc == rc)
	assert(an.qc == qc)
	assert(an.rho == rho)
	assert(an.R1 == R1)
	assert(an.R2 == R2)
	assert(an.CT == CT)
	
	for i, artery in enumerate(an.arteries):
		assert(artery.root_vessel if i==0 else not artery.root_vessel)
		assert(artery.end_vessel if i in an.range_end_arteries\
			   else not artery.end_vessel)
	
	return(an)


def test_define_geometry(an, Nx, Nt, T, N_cycles):
	"""Define geometry on artery network.
	Test correct assignment of parameters.
	"""
	an.define_geometry(Nx, Nt, T, N_cycles)

	assert(an.Nx == Nx)
	assert(an.Nt == Nt)
	assert(near(an.T, T))
	assert(an.N_cycles == N_cycles)


def test_define_solution(an, output_location, q0, theta):
	"""Define solution on artery network.
	Test correct assignment of parameters.
	"""
	an.define_solution(output_location, q0, theta)
	
	assert(an.output_location == output_location)
	assert(an.theta == theta)


def test_daughter_arteries(an):
	"""Test correct finding of daughter vessels.
	"""
	for ip in an.range_parent_arteries:
		i1, i2 = an.daughter_arteries(ip)
		assert(i1 == 2*ip+1)
		assert(i2 == 2*ip+2)


def test_parent_artery(an):
	"""Test correct indices for parent vessels.
	"""
	for i in an.range_daughter_arteries:
		ip = an.parent_artery(i)
		if i % 2 == 1:
			assert(ip == (i-1)//2)
		else:
			assert(ip == (i-2)//2)


def test_sister_artery(an):
	"""Test correct indices for sister vessel.
	"""
	for i in an.range_daughter_arteries:
		s = an.sister_artery(i)
		if i % 2 == 1:
			assert(s == i+1)
		else:
			assert(s == i-1)


def test_flux(an):
	"""Test correct behaviour of flux function.
	"""
	for a in an.arteries:
		for x in np.linspace(0, a.L, 100):
			U = a.Un(x)
			anflux = an.flux(a, U, x)
			F1 = U[1]
			F2 = U[1]**2 + a.f(x)*np.sqrt(a.A0(x)*U[0])
			assert(near(anflux[0], F1))
			assert(near(anflux[1], F2))


def test_source(an):
	"""Test correct behaviour of source function.
	"""
	for a in an.arteries:
		for x in np.linspace(0, a.L, 100):
			U = a.Un(x)
			ansource = an.source(a, U, x)
			S2 = - 2*np.sqrt(np.pi)/a.db/a.Re*U[1]/np.sqrt(U[0]) + (2*np.sqrt(U[0])*(np.sqrt(np.pi)*a.f(x) + np.sqrt(a.A0(x))*a.dfdr(x)) - U[0]*a.dfdr(x))*a.drdx(x)
			assert(near(ansource[0], 0))
			assert(near(ansource[1], S2))


def test_compute_U_half(an):
	"""
	"""


def test_compute_A_out(an):
	"""
	"""


def test_initial_x(an):
	"""
	"""


def test_define_x(an):
	"""
	"""


def test_problem_function(an):
	"""
	"""


def test_jacobian(an):
	"""Test that the analytical expression for the jacobian matrix is close to a numerically computed jacobian.
	"""
	tol = 1.e-4
	reltol = 1.e-4
	for ip in an.range_parent_arteries:
		i1, i2 = an.daughter_arteries(ip)
		p, d1, d2 = an.arteries[ip], an.arteries[i1], an.arteries[i2]
		x = np.ones(18)
		h = 1.e-10
		he = np.zeros(18)
		J = an.jacobian(p, d1, d2, x)
		F = an.problem_function(p, d1, d2, x)
		for j in range(18):
			he[j] += h
			Fph = an.problem_function(p, d1, d2, x+he)
			dF = (Fph-F)/h
			for i in range(18):
				assert(near(J[i, j], dF[i], tol, reltol))
			he[j] = 0

def test_newton(an):
	"""
	"""


def test_adjust_bifurcation_step(an):
	"""
	"""


def test_set_inner_bc(an):
	"""
	"""


def test_set_bcs(an):
	"""
	"""


def test_dump_metadata(an):
	"""
	"""


def test_solve(an):
	"""
	"""


def test_artery_network(config_location):
	"""Test correct functionality of the methods from the artery_network class.
	"""
	order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT, Nt, Nx, T, N_cycles, output_location, theta, Nt_store, N_cycles_store, store_area, store_pressure, q0 = get_parameters(config_location)
	
	an = test_constructor(order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT)

	test_daughter_arteries(an)

	test_parent_artery(an)

	test_sister_artery(an)

	test_define_geometry(an, Nx, Nt, T, N_cycles)

	test_define_solution(an, output_location, q0, theta)

	test_flux(an)

	test_source(an)

	test_compute_U_half(an)

	test_compute_A_out(an)

	test_problem_function(an)

	test_jacobian(an)

	test_newton(an)

	test_adjust_bifurcation_step(an)

	test_set_inner_bc(an)

	test_set_bcs(an)

	test_initial_x(an)

	test_define_x(an)

	test_dump_metadata(an)

	test_solve(an)


if __name__ == '__main__':
	test_artery_network(sys.argv[1])
