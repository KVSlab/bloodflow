import sys
import numpy as np

from configparser import ConfigParser
import fenics as fn

sys.path.insert(0, 'src/')
from utils import *
from utils import is_near as near
from artery import Artery


def get_parameters(config_location):
	"""Read parameters for tests from file.
	:param config_location: Location of config file
	:return: The parameters needed for testing
	"""
	config = ConfigParser()
	config.read(config_location)

	# Constructor parameters
	root_vessel = config.getboolean('Parameters', 'root_vessel')
	end_vessel = config.getboolean('Parameters', 'end_vessel')
	rc = config.getfloat('Parameters', 'rc')
	qc = config.getfloat('Parameters', 'qc')
	Ru = config.getfloat('Parameters', 'Ru')
	Rd = config.getfloat('Parameters', 'Rd')
	L = config.getfloat('Parameters', 'L')
	k1 = config.getfloat('Parameters', 'k1')
	k2 = config.getfloat('Parameters', 'k2')
	k3 = config.getfloat('Parameters', 'k3')
	rho = config.getfloat('Parameters', 'rho')
	nu = config.getfloat('Parameters', 'nu')
	p0 = config.getfloat('Parameters', 'p0')

	# Geometry parameters
	Nt = config.getint('Geometry', 'Nt')
	Nx = config.getint('Geometry', 'Nx')
	T = config.getfloat('Geometry', 'T')
	N_cycles = config.getint('Geometry', 'N_cycles')

	# Solution parameters
	q0 = config.getfloat('Solution', 'q0')
	theta = config.getfloat('Solution', 'theta')
	
	# Adimensionalise parameters
	Ru, Rd, L, k1, k2, k3, Re, nu, p0, R1, R2, CT, q0, T = adimensionalise(
		rc, qc, Ru, Rd, L, k1, k2, k3, rho, nu, p0, R1, R2, CT, q0, T)
	
	return root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu,\
		   p0, Nt, Nx, T, N_cycles, q0, theta


def test_constructor(root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3,
					 rho, Re, nu, p0):
	"""Construct artery.
	Test correct assignment of parameters.
	:param boolean root_vessel: True if the artery is root-vessel (no parent)
	:param boolean end_vessel: True if the artery is end-vessel (no daughter)
	:param rc: Characteristic radius (length)
	:param qc: Characteristic flow
	:param Ru: Upstream radius
	:param Rd: Downstream radius
	:param L: Vessel length
	:param k1: First constant from the relation Eh/r0
	:param k2: Second constant from the relation Eh/r0
	:param k3: Third constant from the relation Eh/R0
	:param rho: Density of blood
	:param Re: Reynolds' number
	:param nu: Blood viscosity
	:param p0: Diastolic pressure
	:return: Intitialised artery object
	"""
	a = Artery(root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re,
			   nu, p0)
	
	assert(a.root_vessel == root_vessel)
	assert(a.end_vessel == end_vessel)
	assert(near(a.rc, rc))
	assert(near(a.qc, qc))
	assert(near(a.Ru, Ru))
	assert(near(a.Rd, Rd))
	assert(near(a.L, L))
	assert(near(a.k1, k1))
	assert(near(a.k2, k2))
	assert(near(a.k3, k3))
	assert(near(a.rho, rho))
	assert(near(a.Re, Re))
	assert(near(a.nu, nu))
	assert(near(a.p0, p0))
	
	return a


def test_define_geometry(a, Nx, Nt, T, N_cycles):
	"""Define geometry on artery.
	Test correct assignment of parameters.
	Test types of FEniCS objects.
	Test correct behaviour of FEniCS expressions.
	:param a: Artery on which define_geometry is to be called and tested
	:param Nx: Spatial refinement
	:param Nt: Number of temporal steps
	:param T: Duration of one cardiac cycle
	:param N_cycles: Number of cardiac cycles
	"""
	X = np.linspace(0, a.L, 100)
	
	# Higher tolerance due to error from FEniCS Expression-interpolation
	tol = 1.e-12
	
	a.define_geometry(Nx, Nt, T, N_cycles)
	
	assert(a.Nx == Nx)
	assert(a.Nt == Nt)
	assert(near(a.T, T))
	assert(a.N_cycles == N_cycles)
	assert(near(a.dx, a.L/Nx))
	assert(near(a.dt, T/Nt))
	assert(near(a.dex, a.L/Nx))
	assert(near(a.db, np.sqrt(a.nu*T/2/np.pi)))

	assert(isinstance(a.mesh, fn.IntervalMesh))
	assert(isinstance(a.elV, fn.FiniteElement))
	assert(isinstance(a.V, fn.FunctionSpace))
	assert(isinstance(a.V2, fn.FunctionSpace))

	for x in X:
		r0 = a.Ru*(a.Rd/a.Ru)**(x/a.L)
		A0 = np.pi*(a.Ru*(a.Rd/a.Ru)**(x/a.L))**2
		f = 4/3*(a.k1*np.exp(a.k2*a.Ru*(a.Rd/a.Ru)**(x/a.L))+a.k3)
		dfdr = 4/3*a.k1*a.k2*np.exp(a.k2*a.Ru*(a.Rd/a.Ru)**(x/a.L))
		drdx = np.log(a.Rd/a.Ru)/a.L*a.Ru*(a.Rd/a.Ru)**(x/a.L)
		assert(near(a.r0(x), r0))
		assert(near(a.A0(x), A0))
		assert(near(a.f(x), f))
		assert(near(a.dfdr(x), dfdr))
		assert(near(a.drdx(x), drdx))


def test_define_solution(a, q0, theta, bc_tol=1.e-14):
	"""Define solution on artery.
	Test correct assignment of parameters.
	Test types of FEniCS objects.
	Test correct behaviour of boundaries.
	Test variational form.
	:param a: Artery on which define_solution is to be called and tested
	:param q0: Initial flow
	:param theta: Crank-Nicolson parameter
	:param bc_tol: Inlet and outlet boundary thickness (tolerance)
	"""
	X = np.linspace(0, a.L, 100)
	
	a.define_solution(q0, theta, bc_tol)
	
	assert(near(a.q0, q0))
	assert(near(a.theta, theta))
	
	assert(type(a.U) == fn.Function)
	assert(type(a.Un) == fn.Function)
	assert(type(a.pn) == fn.Function)
	
	for x in X:
		assert(near(a.Un(x)[0], a.A0(x)))
		assert(near(a.Un(x)[1], q0))


def test_solve(a):
	"""Solve equation on artery for one time-step.
	:param a: Artery on which the equation is to be solved
	"""
	a.solve()


def test_update_solution(a):
	"""Update solution.
	Test equality beween Un and U.
	"""
	a.update_solution()
	Un = a.Un.vector().get_local()
	U = a.U.vector().get_local()
	for i in range(len(Un)):
		assert(near(Un[i], U[i]))


def test_update_pressure(a):
	"""Update pressure.
	Test correct behaviour of pressure function.
	"""
	X = np.linspace(0, a.L, 100)
	a.update_pressure()
	reltol = 1.e-12
	
	for x in X:
		p = a.p0 + a.f(x)*(1-np.sqrt(a.A0(x)/a.Un(x)[0]))
		assert(near(a.pn(x), p, reltol=reltol))


def test_compute_pressure(a):
	"""Test correct value of computed pressure.
	"""
	for x in np.linspace(0, a.L, 100):
		f, A0, A = a.f(x), a.A0(x), a.Un(x)[0]
		p = a.compute_pressure(f, A0, A)
		assert(near(p, a.p0+f*(1-np.sqrt(A0/A))))


def test_compute_outlet_pressure(a):
	"""Test correct value of outlet pressure.
	"""
	A = a.Un(a.L)[0]
	p = a.compute_outlet_pressure(A)
	assert(near(p, a.p0+a.f(a.L)*(1-np.sqrt(a.A0(a.L)/A))))


def test_CFL_term(a):
	"""Test correct value of CFL-term.
	"""
	for x in [0, a.L]:
		A, q = a.Un(x)
		CFL = 1/np.abs(q/A+np.sqrt(a.f(x)/2/a.rho*np.sqrt(a.A0(x)/A)))
		assert(near(a.CFL_term(x, A, q), CFL))


def test_check_CFL(a):
	"""Test CFL-condition-checking.
	"""
	margin = 1.e-10
	
	for x in [0, a.L]:
		A, q = a.Un(x)
		M = a.dt/a.CFL_term(x, A, q)
		a.dex = (1-margin)*M
		assert(not a.check_CFL(x, A, q))
		a.dex = (1+margin)*M
		assert(a.check_CFL(x, A, q))


def test_adjust_dex(a):
	"""Test correct adjustment of dex.
	"""
	for margin in [0.1, 0.05, 1.e-4, 1.e-8, 1.e-10]:
		for x in [0, a.L]:
			A, q = a.Un(x)
			a.adjust_dex(x, A, q, margin)
			assert(near(a.dex, (1+margin)*a.dt/a.CFL_term(x, A, q)))
			assert(a.check_CFL(x, A, q))


def test_artery(config_location):
	"""Test artery class.
	:param config_location: Location of config-file with test-parameters
	"""
	root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0,\
		Nt,	Nx, T, N_cycles, q0, theta = get_parameters(config_location)
	
	a = test_constructor(root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3,
						 rho, Re, nu, p0)
						 
	test_define_geometry(a, Nx, Nt, T, N_cycles)
	
	test_define_solution(a, q0, theta)
	
	test_solve(a)
	
	test_update_solution(a)
	
	test_update_pressure(a)

	test_compute_pressure(a)
	
	test_compute_outlet_pressure(a)
	
	test_CFL_term(a)
	
	test_check_CFL(a)
	
	test_adjust_dex(a)

if __name__ == '__main__':
	test_artery(sys.argv[1])
	
