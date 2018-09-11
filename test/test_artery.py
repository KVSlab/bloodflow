import sys
import numpy as np

import pytest

from configparser import ConfigParser
import fenics as fn

sys.path.insert(0, 'src/')
from utils import *
from utils import is_near as near
from artery import Artery


def test_constructor(artery, param):
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

	root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0,\
		Nt,	Nx, T, N_cycles, q0, theta = param

	assert(artery.root_vessel == root_vessel)
	assert(artery.end_vessel == end_vessel)
	assert(near(artery.rc, rc))
	assert(near(artery.qc, qc))
	assert(near(artery.Ru, Ru))
	assert(near(artery.Rd, Rd))
	assert(near(artery.L, L))
	assert(near(artery.k1, k1))
	assert(near(artery.k2, k2))
	assert(near(artery.k3, k3))
	assert(near(artery.rho, rho))
	assert(near(artery.Re, Re))
	assert(near(artery.nu, nu))
	assert(near(artery.p0, p0))


def test_define_geometry(artery, param):
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
	root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0,\
		Nt,	Nx, T, N_cycles, q0, theta = param

	X = np.linspace(0, artery.L, 100)

	# Higher tolerance due to error from FEniCS Expression-interpolation
	tol = 1.e-12

	artery.define_geometry(Nx, Nt, T, N_cycles)

	assert(artery.Nx == Nx)
	assert(artery.Nt == Nt)
	assert(near(artery.T, T))
	assert(artery.N_cycles == N_cycles)
	assert(near(artery.dx, artery.L/Nx))
	assert(near(artery.dt, T/Nt))
	assert(near(artery.dex, artery.L/Nx))
	assert(near(artery.db, np.sqrt(artery.nu*T/2/np.pi)))

	assert(isinstance(artery.mesh, fn.IntervalMesh))
	assert(isinstance(artery.elV, fn.FiniteElement))
	assert(isinstance(artery.V, fn.FunctionSpace))
	assert(isinstance(artery.V2, fn.FunctionSpace))

	a = artery

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


def test_define_solution(artery, param, bc_tol=1.e-14):
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
	root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0,\
		Nt,	Nx, T, N_cycles, q0, theta = param

	X = np.linspace(0, artery.L, 100)

	artery.define_geometry(Nx, Nt, T, N_cycles)
	artery.define_solution(q0, theta, bc_tol)

	assert(near(artery.q0, q0))
	assert(near(artery.theta, theta))

	assert(type(artery.U) == fn.Function)
	assert(type(artery.Un) == fn.Function)
	assert(type(artery.pn) == fn.Function)

	for x in X:
		assert(near(artery.Un(x)[0], artery.A0(x)))
		assert(near(artery.Un(x)[1], q0))


def test_solve(artery_def):
	"""Solve equation on artery for one time-step.
	:param a: Artery on which the equation is to be solved
	"""
	artery_def.solve()


def test_update_solution(artery_def):
	"""Update solution.
	Test equality beween Un and U.
	"""
	a = artery_def
	a.update_solution()
	Un = a.Un.vector().get_local()
	U = a.U.vector().get_local()
	for i in range(len(Un)):
		assert(near(Un[i], U[i]))


def test_update_pressure(artery_def):
	"""Update pressure.
	Test correct behaviour of pressure function.
	"""
	a = artery_def
	X = np.linspace(0, a.L, 100)
	a.update_pressure()
	reltol = 1.e-12

	for x in X:
		p = a.p0 + a.f(x)*(1-np.sqrt(a.A0(x)/a.Un(x)[0]))
		assert(near(a.pn(x), p, reltol=reltol))


def test_compute_pressure(artery_def):
	"""Test correct value of computed pressure.
	"""
	a = artery_def
	for x in np.linspace(0, a.L, 100):
		f, A0, A = a.f(x), a.A0(x), a.Un(x)[0]
		p = a.compute_pressure(f, A0, A)
		assert(near(p, a.p0+f*(1-np.sqrt(A0/A))))


def test_compute_outlet_pressure(artery_def):
	"""Test correct value of outlet pressure.
	"""
	a = artery_def
	A = a.Un(a.L)[0]
	p = a.compute_outlet_pressure(A)
	assert(near(p, a.p0+a.f(a.L)*(1-np.sqrt(a.A0(a.L)/A))))


def test_CFL_term(artery_def):
	"""Test correct value of CFL-term.
	"""
	a = artery_def
	for x in [0, a.L]:
		A, q = a.Un(x)
		CFL = 1/np.abs(q/A+np.sqrt(a.f(x)/2/a.rho*np.sqrt(a.A0(x)/A)))
		assert(near(a.CFL_term(x, A, q), CFL))


def test_check_CFL(artery_def):
	"""Test CFL-condition-checking.
	"""
	margin = 1.e-10

	a = artery_def

	for x in [0, a.L]:
		A, q = a.Un(x)
		M = a.dt/a.CFL_term(x, A, q)
		a.dex = (1-margin)*M
		assert(not a.check_CFL(x, A, q))
		a.dex = (1+margin)*M
		assert(a.check_CFL(x, A, q))


def test_adjust_dex(artery_def):
	"""Test correct adjustment of dex.
	"""
	a = artery_def
	for margin in [0.1, 0.05, 1.e-4, 1.e-8, 1.e-10]:
		for x in [0, a.L]:
			A, q = a.Un(x)
			a.adjust_dex(x, A, q, margin)
			assert(near(a.dex, (1+margin)*a.dt/a.CFL_term(x, A, q)))
			assert(a.check_CFL(x, A, q))


@pytest.fixture
def config_location():
    fconfig = '/tmp/config.cfg'
    f = open(fconfig, 'w')
    f.write(FCONF)
    f.close()
    return fconfig


@pytest.fixture
def param(config_location):
	config = ConfigParser()
	config.read(config_location)

	# Constructor parameters
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
	N_cycles = config.getint('Geometry', 'N_cycles')

	root_vessel = True
	end_vessel = True
	T = 1.0
	q0 = 2.0
	theta = 0.51

	# Adimensionalise parameters
	R1 = R2 = CT = 0
	Ru, Rd, L, k1, k2, k3, Re, nu, p0, R1, R2, CT, q0, T = \
		adimensionalise_parameters(rc, qc, Ru, Rd, L, k1, k2, k3, rho, nu, p0,
								   R1, R2, CT, q0, T)

	return root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu,\
		   p0, Nt, Nx, T, N_cycles, q0, theta


@pytest.fixture
def artery(param):
	root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0,\
		Nt,	Nx, T, N_cycles, q0, theta = param
	artery = Artery(root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re,
			   nu, p0)
	return artery


@pytest.fixture
def artery_def(artery, param):
	root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0,\
		Nt,	Nx, T, N_cycles, q0, theta = param
	artery_def = artery
	artery.define_geometry(Nx, Nt, T, N_cycles)
	artery.define_solution(q0, theta, bc_tol=1.e-14)
	return artery_def


FCONF = """
[Parameters]
order = 2
rc = 1.0
qc = 10.0
Ru = 0.37
Rd = 0.37
L = 20.8
k1 = 2.0e7
k2 = -22.53
k3 = 8.65e5
rho = 1.06
nu = 0.046
p0 = 119990.0
R1 = 25300.0
R2 = 13900.0
CT = 1.3384e-6

[Geometry]
Nx = 200
Nt = 2000
N_cycles = 1

[Solution]
inlet_flow_location = data/example_inlet.csv
output_location = output/1cycle
theta = 0.55
Nt_store = 200
N_cycles_store = 1
store_area = 1
store_pressure = 1
"""
