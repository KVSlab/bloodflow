import sys

import fenics as fn

sys.path.insert(0, 'src/')
from artery import Artery

def near(a, b, tol=1.e-15):
	"""Check equality between two floats with a certain tolerance.
	:param a: First number
	:param b: Second number
	:param tol: Tolerance for equality
	:return: True if the two numbers are to be considered equal
	"""
	return np.abs(a-b) < tol


def get_parameters(config_location):
	"""Read parameters for tests from file.
	:param config_location: Location of config file
	:return: The parameters needed for testing
	"""
	config = configparser.ConfigParser()
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
	Re = config.getfloat('Parameters', 'Re')
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
	
	return root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, nu,
		   p0, Nt, Nx, T, N_cycles, q0, theta


def test_constructor(root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3,
					 rho, nu, p0):
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
	a = Artery(rc, qc, Ru, Rd, L, k1, k2, k3, rho, nu, p0)
	
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
	"""
	X = np.linspace(0, L, 100)
	
	# Higher tolerance due to error from FEniCS Expression-interpolation
	tol = 1.e-12
	
	a.define_geometry(Nx, Nt, T, N_cycles)
	
	assert(a.Nx == Nx)
	assert(a.Nt == Nt)
	assert(near(a.T, T))
	assert(a.N_cycles == N_cycles)
	assert(near(a.dx, L/Nx))
	assert(near(a.dt, T/Nt))
	assert(near(a.dex, L/Nx))
	assert(near(a.db, np.sqrt(nu*T/2/np.pi)))
	
	assert(type(a.mesh) == ''))
	assert(type(a.elV) == ''))
	assert(type(a.V) == ''))
	assert(type(a.V2) == ''))
	
	for x in X:
		assert(near(a.r0(x), Ru*(Rd/Ru)**(x/L)))
		assert(near(a.A0(x), np.pi*(Ru*(Rd/Ru)**(x/L))**2))
		assert(near(a.f(x), 4/3*(k1*np.exp(k2*Ru*(Rd/Ru)**(x/L))+k3)))
		assert(near(a.dfdr(x), 4/3*k1*k2*np.exp(k2*Ru*(Rd/Ru)**(x/L))))
		assert(near(a.drdx(x), np.ln(Rd/Ru)/L*Ru*(Rd/Ru)**(x[0]/L)))


def test_define_solution(a, q0, theta):
	"""Define solution on artery.
	Test correct assignment of parameters.
	Test types of FEniCS objects.
	Test correct behaviour of boundaries.
	Test variational form.
	:param a: Artery on which define_solution is to be called and tested
	"""
	bc_tol = 1.e-14
	a.define_solution(q0, theta, bc_tol)
	
	assert(near(a.q0, q0))
	assert(near(a.theta, theta))
	
	assert(type(a.U) == type(fn.Function(a.V2)))
	assert(type(a.v1) == '')
	assert(type(a.v2) == '')
	assert(type(a.Un) == '')
	assert(type(a.pn == '')
	
	assert(near(a.Un(x)[0], np.pi*(Ru*(Rd/Ru)**(x/L))**2
	assert(near(a.Un(x)[1], q0))


def test_solve():
def test_update_solution():
def test_update_pressure():
def test_compute_pressure():
def test_compute_outlet_pressure():
def test_CFL_term():
def test_check_CFL():
def test_adjust_dex()


def test_artery(config_location)
	
	root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3, rho, nu, p0, Nt,
		Nx, T, N_cycles, q0, theta = get_parameters(config_location)
	
	a = test_constructor(root_vessel, end_vessel, rc, qc, Ru, Rd, L, k1, k2, k3,
						 rho, nu, p0)
	test_define_geometry(a, Nx, Nt, T, N_cycles)
	test_define_solution(a, q0, theta)
	test_solve()
	test_update_solution()
	test_update_pressure()
	test_compute_pressure()
	test_compute_outlet_pressure()
	test_CFL_term()
	test_check_CFL()
	test_adjust_dex()

if __name__ == '__main__':
	test_artery(sys.argv[1])
	
