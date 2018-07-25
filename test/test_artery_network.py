import sys

sys.path.insert(0, 'src/')
sys.path.insert(0, 'test/')

import test_artery as ta
from artery_network import Artery_Network

def get_parameters(config_location):
	"""Read parameters for tests from file.
	:param config_location: Location of config file
	:return: The parameters needed for testing
	"""
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
	
	return rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, Nt, Nx, T, N_cycles, output_location, theta, Nt_store, N_cycles_store, store_area, store_pressure, q0




def test_constructor():

def test_daughter_vessels():

def test_parent_vessel():

def sister_vessel(self, i):

def test_flux():

def test_source():

def test_compute_U_half():
	
def test_compute_A_out():

def test_initial_x():

def test_define_x():

def test_define_geometry():

def test_define_solution():

def test_problem_function():

def test_jacobian():

def test_newton():

def test_adjust_bifurcation_step():

def test_set_inner_bc():

def test_set_bcs():

def test_dump_metadata():

def test_solve():

def test_artery_network(config_location):
	... = get_parameters(config_location)

if __name__ == '__main__':
	test_artery_network(sys.argv[1])
