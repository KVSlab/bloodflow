import sys

sys.path.insert(0, 'src/')
from artery_network import Artery_Network

def get_parameters(config_location):

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
