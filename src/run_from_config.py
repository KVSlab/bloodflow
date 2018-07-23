import sys
import numpy as np
from scipy.interpolate import interp1d

from fenics import *
import configparser

sys.path.insert(0, 'src/')

from utils import *
from artery_network import Artery_Network


def main(config_location):
	"""
	:param string config_location: Location of config file
	"""
	config = configparser.ConfigParser()
	config.read(config_location)

	# Parameters
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
	nu = config.getfloat('Parameters', 'nu')
	p0 = config.getfloat('Parameters', 'p0')
	R1 = config.getfloat('Parameters', 'R1')
	R2 = config.getfloat('Parameters', 'R2')
	CT = config.getfloat('Parameters', 'CT')

	# Geometry parameters
	Nt = config.getint('Geometry', 'Nt')
	Nx = config.getint('Geometry', 'Nx')
	N_cycles = config.getint('Geometry', 'N_cycles')

	# Solution parameters
	inlet_flow_location = config.get('Solution', 'inlet_flow_location')
	output_location = config.get('Solution', 'output_location')
	theta = config.getfloat('Solution', 'theta')
	Nt_store = config.getint('Solution', 'Nt_store')
	N_cycles_store = config.getint('Solution', 'N_cycles_store')
	store_area = config.getint('Solution', 'store_area')
	store_pressure = config.getint('Solution', 'store_pressure')

	# Import inlet flow data
	data_q = np.genfromtxt(inlet_flow_location, delimiter = ',')
	tt = data_q[:, 0]
	qq = data_q[:, 1]
	T = data_q[-1, 0]

	# Interpolate inlet flow
	q = interp1d(tt, qq)
	t = np.linspace(0, T, Nt)
	q_ins = q(t)
	#q_ins = np.zeros(Nt)

	# Adimensionalise
	Ru, Rd, L, k1, k2, k3, Re, nu, p0, R1, R2, CT, q_ins, T  =\
		adimensionalise(rc, qc, Ru, Rd, L, k1, k2, k3,
						rho, nu, p0, R1, R2, CT, q_ins, T)

	# Create artery network
	an = Artery_Network(order, rc, qc, Ru, Rd, L, k1, k2,
						k3,	rho, Re, nu, p0, R1, R2, CT)
	an.define_geometry(Nx, Nt, T, N_cycles)
	an.define_solution(output_location, q_ins[0], theta)
	an.solve(q_ins, Nt_store, N_cycles_store, store_area, store_pressure)

if __name__ == '__main__':
	main(sys.argv[1])
