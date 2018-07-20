import sys
import numpy as np

import matplotlib.pyplot as plt

sys.path.insert(0, 'src/')

from utils import *


def main(data_location):
	"""
	:param string data_location: Location of data file.
	"""
	order, Nx, Nt, T, L, rc, qc, rho, mesh_locations, locations, names =\
		read_output(data_location)

	T = redimensionalise(rc, qc, rho, T, 'time')

	t = np.linspace(0, T, Nt)

	for i, name in enumerate(names):

		for j in range(2**order-1):
			
			M = XDMF_to_matrix(Nx, Nt, mesh_locations[j],
				'%s/%s_%i.xdmf' % (locations[i], name, j), name)
				
			M = redimensionalise(rc, qc, rho, M, name)
			
			if name == 'pressure':
				M = unit_to_mmHg(M)
			
			x = np.linspace(0, L[j], Nx+1)

			plot_matrix(t, x, M, name, '%s/%s_%i.png' % (locations[i], name, j))


if __name__ == '__main__':
	main(sys.argv[1])
