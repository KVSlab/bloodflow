import sys
import numpy as np

from arteryfe.utils import *


def main(data_location):
    """
    :param string data_location: Location of data file.
    """
    order, Nx, Nt, T0, T, L, rc, qc, rho, mesh_locations, names, locations =\
        read_output(data_location)

    T0 = redimensionalise(rc, qc, rho, T0, 'time')
    T = redimensionalise(rc, qc, rho, T, 'time')

    t = np.linspace(T0, T, Nt)

    for i, name in enumerate(names):

        for j in range(2**order-1):

            if L[j] > 0.0:

                M = XDMF_to_matrix(Nx, Nt, mesh_locations[j],
                    '%s/%s_%i.xdmf' % (locations[i], name, j), name)

                M = redimensionalise(rc, qc, rho, M, name)

                # Convert pressure units
                if name == 'pressure':
                    M = unit_to_mmHg(M)

                x = np.linspace(0, L[j], Nx+1)

                plot_matrix(t, x, M, name, '%s/%s_%i.png' % (locations[i], name, j))


if __name__ == '__main__':
    main(sys.argv[1])
