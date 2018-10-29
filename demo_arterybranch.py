import sys

import numpy as np
import configparser

import arteryfe as af


def main(config_location):
    """Read config-file.
    Run the necessary functions to compute the solution.
    :param string config_location: Location of config file
    """
    param = af.ParamParser(config_location)

    # Constructor parameters
    order = param.param['order']
    rc = param.param['rc']
    qc = param.param['qc']
    Ru = param.param['Ru']
    Rd = param.param['Rd']
    L = param.param['L']
    k1 = param.param['k1']
    k2 = param.param['k2']
    k3 = param.param['k3']
    rho = param.param['rho']
    nu = param.param['nu']
    p0 = param.param['p0']
    R1 = param.param['R1']
    R2 = param.param['R2']
    CT = param.param['CT']

    # Geometry parameters
    Nt = param.geo['Nt']
    Nx = param.geo['Nx']
    N_cycles = param.geo['N_cycles']

    # Solution parameters
    inlet_flow_location = param.solution['inlet_flow_location']
    output_location = param.solution['output_location']
    theta = param.solution['theta']
    Nt_store = param.solution['Nt_store']
    N_cycles_store = param.solution['N_cycles_store']
    store_area = param.solution['store_area']
    store_pressure = param.solution['store_pressure']

    # Import inlet flow data
    T, q_ins = af.read_inlet(inlet_flow_location, Nt)

    # Nondimensionalise data and compute Reynolds' number
    Ru, Rd, L, k1, k2, k3, Re, nu, p0, R1, R2, CT, q_ins, T =\
        af.nondimensionalise_parameters(rc, qc, Ru, Rd, L, k1, k2, k3,
                                   rho, nu, p0, R1, R2, CT, q_ins, T)

    # Create artery network
    an = af.ArteryNetwork(order, rc, qc, Ru, Rd, L, k1, k2,
                        k3,	rho, Re, nu, p0, R1, R2, CT)
    an.define_geometry(Nx, Nt, T, N_cycles)
    an.define_solution(output_location, q_ins[0], theta)

    # Solve problem and store data
    an.solve(q_ins, Nt_store, N_cycles_store, store_area, store_pressure)


if __name__ == '__main__':
    main(sys.argv[1])
