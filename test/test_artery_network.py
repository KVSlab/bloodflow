import sys
sys.path.insert(0, 'arteryfe/')

import numpy as np
from configparser import SafeConfigParser

import pytest

import test_artery as ta
from arteryfe.artery_network import ArteryNetwork
from arteryfe.utils import *
from arteryfe.utils import is_near as near


def test_constructor(arterynetwork, param):
    """Construct artery network.
    Test correct assignment of parameters.
    Test correct structure of network.
    :param arterynetwork: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    assert(an.order == order)
    assert(len(an.arteries) == 2**order-1)
    assert(an.range_arteries == range(2**order-1))
    assert(an.range_parent_arteries == range(2**(order-1)-1))
    assert(an.range_daughter_arteries == range(1, 2**order-1))
    assert(an.range_end_arteries == range(2**(order-1)-1, 2**order-1))
    assert(near(an.rc, rc))
    assert(near(an.qc, qc))
    assert(near(an.rho, rho))
    assert(near(an.R1, R1))
    assert(near(an.R2, R2))
    assert(near(an.CT, CT))

    for i, artery in enumerate(an.arteries):
        assert(artery.root_vessel if i==0 else not artery.root_vessel)
        assert(artery.end_vessel if i in an.range_end_arteries\
               else not artery.end_vessel)

    return(an)


def test_define_geometry(arterynetwork, param):
    """Define geometry on artery network.
    Test correct assignment of parameters.
    :param arterynetwork: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    an.define_geometry(Nx, Nt, T, N_cycles)

    assert(an.Nx == Nx)
    assert(an.Nt == Nt)
    assert(near(an.T, T))
    assert(an.N_cycles == N_cycles)


def test_define_solution(arterynetwork, param):
    """Define solution on artery network.
    Test correct assignment of parameters.
    :param arterynetwork: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    an.define_geometry(Nx, Nt, T, N_cycles)
    an.define_solution(output_location, q0, theta)

    assert(an.output_location == output_location)
    assert(near(an.theta, theta))


def test_daughter_arteries(arterynetwork_def, param):
    """Test correct finding of daughter vessels.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for ip in an.range_parent_arteries:
        i1, i2 = an.daughter_arteries(ip)
        assert(i1 == 2*ip+1)
        assert(i2 == 2*ip+2)


def test_parent_artery(arterynetwork_def, param):
    """Test correct indices for parent vessels.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for i in an.range_daughter_arteries:
        ip = an.parent_artery(i)
        if i % 2 == 1:
            assert(ip == (i-1)//2)
        else:
            assert(ip == (i-2)//2)


def test_sister_artery(arterynetwork_def, param):
    """Test correct indices for sister vessel.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for i in an.range_daughter_arteries:
        s = an.sister_artery(i)
        if i % 2 == 1:
            assert(s == i+1)
        else:
            assert(s == i-1)


def test_flux(arterynetwork_def, param):
    """Test correct behaviour of flux function.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for a in an.arteries:
        for x in np.linspace(0, a.L, 100):
            U = a.Un(x)
            anflux = an.flux(a, U, x)
            F1 = U[1]
            F2 = U[1]**2 + a.f(x)*np.sqrt(a.A0(x)*U[0])
            assert(near(anflux[0], F1))
            assert(near(anflux[1], F2))


def test_source(arterynetwork_def, param):
    """Test correct behaviour of source function.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for a in an.arteries:
        for x in np.linspace(0, a.L, 100):
            U = a.Un(x)
            ansource = an.source(a, U, x)
            S2 = - 2*np.sqrt(np.pi)/a.db/a.Re*U[1]/np.sqrt(U[0]) + (2*np.sqrt(U[0])*(np.sqrt(np.pi)*a.f(x) + np.sqrt(a.A0(x))*a.dfdr(x)) - U[0]*a.dfdr(x))*a.drdx(x)
            assert(near(ansource[0], 0))
            assert(near(ansource[1], S2))


def test_compute_U_half(arterynetwork_def, param):
    """Test correct behaviour of compute_U_half.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for a in an.arteries:
        for x in [[0, a.dex], [a.L-a.dex, a.L]]:
            U0, U1 = a.Un(x[0]), a.Un(x[1])
            an_U_half = an.compute_U_half(a, x[0], x[1], U0, U1)
            F0, S0 = an.flux(a, U0, x[0]), an.source(a, U0, x[0])
            F1, S1 = an.flux(a, U1, x[1]), an.source(a, U1, x[1])
            U_half = (U0+U1)/2 - a.dt/a.dex*(F1-F0) + a.dt/4*(S0+S1)
            assert(near(an_U_half[0], U_half[0]))
            assert(near(an_U_half[1], U_half[1]))


def test_compute_A_out(arterynetwork_def, param):
    """Test correct behaviour of compute_A_out.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for a in an.arteries:
        pn = a.compute_outlet_pressure(a.Un(a.L)[0])
        A_out = an.compute_A_out(a)
        p = a.compute_outlet_pressure(A_out)

        # Relative pressure variation for one iteration should be reasonable
        assert(near(p, pn, reltol=1.e-1))


def test_initial_x(arterynetwork_def, param):
    """Test correct assignment of parameters.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for ip in an.range_parent_arteries:
        i1, i2 = an.daughter_arteries(ip)
        p, d1, d2 = an.arteries[ip], an.arteries[i1], an.arteries[i2]
        x = an.initial_x(p, d1, d2)
        for xi in x[:3]: assert(near(xi, p.q0))
        for xi in x[3:6]: assert(near(xi, d1.q0))
        for xi in x[6:9]: assert(near(xi, d2.q0))
        for xi in x[9:12]: assert(near(xi, p.A0(p.L)))
        for xi in x[12:15]: assert(near(xi, d1.A0(0)))
        for xi in x[15:]: assert(near(xi, d2.A0(0)))


def test_define_x(arterynetwork_def, param):
    """Test correct value for x.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    an.define_x()
    for ip in an.range_parent_arteries:
        i1, i2 = an.daughter_arteries(ip)
        p, d1, d2 = an.arteries[ip], an.arteries[i1], an.arteries[i2]
        x = an.initial_x(p, d1, d2)
        for i in range(18):
            assert(near(an.x[ip, i], x[i]))


def test_problem_function(arterynetwork_def, param):
    """Test correct behaviour of problem_function.
    For the right (analytical) value of x, the function should take zero-values.
    By perturbing a given x, certain components should be zero.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    x = np.ones(18)
    for ip in an.range_parent_arteries:
        i1, i2 = an.daughter_arteries(ip)
        p, d1, d2 = an.arteries[ip], an.arteries[i1], an.arteries[i2]

        A0p, A01, A02 = p.A0(p.L), d1.A0(0), d2.A0(0)
        fp, f1, f2 = p.f(p.L), d1.f(0), d2.f(0)

        Um1p, Um0p = p.Un(p.L-p.dex), p.Un(p.L)
        U0d1, U1d1 = d1.Un(0), d1.Un(d1.dex)
        U0d2, U1d2 = d2.Un(0), d2.Un(d2.dex)

        U_half_p = an.compute_U_half(p, p.L-p.dex, p.L, Um1p, Um0p)
        U_half_d1 = an.compute_U_half(d1, 0, d1.dex, U0d1, U1d1)
        U_half_d2 = an.compute_U_half(d2, 0, d2.dex, U0d2, U1d2)

        F_half_p = an.flux(p, U_half_p, p.L-p.dex/2)
        F_half_d1 = an.flux(d1, U_half_d1, d1.dex/2)
        F_half_d2 = an.flux(d2, U_half_d2, d2.dex/2)
        S_half_p = an.source(p, U_half_p, p.L-p.dex/2)
        S_half_d1 = an.source(d1, U_half_d1, d1.dex/2)
        S_half_d2 = an.source(d2, U_half_d2, d2.dex/2)

        # Abbreviation
        def F(x):
            return an.problem_function(p, d1, d2, x)

        # 0
        x[1] = (U_half_p[1] + x[2])/2
        assert(near(F(x)[0], 0))
        x[1] = 1

        # 1
        x[4] = (x[5] + U_half_d1[1])/2
        assert(near(F(x)[1], 0))
        x[4] = 1

        # 2
        x[7] = (x[8] + U_half_d2[1])/2
        assert(near(F(x)[2], 0))
        x[7] = 1

        # 3
        x[10] = (U_half_p[0] + x[11])/2
        assert(near(F(x)[3], 0))
        x[10] = 1

        # 4
        x[13] = (x[14] + U_half_d1[0])/2
        assert(near(F(x)[4], 0))
        x[13] = 1

        # 5
        x[16] = (x[17] + U_half_d2[0])/2
        assert(near(F(x)[5], 0))
        x[16] = 1

        # 6
        x[0] = x[3] + x[6]
        assert(near(F(x)[6], 0))
        x[0] = 1

        # 7
        x[1] = x[4] + x[7]
        assert(near(F(x)[7], 0))
        x[1] = 1

        # 8
        x[13] = A01/(1 - fp/f1*(1-np.sqrt(A0p/x[10])))**2
        assert(near(F(x)[8], 0))
        x[13] = 1

        # 9
        x[16] = A02/(1 - fp/f2*(1-np.sqrt(A0p/x[10])))**2
        assert(near(F(x)[9], 0))
        x[16] = 1

        # 10
        x[12] = A01/(1 - fp/f1*(1-np.sqrt(A0p/x[9])))**2
        assert(near(F(x)[10], 0))
        x[12] = 1

        # 11
        x[15] = A02/(1 - fp/f2*(1-np.sqrt(A0p/x[9])))**2
        assert(near(F(x)[11], 0))
        x[15] = 1

        # Ghost half terms
        Fp = an.flux(p, np.array([x[11], x[2]]), p.L + p.dex/2)
        F1 = an.flux(d1, np.array([x[14], x[5]]), -d1.dex/2)
        F2 = an.flux(d2, np.array([x[17], x[8]]), -d2.dex/2)
        Sp = an.source(p, np.array([x[11], x[2]]), p.L + p.dex/2)
        S1 = an.source(d1, np.array([x[14], x[5]]), -d1.dex/2)
        S2 = an.source(d2, np.array([x[17], x[8]]), -d2.dex/2)

        # 12
        x[0] = Um0p[1] - p.dt/p.dex*(Fp[1]-F_half_p[1])\
             + p.dt/2*(Sp[1]+S_half_p[1])
        assert(near(F(x)[12], 0))
        x[0] = 1

        # 13
        x[3] = U0d1[1] - d1.dt/d1.dex*(F_half_d1[1]-F1[1])\
             + d1.dt/2*(S_half_d1[1]+S1[1])
        assert(near(F(x)[13], 0))
        x[3] = 1

        # 14
        x[6] = U0d2[1] - d2.dt/d2.dex*(F_half_d2[1]-F2[1])\
             + d2.dt/2*(S_half_d2[1]+S2[1])
        assert(near(F(x)[14], 0))
        x[6] = 1

        # 15
        x[9] = Um0p[0] - p.dt/p.dex*(Fp[0]-F_half_p[0])\
             + p.dt/2*(Sp[0]+S_half_p[0])
        assert(near(F(x)[15], 0))
        x[9] = 1

        # 16
        x[12] = U0d1[0] - d1.dt/d1.dex*(F_half_d1[0]-F1[0])\
              + d1.dt/2*(S_half_d1[0]+S1[0])
        assert(near(F(x)[16], 0))
        x[12] = 1

        # 17
        x[15] = U0d2[0] - d2.dt/d2.dex*(F_half_d2[0]-F2[0])\
              + d2.dt/2*(S_half_d2[0]+S2[0])
        assert(near(F(x)[17], 0))
        x[15] = 1


def test_jacobian(arterynetwork_def, param):
    """Test that the analytical expression for the jacobian matrix is close to a numerically computed jacobian.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    # Higher tolerance since the numerical jacobian is inaccurate
    tol = 1.e-3
    reltol = 1.e-3

    # h gets absorbed in x if it is too small
    h = 1.e-8

    for ip in an.range_parent_arteries:
        i1, i2 = an.daughter_arteries(ip)
        p, d1, d2 = an.arteries[ip], an.arteries[i1], an.arteries[i2]
        x = np.ones(18)
        he = np.zeros(18)
        J = an.jacobian(p, d1, d2, x)
        F = an.problem_function(p, d1, d2, x)
        for j in range(18):
            he[j] = h
            Fph = an.problem_function(p, d1, d2, x+he)
            dF = (Fph-F)/h
            for i in range(18):
                assert(near(J[i, j], dF[i], tol, reltol))
            he[j] = 0


def test_newton(arterynetwork_def, param):
    """Test correct results from newton function.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for ip in an.range_parent_arteries:
        i1, i2 = an.daughter_arteries(ip)
        p, d1, d2 = an.arteries[ip], an.arteries[i1], an.arteries[i2]
        x = an.initial_x(p, d1, d2)
        x = an.newton(p, d1, d2, x, 1000, 1.e-14)
        FF = an.problem_function(p, d1, d2, x)

        # After solving, all components of FF should be zero
        for F in FF: assert(near(F, 0, 1.e-11))

        # x represents areas and flows, that should be strictly positive
        for xi in x: assert(xi > 1.e-12)


def test_adjust_bifurcation_step(arterynetwork_def, param):
    """Test correct behaviour of adjust_bifurcation_step function.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for ip in an.range_parent_arteries:
        i1, i2 = an.daughter_arteries(ip)
        p, d1, d2 = an.arteries[ip], an.arteries[i1], an.arteries[i2]
        for margin in [0.1, 0.05, 1.e-4, 1.e-8]:
            an.adjust_bifurcation_step(p, d1, d2, margin)
            assert(p.check_CFL(p.L, p.Un(p.L)[0], p.Un(p.L)[1]))
            assert(d1.check_CFL(0, d1.Un(0)[0], d1.Un(0)[1]))
            assert(d2.check_CFL(0, d2.Un(0)[0], d2.Un(0)[1]))


def test_set_inner_bc(arterynetwork_def, param):
    """Test correct assignment of inner boundary conditions.
    Test that the CFL condition is verified.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    for ip in an.range_parent_arteries:
        i1, i2 = an.daughter_arteries(ip)
        p, d1, d2 = an.arteries[ip], an.arteries[i1], an.arteries[i2]

        an.define_x()
        an.set_inner_bc(ip, i1, i2)

        # Test that the CFL-condition is verified
        assert(p.check_CFL(p.L, p.Un(p.L)[0], p.Un(p.L)[1]))
        assert(d1.check_CFL(0, d1.Un(0)[0], d1.Un(0)[1]))
        assert(d2.check_CFL(0, d2.Un(0)[0], d2.Un(0)[1]))

        x = an.initial_x(p, d1, d2)
        x = an.newton(p, d1, d2, x, 1000, 1.e-14)

        assert(near(p.U_out[0], x[9]))
        assert(near(p.U_out[1], x[0]))
        assert(near(d1.U_in[0], x[12]))
        assert(near(d1.U_in[1], x[3]))
        assert(near(d2.U_in[0], x[15]))
        assert(near(d2.U_in[1], x[6]))


def test_set_bcs(arterynetwork_def, param):
    """Test correct assignment of boundary conditions.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    q_in = an.arteries[0].q0
    an.set_bcs(q_in)

    assert(near(an.arteries[0].q_in, q_in))


def test_dump_metadata(arterynetwork_def, param):
    """Test correct execution of dump_metadata.
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    an.dump_metadata(Nt_store, N_cycles_store, store_area, store_pressure)

    order, Nx, Nt, T0, T, L, rc, qc, rho, mesh_locations, names, locations = \
        read_output(an.output_location+'/data.cfg')

    assert(order == an.order)
    assert(Nx == an.Nx)
    assert(Nt == Nt_store*N_cycles_store)
    assert(near(T0, an.T*(an.N_cycles-N_cycles_store)))
    assert(near(T, an.T*an.N_cycles))
    for i in range(len(L)): assert(near(L[i], an.arteries[i].L))
    assert(near(rc, an.rc))
    assert(near(qc, an.qc))
    assert(near(rho, an.rho))
    for i in range(len(mesh_locations)):
        assert(mesh_locations[i] ==\
               ('%s/mesh_%i.xml.gz' % (an.output_location, i)))

    i = 0
    assert(names[i] == 'flow')
    if store_area:
        i += 1
        assert(names[i] == 'area')
    if store_pressure:
        i += 1
        assert(names[i] == 'pressure')

    for i in range(len(locations)):
        assert(locations[i] == ('%s/%s' % (an.output_location, names[i])))


def test_solve(arterynetwork_def, param):
    """
    :param arterynetwork_def: Artery network object
    :param param: Config parameters
    """
    an = arterynetwork_def
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param

    q_first = np.linspace(q0, q_half, an.Nt//2)
    q_second = np.linspace(q_half, q0, an.Nt//2)

    q_ins = np.concatenate([q_first, q_second])
    an.solve(q_ins, Nt_store, N_cycles_store, store_area, store_pressure)

    for artery in arterynetwork_def.arteries:

        if artery.root_vessel:
            assert(near(artery.U(0)[1], artery.q_in))
        else:
            assert(near(artery.U(0)[0], artery.U_in[0]))
            assert(near(artery.U(0)[1], artery.U_in[1]))

        if artery.end_vessel:
            assert(near(artery.U(artery.L)[0], artery.A_out))
        else:
            assert(near(artery.U(artery.L)[0], artery.U_out[0]))
            assert(near(artery.U(artery.L)[1], artery.U_out[1]))


@pytest.fixture
def param(config_location):
    config = SafeConfigParser()
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
    q_half = config.getfloat('Solution', 'q_half')

    # Nondimensionalise parameters
    Ru, Rd, L, k1, k2, k3, Re, nu, p0, R1, R2, CT, q0, T =\
        nondimensionalise_parameters(rc, qc, Ru, Rd, L, k1, k2, k3,
                                     rho, nu, p0, R1, R2, CT, q0, T)
    q_half = nondimensionalise(rc, qc, rho, q_half, 'flow')

    param = order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
            Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
            N_cycles_store, store_area, store_pressure, q0, q_half

    return param


@pytest.fixture
def config_location():
    fconfig = '/tmp/config.cfg'
    f = open(fconfig, 'w')
    f.write(FCONF)
    f.close()
    return fconfig


@pytest.fixture
def arterynetwork(param):
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param
    an = ArteryNetwork(order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0,
                        R1, R2, CT)
    return an


@pytest.fixture
def arterynetwork_def(arterynetwork, param):
    order, rc, qc, Ru, Rd, L, k1, k2, k3, rho, Re, nu, p0, R1, R2, CT,\
        Nt, Nx, T, N_cycles, output_location, theta, Nt_store,\
        N_cycles_store, store_area, store_pressure, q0, q_half = param
    arterynetwork_def = arterynetwork
    arterynetwork_def.define_geometry(Nx, Nt, T, N_cycles)
    arterynetwork_def.define_solution(output_location, q0, theta)
    return arterynetwork_def


FCONF = """
[Parameters]
order = 2
rc = 1
qc = 10
Ru = 0.37,0.177,0.177
Rd = 0.37,0.17,0.17
L = 20.8,17.7,17.6
k1 = 2.0e7
k2 = -22.53
k3 = 8.65e5
rho = 1.06
nu = 0.046
p0 = 119990.131579
R1 = 25300
R2 = 13900
CT = 1.3384e-6

[Geometry]
Nx = 400
Nt = 4000
T = 1
N_cycles = 1

[Solution]
inlet_flow_location = data/example_inlet.csv
output_location = /tmp/
theta = 0.51
Nt_store = 200
N_cycles_store = 1
store_area = 1
store_pressure = 1
q0 = 2.0
q_half = 5.0
"""
