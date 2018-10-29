import sys
import numpy as np
import numpy.linalg as npl

import configparser
from dolfin import *

from arteryfe.artery import Artery
from arteryfe.utils import *


class ArteryNetwork(object):
    """
    Builds an artery network from the given parameters. Arteries in the network
    are assigned indices from left to right and top to bottomself.

    Arguments
    ---------
    order : int
        Number of arterial levels
    rc : float
        Characteristic radius (length)
    qc : float
        Characteristic flow
    Ru : list of float
        Upstream radii of each artery in the network
    Rd : list of float
        Downstream radii of each artery in the network
    L : list of float
        Vessel lengths
    k1 : float
     	First constant from the relation Eh/r0
    k2 : float
        Second constant from the relation Eh/r0
    k3 : float
        Third constant from the relation Eh/R0
    rho : float
        Density of blood
    Re : float
        Reynolds' number
    nu : float
        Viscosity of blood
    p0 : float
        Diastolic pressure
    R1 : float
        First resistance from Windkessel model
    R2 : float
        Second resistance from Windkessel model
    CT : float
        Compliance from Windkessel model
    """


    def __init__(self, order, rc, qc, Ru, Rd, L, k1, k2, k3,
                                                rho, Re, nu, p0, R1, R2, CT):
        set_log_level(30)
        self.order = order
        self.arteries = [0] * (2**self.order-1)
        self.range_arteries = range(2**self.order-1)
        self.range_parent_arteries = range(2**(self.order-1)-1)
        self.range_daughter_arteries = range(1, 2**self.order-1)
        self.range_end_arteries = range(2**(self.order-1)-1, 2**self.order-1)
        self.rc, self.qc, self.rho = rc, qc, rho
        self.R1, self.R2, self.CT = R1, R2, CT
        for i in self.range_arteries:
            root_vessel = (i==0)
            end_vessel = (i in self.range_end_arteries)
            self.arteries[i] = Artery(root_vessel, end_vessel, rc, qc, Ru[i],
                                      Rd[i], L[i], k1, k2, k3, rho, Re, nu, p0)


    def daughter_arteries(self, i):
        """
        Find and return the indices of the daughter arteries of artery i.

        Arguments
        ---------
        i : int
            Index of the parent artery

        Returns
        -------
        return : int
            Daughter artery indices
        """
        return 2*i+1, 2*i+2


    def parent_artery(self, i):
        """
        Find and return the index of the partent artery of artery i.

        Arguments
        ---------
        i : int
            Index of the daughter artery

        Returns
        -------
        return : int
            Parent artery index
        """
        #if i <= 0 or i >= 2**self.order:
        #	raise Exception('Vessel index out of range')
        return (i-1)//2  # d1 is odd, d2=d1+1 is even


    def sister_artery(self, i):
        """
        Find and return the index of the sister artery of artery i.

        Arguments
        ---------
        i : int
            Index of the artery

        Returns
        -------
        return : int
            Sister artery index
        """
        if i%2 == 0:
            return i-1
        else:
            return i+1


    def define_geometry(self, Nx, Nt, T, N_cycles):
        """
        Calls define_geometry() for each artery in the network.

        Arguments
        ---------
        Nx : int
            Number of spatial points per artery
        Nt : int
            Number of time steps per cardiac cycle
        T : float
            Duration of one cardiac cycle
        N_cycles: int
            Number of cardiac cycles in the simulation
        """
        self.Nx = Nx
        self.Nt = Nt
        self.T = T
        self.N_cycles = N_cycles
        self.dt = T/Nt
        for i in self.range_arteries:
            self.arteries[i].define_geometry(Nx, Nt, T, N_cycles)


    def define_solution(self, output_location, q0, theta=0.5):
        """
        Calls define_solution() for each artery in the network.

        Arguments
        ---------
        output_location : string
            Output data directory
        q0 : float
            Initial flow rate in the root vessel
        theta : float
            Weighting parameter for the Crank-Nicolson method, in the interval
            [0, 1]
        """
        self.output_location = output_location
        self.theta = theta
        self.arteries[0].define_solution(q0, theta)
        for i in self.range_daughter_arteries:
            p = self.parent_artery(i)
            s = self.sister_artery(i)
            q0 = self.arteries[i].A0(0)/(self.arteries[i].A0(0)\
                                        +self.arteries[s].A0(0))\
               * self.arteries[p].q0
            self.arteries[i].define_solution(q0, theta)

        self.x = np.ones([len(self.range_parent_arteries), 18])


    def flux(self, a, U, x):
        """
        Computes the flux term F(U) for a given artery a in the network.

        Arguments
        ---------
        a : Artery
            Artery on which the flux term is computed
        U : numpy.array
            Value of solution
        x : float
            Point of evaluation

        Returns
        -------
        return : numpy.array
            Flux term F(U) for artery a at point x
        """
        return np.array([U[1], U[1]**2 + a.f(x)*np.sqrt(a.A0(x)*U[0])])


    def source(self, a, U, x):
        """
        Computes the source term S(U) for a given artery a in the network.

        Arguments
        ---------
        a : Artery
            Artery on which the flux term is computed
        U : numpy.array
            Value of solution
        x : float
            Point of evaluation

        Returns
        -------
        return : numpy.array
            Source term S(U) for artery a at point x
        """
        S1 = 0
        S2 = -2*np.sqrt(np.pi)/a.db/a.Re*U[1]/np.sqrt(U[0])\
           + (2*np.sqrt(U[0])*(np.sqrt(np.pi)*a.f(x)\
                              +np.sqrt(a.A0(x))*a.dfdr(x))\
             -U[0]*a.dfdr(x))*a.drdx(x)
        return np.array([S1, S2])


    def compute_U_half(self, a, x0, x1, U0, U1):
        """
        Computes the solution for a given artery a in the network at a half
        time step.

        Arguments
        ---------
        a : Artery
            Artery on which the flux term is computed
        x0 : float
            Left point
        x1 : float
            Right point
        U0 : numpy.array
            Value of solution at x0
        U1 : numpy.array
            Value of solution at x1
        x : float
            Point of evaluation

        Returns
        -------
        return : numpy.array
            Solution for artery a at the middle point after half a time step
        """
        # Value of terms at time t_n
        F0, S0 = self.flux(a, U0, x0), self.source(a, U0, x0)
        F1, S1 = self.flux(a, U1, x1), self.source(a, U1, x1)

        return (U0+U1)/2 - a.dt/(x1-x0)*(F1-F0) + a.dt/4*(S0+S1)


    def compute_A_out(self, a, k_max=100, tol=1.0e-12):
        """
        Computes the area for artery a at the outlet

        Arguments
        ---------
        a : Artery
            Artery on which the flux term is computed
        k_max : int
            Maximum number of iterations in Piccards scheme
        tol : float
            Tolerance for Piccards fixed point iteration scheme

        Returns
        -------
        return : float
            Outlet boundary value of A at time step t_(n+1)
        """
        a.adjust_dex(a.L, a.Un(a.L)[0], a.Un(a.L)[1])

        # Spatial step, scaled to satisfy the CFL condition
        x2, x1, x0 = a.L-2*a.dex, a.L-a.dex, a.L
        x21, x10 = a.L-1.5*a.dex, a.L-0.5*a.dex
        Um2, Um1, Um0 = a.Un(x2), a.Un(x1), a.Un(x0)

        # Values at time step n
        Fm2, Sm2 = self.flux(a, Um2, x2), self.source(a, Um2, x2)
        Fm1, Sm1 = self.flux(a, Um1, x1), self.source(a, Um1, x1)
        Fm0, Sm0 = self.flux(a, Um0, x0), self.source(a, Um0, x0)

        # Values at time step n+1/2
        U_half_21 = self.compute_U_half(a, x2, x1, Um2, Um1)
        U_half_10 = self.compute_U_half(a, x1, x0, Um1, Um0)
        F_half_21 = self.flux(a, U_half_21, x21)
        S_half_21 = self.source(a, U_half_21, x21)
        F_half_10 = self.flux(a, U_half_10, x10)
        S_half_10 = self.source(a, U_half_10, x10)

        # Value at time step n+1
        qm1 = Um1[1]\
            - a.dt/a.dex*(F_half_10[1]-F_half_21[1])\
            + a.dt/2*(S_half_10[1]+S_half_21[1])

        # Fixed point iteration
        pn = a.compute_outlet_pressure(Um0[0])
        p = pn
        for k in range(k_max):
            p_old = p
            qm0 = Um0[1]\
                + (p-pn)/self.R1\
                + self.dt/self.R1/self.R2/self.CT*pn\
                - self.dt*(self.R1+self.R2)/self.R1/self.R2/self.CT*Um0[1]
            Am0 = Um0[0] - self.dt/a.dex*(qm0-qm1)
            p = a.compute_outlet_pressure(Am0)
            if abs(p-p_old) < tol:
                break

        return Am0


    def initial_x(self, p, d1, d2):
        """
        Computes an initial guess for x at a bifurcation point. Set same value
        at time t_(n+1) and t_(n+1/2) as time t_n. At point M+-1/2, set same
        value as in point M.

        Arguments
        ---------
        p : Artery
            Parent artery in the bifurcation
        d1 : Artery
            First daughter artery in the bifurcation
        d2 : Artery
            Second daughter artery in the bifurcation

        Returns
        -------
        return : numpy.array
            Initial guess for the 18 variables at a bifurcation
        """
        x = np.zeros(18)
        x[:3] = p.q0*np.ones(3)
        x[3:6] = d1.q0*np.ones(3)
        x[6:9] = d2.q0*np.ones(3)
        x[9:12] = p.A0(p.L)*np.ones(3)
        x[12:15] = d1.A0(0)*np.ones(3)
        x[15:] = d2.A0(0)*np.ones(3)
        return x


    def define_x(self):
        """
        Calls initial_x() for each bifurcation in the artery network
        """
        for ip in self.range_parent_arteries:
            i1, i2 = self.daughter_arteries(ip)
            p, d1, d2 = self.arteries[ip], self.arteries[i1], self.arteries[i2]
            self.x[ip] = self.initial_x(p, d1, d2)


    def problem_function(self, p, d1, d2, x):
        """
        Computes the solution f(x) = y of the system of equations at a
        bifurcation.

        Arguments
        ---------
        p : Artery
            Parent artery in the bifurcation
        d1 : Artery
            First daughter artery in the bifurcation
        d2 : Artery
            Second daughter artery in the bifurcation
        x : numpy.array
            Current solution

        Returns
        -------
        return : numpy.array
            Function value f(x)
        """
        # Abbreviations
        A0p, A01, A02 = p.A0(p.L), d1.A0(0), d2.A0(0)
        fp, f1, f2 = p.f(p.L),  d1.f(0), d2.f(0)

        # Ghost half terms
        Fp = self.flux(p, np.array([x[11], x[2]]), p.L+p.dex/2)
        F1 = self.flux(d1, np.array([x[14], x[5]]), -d1.dex/2)
        F2 = self.flux(d2, np.array([x[17], x[8]]), -d2.dex/2)
        Sp = self.source(p, np.array([x[11], x[2]]), p.L+p.dex/2)
        S1 = self.source(d1, np.array([x[14], x[5]]), -d1.dex/2)
        S2 = self.source(d2, np.array([x[17], x[8]]), -d2.dex/2)

        # Compute half-time-step-values in M-1/2 for p and 1/2 for d1 and d2
        Um1p, Um0p = p.Un(p.L-p.dex), p.Un(p.L)
        U0d1, U1d1 = d1.Un(0), d1.Un(d1.dex)
        U0d2, U1d2 = d2.Un(0), d2.Un(d2.dex)

        U_half_p = self.compute_U_half(p, p.L-p.dex, p.L, Um1p, Um0p)
        U_half_1 = self.compute_U_half(d1, 0, d1.dex, U0d1, U1d1)
        U_half_2 = self.compute_U_half(d2, 0, d2.dex, U0d2, U1d2)

        F_half_p = self.flux(p, U_half_p, p.L-p.dex/2)
        S_half_p = self.source(p, U_half_p, p.L-p.dex/2)
        F_half_1 = self.flux(d1, U_half_1, d1.dex/2)
        S_half_1 = self.source(d1, U_half_1, d1.dex/2)
        F_half_2 = self.flux(d2, U_half_2, d2.dex/2)
        S_half_2 = self.source(d2, U_half_2, d2.dex/2)

        # Function return array
        y = np.zeros(18)

        # Entries from equation (20)
        y[0] = 2*x[1] - U_half_p[1] - x[2]
        y[1] = 2*x[4] - x[5] - U_half_1[1]
        y[2] = 2*x[7] - x[8] - U_half_2[1]

        # Entries from equation (21)
        y[3] = 2*x[10] - U_half_p[0] - x[11]
        y[4] = 2*x[13] - x[14] - U_half_1[0]
        y[5] = 2*x[16] - x[17] - U_half_2[0]

        # Entries from equation (22)
        y[6] = x[0] - x[3] - x[6]
        y[7] = x[1] - x[4] - x[7]

        # Entries from equation (23)
        y[8] = fp*(1-np.sqrt(A0p/x[10])) - f1*(1-np.sqrt(A01/x[13]))
        y[9] = fp*(1-np.sqrt(A0p/x[10])) - f2*(1-np.sqrt(A02/x[16]))
        y[10] = fp*(1-np.sqrt(A0p/x[9])) - f1*(1-np.sqrt(A01/x[12]))
        y[11] = fp*(1-np.sqrt(A0p/x[9])) - f2*(1-np.sqrt(A02/x[15]))

        # Entries from equation (26)
        y[12] = x[0] - Um0p[1] + p.dt/p.dex*(Fp[1] - F_half_p[1])\
              - p.dt/2*(Sp[1] + S_half_p[1])
        y[13] = x[3] - U0d1[1] + d1.dt/d1.dex*(F_half_1[1] - F1[1])\
              - d1.dt/2*(S_half_1[1] + S1[1])
        y[14] = x[6] - U0d2[1] + d2.dt/d2.dex*(F_half_2[1] - F2[1])\
              - d2.dt/2*(S_half_2[1] + S2[1])

        # Entries from equation (27)
        y[15] = x[9] - Um0p[0] + p.dt/p.dex*(Fp[0] - F_half_p[0])
        y[16] = x[12] - U0d1[0] + d1.dt/d1.dex*(F_half_1[0] - F1[0])
        y[17] = x[15] - U0d2[0] + d2.dt/d2.dex*(F_half_2[0] - F2[0])

        return y


    def jacobian(self, p, d1, d2, x):
        """
        Computes the analytical Jacobian of the system of equations at a
        bifurcation.

        Arguments
        ---------
        p : Artery
            Parent artery in the bifurcation
        d1 : Artery
            First daughter artery in the bifurcation
        d2 : Artery
            Second daughter artery in the bifurcation
        x : numpy.array
            Current solution

        Returns
        -------
        return : numpy.array
            Jacobian matrix
        """
        # Abbreviations
        A0p, A01, A02 = p.A0(p.L), d1.A0(0), d2.A0(0)
        A0hp, A0h1, A0h2 = p.A0(p.L+p.dex), d1.A0(-d1.dex), d2.A0(-d2.dex)
        fp, f1, f2 = p.f(p.L),  d1.f(0), d2.f(0)
        fhp, fh1, fh2 = p.f(p.L+p.dex), d1.f(-d1.dex), d2.f(-d2.dex)
        dbp, db1, db2 = p.db, d1.db, d2.db
        Rep, Re1, Re2 = p.Re, d1.Re, d2.Re
        dfdrp, dfdr1, dfdr2 = p.dfdr(p.L+p.dex),\
                              d1.dfdr(-d1.dex), d2.dfdr(-d2.dex)
        drdxp, drdx1, drdx2 = p.drdx(p.L+p.dex),\
                              d1.drdx(-d1.dex), d2.drdx(-d2.dex)
        rpi = np.sqrt(np.pi)

        # Jacobian matrix
        J = np.zeros([18, 18])

        # Entries from equation (20)
        J[0, 1] = 2
        J[0, 2] = -1
        J[1, 4] = 2
        J[1, 5] = -1
        J[2, 7] = 2
        J[2, 8] = -1

        # Entries from equation (21)
        J[3, 10] = 2
        J[3, 11] = -1
        J[4, 13] = 2
        J[4, 14] = -1
        J[5, 16] = 2
        J[5, 17] = -1

        # Entries from equation (22)
        J[6, 0] = 1
        J[6, 3] = -1
        J[6, 6] = -1
        J[7, 1] = 1
        J[7, 4] = -1
        J[7, 7] = -1

        # Entries from equation (23)
        J[8, 10] = fp*np.sqrt(A0p)/2/x[10]**(3/2)
        J[8, 13] = -f1*np.sqrt(A01)/2/x[13]**(3/2)
        J[9, 10] = fp*np.sqrt(A0p)/2/x[10]**(3/2)
        J[9, 16] = -f2*np.sqrt(A02)/2/x[16]**(3/2)
        J[10, 9] = fp*np.sqrt(A0p)/2/x[9]**(3/2)
        J[10, 12] = -f1*np.sqrt(A01)/2/x[12]**(3/2)
        J[11, 9] = fp*np.sqrt(A0p)/2/x[9]**(3/2)
        J[11, 15] = -f2*np.sqrt(A02)/2/x[15]**(3/2)

        # Entries from equation (26)
        J[12, 0] = 1
        J[12, 2] = p.dt/p.dex*2*x[2]/x[11] + p.dt*rpi/dbp/Rep/np.sqrt(x[11])
        J[12, 11] = p.dt/p.dex*(-(x[2]/x[11])**2 + fhp/2*np.sqrt(A0hp/x[11]))\
                  - p.dt/2*(rpi/dbp/Rep*x[2]/x[11]**(3/2)\
                           + (1/np.sqrt(x[11])*(rpi*fhp+np.sqrt(A0hp)*dfdrp)\
                                               - dfdrp)*drdxp)
        J[13, 3] = 1
        J[13, 5] = -d1.dt/d1.dex*2*x[5]/x[14] + d1.dt*rpi/db1/Re1/np.sqrt(x[14])
        J[13, 14] = d1.dt/d1.dex*((x[5]/x[14])**2 - fh1/2*np.sqrt(A0h1/x[14]))\
                  - d1.dt/2*(rpi/db1/Re1*x[5]/x[14]**(3/2)\
                            + (1/np.sqrt(x[14])*(rpi*fh1+np.sqrt(A0h1)*dfdr1)\
                                                - dfdr1)*drdx1)
        J[14, 6] = 1
        J[14, 8] = -d2.dt/d2.dex*2*x[8]/x[17] + d2.dt*rpi/db2/Re2/np.sqrt(x[17])
        J[14, 17] = d2.dt/d2.dex*((x[8]/x[17])**2 - fh2/2*np.sqrt(A0h2/x[17]))\
                  - d2.dt/2*(rpi/db2/Re2*x[8]/x[17]**(3/2)\
                            + (1/np.sqrt(x[17])*(rpi*fh2+np.sqrt(A0h2)*dfdr2)\
                                                - dfdr2)*drdx2)

        # Entries from equation (27)
        J[15, 2] = p.dt/p.dex
        J[15, 9] = 1
        J[16, 5] = - d1.dt/d1.dex
        J[16, 12] = 1
        J[17, 8] = - d1.dt/d1.dex
        J[17, 15] = 1

        return J


    def newton(self, p, d1, d2, x=np.ones(18), k_max=30, tol=1.e-10):
        """
        Computes the solution of the system of equations at a bifurcation
        using Newton's method.

        Arguments
        ---------
        p : Artery
            Parent artery in the bifurcation
        d1 : Artery
            First daughter artery in the bifurcation
        d2 : Artery
            Second daughter artery in the bifurcation
        x : numpy.array
            Current solution
        k_max: int
            Maximum number of iterations
        tol : float
            Tolerance for the solution

        Returns
        -------
        return : numpy.array
            Solution to the system of equations
        """
        for k in range(k_max):

            J = self.jacobian(p, d1, d2, x)
            func = self.problem_function(p, d1, d2, x)

            if npl.norm(func) < tol:
                break

            try:
                x -= npl.solve(J, func)
            except npl.LinAlgError:
                print('Singular')
                eps = 1.e-6  # Perturbation value
                J += eps*np.eye(18)
                func[0] += eps
                x -= npl.solve(J, func)

        return x


    def adjust_bifurcation_step(self, p, d1, d2, margin=0.05):
        """
        Chooses spatial step at a bifurcation to respect the CFL condition
        for all three arteries

        Arguments
        ---------
        p : Artery
            Parent artery in the bifurcation
        d1 : Artery
            First daughter artery in the bifurcation
        d2 : Artery
            Second daughter artery in the bifurcation
        margin : float
            Margin of CFL number
        """
        Mp = p.CFL_term(p.L, p.Un(p.L)[0], p.Un(p.L)[1])
        M1 = d1.CFL_term(0, d1.Un(0)[0], d1.Un(0)[1])
        M2 = d2.CFL_term(0, d2.Un(0)[0], d2.Un(0)[1])

        # dex is chosen to respect all three CFL-conditions
        p.dex = d1.dex = d2.dex = (1+margin)*self.dt/min([Mp, M1, M2])


    def set_inner_bc(self, ip, i1, i2):
        """
        Calls newton() for each bifurcation to calculate the boundary values
        for each artery at the bifurcation.

        Arguments
        ---------
        ip : int
            Parent artery index in the bifurcation
        i1 : int
            First daughter artery index in the bifurcation
        i2 : int
            Second daughter artery index in the bifurcation
        """
        p, d1, d2 = self.arteries[ip], self.arteries[i1], self.arteries[i2]

        self.adjust_bifurcation_step(p, d1, d2)

        self.x[ip] = self.newton(p, d1, d2, self.x[ip])

        p.U_out = [self.x[ip, 9], self.x[ip, 0]]
        d1.U_in = [self.x[ip, 12], self.x[ip, 3]]
        d2.U_in = [self.x[ip, 15], self.x[ip, 6]]


    def set_bcs(self, q_in):
        """
        Updates all boundary values using the appropriate boundary conditions.

        Arguments
        ---------
        q_in : float
            Inflow rate in the root vessel at time t_(n+1)
        """
        # Update inlet boundary conditions
        self.arteries[0].q_in = q_in

        # Update bifurcation boundary conditions
        for ip in self.range_parent_arteries:
            i1, i2 = self.daughter_arteries(ip)
            self.set_inner_bc(ip, i1, i2)

        # Update outlet boundary conditions
        for i in self.range_end_arteries:
            self.arteries[i].A_out = self.compute_A_out(self.arteries[i])


    def dump_metadata(self, Nt_store, N_cycles_store, store_area,
      store_pressure):
        """
        Save metadata for the interpretation of XDMF files.

        Arguments
        ---------
        Nt_store : int
            Number of time steps to store
        N_cycles_store : int
            Number of cardiac cycles to store
        store_area : boolean
            Store area if True
        store_pressure : boolean
            Store pressure if True
        """
        # Assemble strings
        mesh_locations = ''
        for i in self.range_arteries:
            mesh_location = self.output_location + '/mesh_%i.xml.gz' % (i)
            File(mesh_location) << self.arteries[i].mesh  # Save mesh
            if i > 0:
                mesh_locations += ','
            mesh_locations += mesh_location
        L = ''
        for artery in self.arteries[:-1]:
            L += str(artery.L)+','
        L += str(self.arteries[-1].L)
        names = ''
        locations = ''
        names += 'flow'
        locations += self.output_location + '/flow'
        if store_area:
            names += ',area'
            locations += ',' + self.output_location + '/area'
        if store_pressure:
            names += ',pressure'
            locations += ',' + self.output_location + '/pressure'

        # Save metadata
        config = configparser.RawConfigParser()
        config.add_section('data')
        config.set('data', 'order', str(self.order))
        config.set('data', 'Nx', str(self.Nx))
        config.set('data', 'Nt', str(Nt_store*N_cycles_store))
        config.set('data', 'T0', str(self.T*(self.N_cycles-N_cycles_store)))
        config.set('data', 'T', str(self.T*self.N_cycles))
        config.set('data', 'L', L)
        config.set('data', 'rc', str(self.rc))
        config.set('data', 'qc', str(self.qc))
        config.set('data', 'rho', str(self.rho))
        config.set('data', 'mesh_locations', mesh_locations)
        config.set('data', 'names', names)
        config.set('data', 'locations', locations)
        with open(self.output_location+'/data.cfg', 'w') as configfile:
            config.write(configfile)


    def solve(self, q_ins, Nt_store, N_cycles_store=1, store_area=False,
      store_pressure=True):
        """
        Call solve for each artery in the network.

        Arguments
        ---------
        q_ins : numpy.array
            Array containing inflow rates for the root artery at every time step
        Nt_store : int
            Number of time steps to store
        N_cycles_store : int
            Number of cardiac cycles to store
        store_area : boolean
            Store area if True
        store_pressure : boolean
            Store pressure if True
        """
        self.define_x()

        # Store parameters necessary for postprocessing
        self.dump_metadata(Nt_store, N_cycles_store, store_area, store_pressure)

        # Setup storage files
        xdmffile_flow = [0] * len(self.range_arteries)
        if store_area:
            xdmffile_area = [0] * len(self.range_arteries)
        if store_pressure:
            xdmffile_pressure = [0] * len(self.range_arteries)
        for i in self.range_arteries:
            xdmffile_flow[i] = XDMFFile('%s/flow/flow_%i.xdmf'\
                                        % (self.output_location, i))
            if store_area:
                xdmffile_area[i] = XDMFFile('%s/area/area_%i.xdmf'\
                                            % (self.output_location, i))
            if store_pressure:
                xdmffile_pressure[i] = XDMFFile('%s/pressure/pressure_%i.xdmf'\
                                                % (self.output_location, i))

        # Initialise time
        t = 0

        # Cardiac cycle iteration
        for n_cycle in range(self.N_cycles):

            # Time-stepping for one period
            for n in range(self.Nt):

                print_progress(n_cycle, n, n_cycle*self.Nt+n)

                # Apply boundary conditions for time t_(n+1)
                self.set_bcs(q_ins[(n+1) % (self.Nt)])

                # Solve equation on each artery
                for i, artery in enumerate(self.arteries):

                    # Store solution at time t_n
                    cycle_store = (n_cycle >= self.N_cycles-N_cycles_store)

                    if cycle_store and n % (self.Nt/Nt_store) == 0:

                        # Split solution for storing, with deepcopy
                        area, flow = artery.Un.split(True)

                        write_file(xdmffile_flow[i], flow, 'flow', t)

                        if store_area:
                            write_file(xdmffile_area[i], area, 'area', t)

                        if store_pressure:
                            artery.update_pressure()
                            write_file(xdmffile_pressure[i], artery.pn, 'pressure', t)

                    # Solve problem on artery for time t_(n+1)
                    artery.solve()

                    # Update current solution on artery
                    artery.update_solution()

                t += self.dt
