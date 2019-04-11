import sys
import numpy as np
import copy

from dolfin import *

# Compiler parameters
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)


class Artery(object):
    """
    Represents an artery whose flow rate and area are calculated using the
    1D system of blood flow equations in conservation form.

    Arguments
    ---------

    root : boolean
        True if the artery is a root vessel in the artery network (has no
        parent)
    leaf : boolean
        True if the artery is a terminal vessel in the artery network (has no
        daughters)
    rc : float
        Characteristic radius (length)
    qc : float
        Characteristic flow
    Ru : float
        Upstream radius
    Rd : float
        Downstream radius
    L : float
        Vessel length
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
    """


    def __init__(self, i, T, param, root=False, leaf=False):
        self.i = i
        self.root = root
        self.leaf = leaf
        self.T = T
        self.param = copy.deepcopy(param)
        self.param['Ru'] = param['Ru'][i]
        self.param['Rd'] = param['Rd'][i]
        self.param['L'] = param['L'][i]


    def define_geometry(self, geo):
        """
        Initialises the artery geometry by creating the spatial refinement,
        temporal refinement and FEniCS objects.

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
        Nx = geo['Nx']
        Nt = geo['Nt']
        N_cycles = geo['N_cycles']
        L = self.param['L']
        Ru = self.param['Ru']
        Rd = self.param['Rd']
        k1 = self.param['k1']
        k2 = self.param['k2']
        k3 = self.param['k3']

        self.dx = L/Nx
        self.dt = self.T/Nt

        # Step for boundary condition computations, starting at normal size
        self.dex = self.dx

        self.db = np.sqrt(self.param['nu']*self.T/2/np.pi)

        self.mesh = IntervalMesh(Nx, 0, L)
        self.elV = FiniteElement('CG', self.mesh.ufl_cell(), 1)
        self.V = FunctionSpace(self.mesh, self.elV)
        self.V2 = FunctionSpace(self.mesh, self.elV*self.elV)

        # Initial vessel-radius and deduced quantities
        self.r0 = Expression('Ru*pow(Rd/Ru, x[0]/L)',
                             degree=2, Ru=Ru, Rd=Rd, L=L)
        self.A0 = Expression('pi*pow(r0, 2)', degree=2, r0=self.r0)
        self.f = Expression('4.0/3.0*(k1*exp(k2*r0) + k3)', degree=2,
                            k1=k1, k2=k2, k3=k3, r0=self.r0)
        self.dfdr = Expression('4.0/3.0*k1*k2*exp(k2*r0)', degree=2,
                               k1=k1, k2=k2, r0=self.r0)
        self.drdx = Expression('logRdRu/L*r0', degree=2,
                               logRdRu=np.log(Rd/Ru), L=L, r0=self.r0)


    def define_solution(self, q0, theta=0.5, bc_tol=1.e-14):
        """
        Defines FEniCS Function objects, boundary conditions and variational
        form of the problem.

        Arguments
        ---------
        q0 : float
            Initial flow rate in the root vessel
        theta : float
            Weighting parameter for the Crank-Nicolson method, in the interval
            [0, 1]
        bc_tol : float
            Inlet and outlet boundary thickness (tolerance)
        """
        L = self.param['L']
        p0 = self.param['p0']
        Re = self.param['Re']

        # Initial flow value
        self.q0 = q0

        # Crank-Nicolson parameter
        self.theta = theta

        # Trial function
        self.U = Function(self.V2)
        A, q = split(self.U)

        # Test functions
        v1, v2 = TestFunctions(self.V2)

        # Current solution, initialised
        self.Un = Function(self.V2)
        self.Un.assign(Expression(('A0', 'q0'), degree=2,
                                  A0=self.A0, q0=self.q0))

        # Current pressure, initialised
        self.pn = Function(self.V)
        self.pn.assign(Expression('p0', degree=2, p0=p0))

        # Boundary conditions (spatial)
        def inlet_bdry(x, on_boundary):
            return on_boundary and near(x[0], 0, bc_tol)

        def outlet_bdry(x, on_boundary):
            return on_boundary and near(x[0], L, bc_tol)

        # Inlet boundary conditions
        if self.root:
            self._q_in = Expression('value', degree=0, value=self.q0)
            bc_inlet = DirichletBC(self.V2.sub(1), self._q_in, inlet_bdry)
        else:
            self._U_in = Expression(('A', 'q'), degree=0,
                                    A=self.A0(0), q=self.q0)
            bc_inlet = DirichletBC(self.V2, self._U_in, inlet_bdry)

        # Outlet boundary conditions
        if self.leaf:
            self._A_out = Expression('value', degree=0, value=self.A0(L))
            bc_outlet = DirichletBC(self.V2.sub(0), self._A_out, outlet_bdry)
        else:
            self._U_out = Expression(('A', 'q'), degree=0,
                                     A=self.A0(L), q=self.q0)
            bc_outlet = DirichletBC(self.V2, self._U_out, outlet_bdry)

        self.bcs = [bc_inlet, bc_outlet]

        # Terms for variational form
        U_v_dx = A*v1*dx + q*v2*dx

        Un_v_dx = self.Un[0]*v1*dx + self.Un[1]*v2*dx


        F2_v2_ds = (pow(q, 2)/(A+DOLFIN_EPS)\
                   +self.f*sqrt(self.A0*(A+DOLFIN_EPS)))*v2*ds

        F2_dv2_dx = (pow(q, 2)/(A+DOLFIN_EPS)\
                    +self.f*sqrt(self.A0*(A+DOLFIN_EPS)))*grad(v2)[0]*dx

        dF_v_dx = grad(q)[0]*v1*dx + F2_v2_ds - F2_dv2_dx


        Fn_v_ds = (pow(self.Un[1], 2)/(self.Un[0])\
                  +self.f*sqrt(self.A0*(self.Un[0])))*v2*ds

        Fn_dv_dx = (pow(self.Un[1], 2)/(self.Un[0])\
                   +self.f*sqrt(self.A0*(self.Un[0])))*grad(v2)[0]*dx

        dFn_v_dx = grad(self.Un[1])[0]*v1*dx + Fn_v_ds - Fn_dv_dx


        S_v_dx = - 2*sqrt(pi)/self.db/Re*q/sqrt(A+DOLFIN_EPS)*v2*dx\
               + (2*sqrt(A+DOLFIN_EPS)*(sqrt(pi)*self.f
                                       +sqrt(self.A0)*self.dfdr)\
                 -(A+DOLFIN_EPS)*self.dfdr)*self.drdx*v2*dx

        Sn_v_dx = -2*sqrt(pi)/self.db/Re*self.Un[1]/sqrt(self.Un[0])*v2*dx\
                + (2*sqrt(self.Un[0])*(sqrt(pi)*self.f+sqrt(self.A0)*self.dfdr)\
                  -(self.Un[0])*self.dfdr)*self.drdx*v2*dx

        # Variational form
        self.variational_form = U_v_dx\
            - Un_v_dx\
            + self.dt*self.theta*dF_v_dx\
            + self.dt*(1-self.theta)*dFn_v_dx\
            - self.dt*self.theta*S_v_dx\
            - self.dt*(1-self.theta)*Sn_v_dx


    def solve(self):
        """
        Calls FEniCS's solve() function for the variational form.
        """
        F = self.variational_form
        J = derivative(F, self.U)
        prob = NonlinearVariationalProblem(F, self.U, self.bcs, J=J)
        sol = NonlinearVariationalSolver(prob)
        sol.parameters['newton_solver']['linear_solver'] = 'mumps'
        sol.parameters['newton_solver']['lu_solver']['reuse_factorization'] = True
        sol.parameters['newton_solver']['maximum_iterations'] = 1000
        sol.solve()


    def update_solution(self):
        """
        Stores current solution U as previous solution Un.
        """
        self.Un.assign(self.U)


    def update_pressure(self):
        """
        Calculates pressure.
        """
        self.pn.assign(Expression('p0 + f*(1-sqrt(A0/A))', degree=2,
                                    p0=self.param['p0'],
                                    f=self.f, A0=self.A0,
                                    A=self.Un.split(True)[0]))


    def compute_pressure(self, f, A0, A):
        """
        Computes pressure.

        Arguments
        ---------
        f : numpy.array
            Elasticity relation
        A0: numpy.array
            Area at rest
	    A : numpy.array
            Area

        Returns
        -------
        return : float
            Pressure at the outlet
        """
        return self.p0 + f*(1-np.sqrt(A0/A))


    def compute_outlet_pressure(self, A):
        """
        Computes pressure at the outlet.

        Arguments
        ---------
        A : float
            Area value at the outlet

        Returns
        -------
        return : float
            Pressure at the outlet
        """
        L = self.param['L']
        p0 = self.param['p0']
        return p0 + self.f(L)*(1-np.sqrt(self.A0(L)/A))


    def CFL_term(self, x, A, q):
        """
        Computes the CFL number.

        Arguments
        ---------
        x : float
            Point at which the condition is to be checked
        A : float
            Area value at x
        q : float
            Flow rate value at x

        Returns
        -------
        return : float
            CFL number
        """
        return 1/np.abs(q/A+np.sqrt(self.f(x)/2/self.param['rho']\
                                   *np.sqrt(self.A0(x)/A)))


    def check_CFL(self, x, A, q):
        """
        Checks the CFL condition.

        Arguments
        ---------
        x : float
            Point at which the condition is to be checked
        A : float
            Area value at x
        q : float
            Flow rate value at x

        Returns
        -------
        return : boolean
            True if the CFL condition is fulfilled
        """
        return self.dt/self.dex < self.CFL_term(x, A, q)


    def adjust_dex(self, x, A, q, margin=0.05):
        """
        Adjusts spatial step at a bifurcation to respect the CFL condition.

        Arguments
        ---------
        x : float
            Point at which the condition is to be checked
        A : float
            Area value at x
        q : float
            Flow rate value at x
        margin : float
            Margin of CFL number
        """
        M = self.CFL_term(x, A, q)
        self.dex = (1+margin)*self.dt/M


    @property
    def q_in(self):
        """
        Inlet flow rate (only for a root artery)
        """
        return self._q_in.value

    @q_in.setter
    def q_in(self, value):
        self._q_in.value = value


    @property
    def U_in(self):
        """
        Inlet boundary conditions (only for non-root arteries (daughters))
        """
        return np.array([self._U_in.A, self._U_in.q])

    @U_in.setter
    def U_in(self, U):
        self._U_in.A = U[0]
        self._U_in.q = U[1]


    @property
    def A_out(self):
        """
        Outlet area (only in use for end arteries)
        """
        return self._A_out.value

    @A_out.setter
    def A_out(self, value):
        self._A_out.value = value


    @property
    def U_out(self):
        """
        Outlet boundary conditions (only for parent arteries)
        """
        return np.array([self._U_out.A, self._U_out.q])

    @U_out.setter
    def U_out(self, U):
        self._U_out.A = U[0]
        self._U_out.q = U[1]
