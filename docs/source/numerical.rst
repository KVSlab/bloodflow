Numerical implementation
========================

This page explains the numerical implementation of the equations derived in :ref:`theory`. A solution is calculated for each vessel/artery within the network separately, and the relationship between parent and daughter vessels is governed via boundary conditions. Thus, time advancement is handled by the class `ArteryNetwork`, while the computation of a solution per time step is handled by the class `Artery`.

Time discretisation
-------------------

Time is discretised using a finite difference :math:`\theta`-rule based algorithm, that is given the previous solution :math:`U^n` one can find the next solution :math:`U^{n+1}` using

.. math::
  U^{n+1} = U^n + k \left[ \theta U^{n+1} + (1-\theta) U^{n} \right]

with :math:`k = \delta t` and :math:`0 \leq \theta \leq 1`, where :math:`\theta = 0.5` is used by default and corresponds to the Crank-Nicolson method, which is second order and unconditionally stable. The method is implemented in :meth:`arteryfe.Artery.define_solution`::

  # Crank-Nicolson parameter
  self.theta = theta

  # Variational form
  self.variational_form = U_v_dx\
        - Un_v_dx\
        + self.dt*self.theta*dF_v_dx\
        + self.dt*(1-self.theta)*dFn_v_dx\
        - self.dt*self.theta*S_v_dx\
        - self.dt*(1-self.theta)*Sn_v_dx

Finite element discretisation
-----------------------------

Arteries are discretised using the finite element (FE) method implemented in FEniCS_. Because we model arteries in 1D using their cross-sectional area along the longitudinal axis, we can use FEniCS built-in interval mesh to initialise the geometry using :meth:`arteryfe.Artery.define_geometry`::

  self.mesh = IntervalMesh(self.Nx, 0, self.L)
  self.elV = FiniteElement('CG', self.mesh.ufl_cell(), 1)
  self.V = FunctionSpace(self.mesh, self.elV)
  self.V2 = FunctionSpace(self.mesh, self.elV*self.elV)

  # Initial vessel-radius and deduced quantities
  self.r0 = Expression('Ru*pow(Rd/Ru, x[0]/L)',
                         degree=2, Ru=self.Ru, Rd=self.Rd, L=self.L)
  self.A0 = Expression('pi*pow(r0, 2)', degree=2, r0=self.r0)
  self.f = Expression('4.0/3.0*(k1*exp(k2*r0) + k3)', degree=2,
                        k1=self.k1, k2=self.k2, k3=self.k3, r0=self.r0)
  self.dfdr = Expression('4.0/3.0*k1*k2*exp(k2*r0)', degree=2,
                           k1=self.k1, k2=self.k2, r0=self.r0)
  self.drdx = Expression('logRdRu/L*r0', degree=2,
                           logRdRu=np.log(self.Rd/self.Ru), L=self.L,
                      r0=self.r0)

The FE method requires rewriting the system of equations in its variational form. For a basic introduction to deriving the variational form of a system of equations you can refer to the `FEniCS tutorial`_. Recall from :ref:`theory` that the final system of equations is

.. _FEniCS: https://fenicsproject.org/documentation/
.. _FEniCS tutorial: https://fenicsproject.org/pub/tutorial/html/._ftut1004.html#ch:poisson0:varform

.. math::
  \begin{split}
  &\dfrac{\partial}{\partial t} \begin{pmatrix} A(z,t) \\ q(z,t) \end{pmatrix} + \dfrac{\partial}{\partial z} \begin{pmatrix} q(z,t)\\ \dfrac{q(z,t)^2}{A(z,t)} + f(r_0) \sqrt{A_0(z) A(z,t)} \end{pmatrix} =\\
  &\begin{pmatrix} 0 \\ -\dfrac{2 \pi R(z,t)}{\delta_b \mathcal{Re}} \dfrac{q(z,t)}{A(z,t)} +\left( 2 \sqrt{A(z,t)} \left( \sqrt{\pi} f(r_0) + \sqrt{A_0(z)} \frac{df(r_0)}{dr_0 } \right) - A(z,t) \dfrac{df(r_0)}{dr_0} \right) \dfrac{dr_0(z)}{dz} \end{pmatrix},
  \end{split}

which is a conservation system of equations.

.. math::
  \dfrac{\partial}{\partial t} \boldsymbol{U} + \dfrac{\partial}{\partial z} \boldsymbol{F} =
  \boldsymbol{S}.

The terms of the governing equations are hence implemented in terms of the vectors :math:`\boldsymbol{U}, \boldsymbol{F}, \boldsymbol{S}` in :meth:`arteryfe.ArteryNetwork.flux`::

  def flux(self, a, U, x):
    return np.array([U[1], U[1]**2 + a.f(x)*np.sqrt(a.A0(x)*U[0])])

and :meth:`arteryfe.ArteryNetwork.source`::

  def source(self, a, U, x):
    S1 = 0
    S2 = -2*np.sqrt(np.pi)/a.db/a.Re*U[1]/np.sqrt(U[0])\
        + (2*np.sqrt(U[0])*(np.sqrt(np.pi)*a.f(x)\
                            +np.sqrt(a.A0(x))*a.dfdr(x))\
           -U[0]*a.dfdr(x))*a.drdx(x)
    return np.array([S1, S2])

The variational form is implemented in :meth:`arteryfe.Artery.define_solution` as::

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
  self.pn.assign(Expression('p0', degree=2, p0=self.p0))

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
  S_v_dx = - 2*sqrt(pi)/self.db/self.Re*q/sqrt(A+DOLFIN_EPS)*v2*dx\
           + (2*sqrt(A+DOLFIN_EPS)*(sqrt(pi)*self.f
                                   +sqrt(self.A0)*self.dfdr)\
             -(A+DOLFIN_EPS)*self.dfdr)*self.drdx*v2*dx
  Sn_v_dx = -2*sqrt(pi)/self.db/self.Re*self.Un[1]/sqrt(self.Un[0])*v2*dx\
            + (2*sqrt(self.Un[0])*(sqrt(pi)*self.f+sqrt(self.A0)*self.dfdr)\
              -(self.Un[0])*self.dfdr)*self.drdx*v2*dx

  # Variational form
  self.variational_form = U_v_dx\
        - Un_v_dx\
        + self.dt*self.theta*dF_v_dx\
        + self.dt*(1-self.theta)*dFn_v_dx\
        - self.dt*self.theta*S_v_dx\
        - self.dt*(1-self.theta)*Sn_v_dx

The variable `self.variational_form` is solved in :meth:`arteryfe.Artery.solve` using a nonlinear variational solver from FEniCS_::

  F = self.variational_form
  J = derivative(F, self.U)
  solve(F == 0, self.U, self.bcs, J=J)

Boundary conditions
-------------------

We prescribe the flow rate directly at the inlet of the root vessel, which is implemented in :meth:`arteryfe.ArteryNetwork.set_bcs` using the inlet file provided in the .cfg file::

  # Update inlet boundary conditions
  self.arteries[0].q_in = q_in

At the outlet of the terminal vessels a three-element Windkessel model is applied. This type of model is also called a lumped model and uses an electric circuit analog with specific resistance and compliance parameters to represent the downstream artery tree. Because the Windkessel model provides an estimate for pressure instead of cross-sectional area the boundary condition cannot be calculated directly. A fixed-point iterative scheme is implemented in :meth:`arteryfe.ArteryNetwork.compute_A_out` with an initial guess for the outlet pressure::

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

The boundary conditions at a bifurcation are somewhat more complex. Three arteries are involved in a bifurcation, for which the current and next time step for each of the three variables need to be calculated. Thus, we arrive at a system of 18 equations for 18 variables. The solution of this system is implemented using Newton's method in :meth:`arteryfe.ArteryNetwork.newton`::

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
