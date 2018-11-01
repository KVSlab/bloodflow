.. _theory:

Blood flow dynamics in 1D
=========================

artery.fe implements the 1D system of equations derived by [Olufsen:2000] that is commonly used to numerically model blood flow dynamics. We assume that the reader is familiar with the general methods of modelling laminar fluid flow in 1D and provide a derivation for completeness.

Nomenclature
------------

We use the nomenclature listed in the table below.

===================================  ====================
Symbol                               Interpretation
===================================  ====================
:math:`r`                            radial direction
:math:`z`                            axial direction
:math:`t`                            time
:math:`\boldsymbol{u} = (u_z, u_r)`  velocity
:math:`A`                            cross-sectional area
:math:`R`                            radius
:math:`q`                            flow rate
:math:`p`                            pressure
:math:`\rho`                         density
:math:`\nu`                          viscosity
:math:`Q`                            characteristic flow rate
:math:`T`                            cardiac cycle length
:math:`\delta`                       artery boundary layer
:math:`\bar{\cdot}`                  average
:math:`E`                            Young's modulus
:math:`h`                            artery wall thickness
:math:`k_i`                          elastic parameters
:math:`f`                            elastic relation
:math:`\mathcal{Re}`                 Reynold's number
:math:`\Delta t`                     discrete time step size
:math:`\Delta x`                     discrete spatial step size
:math:`m`                            grid location
:math:`M`                            outlet grid location
:math:`\mathcal{M}`                  bifurcation grid location, that is :math:`M` for the parent vessel and 0 for its daughter vessels
:math:`n`                            time point
:math:`R_i`                          resistance parameters
:math:`C`                            compliance parameters
===================================  ====================

Governing equations
-------------------

Consider an arterial segment that we model as an axisymmetric tube in a cylindrical coordinate system with radial direction :math:`r` and axial direction :math:`z`. Therefore, the governing equation reduces to

.. math::

  \frac{\partial u_z(r,z,t)}{\partial z} + \frac{1}{r} \frac{\partial(ru_r(r,z,t))}{\partial r} = 0 \qquad (1),

where :math:`\boldsymbol{u} = (u_z(r,z,t), u_r(r,z,t))`. Integration of (1) over the cross-sectional area with :math:`A(z,t) = \pi R(z,t)^2` yields

.. math::
  2 \pi \int_0^{R(z,t)} \left( \frac{\partial u_z(r,z,t)}{\partial z} + \frac{1}{r} \frac{\partial(ru_r(r,z,t))}{\partial r} \right) r dr = 0,\\
  2 \pi \int_0^{R(z,t)} \frac{\partial u_z(r,z,t)}{\partial z} r dr + 2 \pi \left[ r u_r(r,z,t) \right]_0^{R(z,t)} = 0 \qquad (2)

where :math:`R(z,t)` describes the vessel radius.

The first term of (2) can be evaluated using Leibniz' integral rule ("differentiation under the integral sign"), which states

.. math::

  \frac{\partial}{\partial x} \int_{a(x)}^{b(x)} f(x, y) dy = \int_{a(x)}^{b(x)} \frac{\partial f(x,y)}{\partial x} dy + f(x,b(x)) \frac{\partial b(x)}{\partial x} - f(x,a(x)) \frac{\partial a(x)}{\partial x}.

Applying this rule to (2) results in

.. math::
  2 \pi \frac{\partial}{\partial z} \int_0^{R(z,t)} u_z(r,z,t) r dr - 2 \pi \frac{\partial R(z,t)}{\partial z}\left[ r u_z(r,z,t) \right]_{R(z,t)} + 2 \pi \left[ r u_r(r,z,t) \right]_0^{R(z,t)} = 0.

Defining flux through the vessel as

.. math::
  q(r,z,t) = 2\pi \int_0^{R(z,t)} u_z(r,z,t) r dr


yields

.. math::
  \frac{\partial q(r,z,t)}{\partial z} - 2 \pi \frac{\partial R(z,t)}{\partial z}\left[ r u_z(r,z,t) \right]_{R(z,t)} + 2 \pi \left[ r u_r(r,z,t) \right]_0^{R(z,t)} = 0. \qquad (3)

Due to no-slip :math:`u_z` vanishes at the boundary :math:`r = R(z,t)`

.. math::
  \left. u_z(r,z,t) \right|_{R(z,t)} = 0. \qquad (4)

An arterial segment is embedded within an arterial tree and thus stretch along the :math:`z`-direction is restricted (tether). Thus

.. math::

  \left. u_r(r,z,t) \right|_{R(z,t)} = \frac{\partial R(z,t)}{\partial t}, \qquad(5)

which allows us to write the continuity equation (1) in terms of flow rate :math:`q` and cross-sectional area :math:`A`

.. math::

  \frac{\partial q(r,z,t)}{\partial z} + \frac{\partial A(z,t)}{\partial t} = 0. \qquad (6)

We treat the momentum equation in the same coordinate system in a similar fashion. For Poiseuille flow it reads

.. math::

  \begin{split}
    \frac{\partial u_z(r,z,t)}{\partial t} + u_z(r,z,t) \frac{\partial u_z(r,z,t)}{\partial z} + u_r(r,z,t) \frac{\partial u_z(r,z,t)}{\partial r} +& \frac{1}{\rho} \frac{\partial p(z,t)}{\partial z} =\\
    & \frac{\nu}{r} \frac{\partial}{\partial r} \left( r \frac{\partial u_z(r,z,t)}{\partial r} \right), \qquad (7)
  \end{split}

where :math:`p(z,t)` denotes pressure and :math:`\nu` kinematic viscosity. Nondimensionalisation of (7) shows that the longitudinal viscous term :math:`\nu \partial^2 u_z / \partial z^2` is much smaller than the radial viscous term, due to the much longer length scale of arteries compared to the radius, and was therefore eliminated. Integrating over cross-sectional area, whilst recognising that :math:`p(z,t)` is constant over this area, yields

.. math::

  \begin{split}
  2\pi \int_0^{R(z,t)} \frac{\partial u_z(r,z,t)}{\partial t} r dr + 2\pi \int_0^{R(z,t)} u_z(r,z,t) \frac{\partial u_z(r,z,t)}{\partial z} r dr &\\
  + 2\pi \int_0^{R(z,t)} u_r(r,z,t) \frac{\partial u_z(r,z,t)}{\partial r} r dr + \frac{A(z,t)}{\rho} \frac{\partial p(z,t)}{\partial z} & =\\
  2 \pi \nu R(z,t) \frac{\partial u_z(r,z,t)}{\partial r} & \left. \right|_{r = R(z,t)}. \qquad (8)
  \end{split}

Application of Leibniz' integral rule to the first term on the left-hand side (LHS) in (8) and using (4) gives

.. math::

  \frac{\partial}{\partial t} \int_0^{R(z,t)} u_z(r,z,t) r dr = \frac{\partial}{\partial t} \int_0^{R(z,t)} u_z(r,z,t) r dr - \frac{\partial R(z,t)}{\partial t} \left[ u_z(r,z,t) r \right]_{R(z,t)} = \frac{\partial q(r,z,t)}{\partial t}.

Integration by parts of the third LHS term of (8) results in

.. math::

  \begin{split}
  2\pi \int_0^{R(z,t)} u_r(r,z,t) & \frac{\partial u_z(r,z,t)}{\partial r} r dr =\\
  2\pi &[u_r(r,z,t) u_z(r,z,t) r]_0^{R(z,t)} - 2\pi \int_{R(z,t)} u_z(r,z,t) \frac{\partial r u_r(r,z,t)}{\partial r} dr
  \end{split}

and using (1) and (4) leads to

.. math::

  \begin{split}
  2\pi \int_0^{R(z,t)} u_r(r,z,t) \frac{\partial u_z(r,z,t)}{\partial r} r dr = 2 \pi \int_0^{R(z,t)} u_z(r,z,t) & \frac{\partial u_z(r,z,t)}{\partial z} r dr =\\
  \pi & \int_0^{R(z,t)}\frac{\partial u_z(r,z,t)^2}{\partial z} r dr.
  \end{split}

Using these results in (8) gives

.. math::

  \begin{split}
  \frac{\partial q(r,z,t)}{\partial t} + 2\pi \frac{\partial}{\partial z} \int_0^{R(z,t)} u_z(r,z,t)^2 r dr + \frac{A(z,t)}{\rho} \frac{\partial p(z,t)}{\partial z} &=\\
  2\pi & \nu R(z,t) \left. \frac{\partial u_z(r,z,t)}{\partial r} \right|_{R(z,t)}. \qquad (9)
  \end{split}

To solve the remaining terms it is necessary to make assumptions about the velocity profile of blood flow through an artery. Blood flow is considered positively pulsatile and laminar, and vessels can be considered slightly tapered, therefore the velocity profile is assumed to be mostly flat with a thin boundary layer with cardiac cycle length :math:`T` and width :math:`\delta = (\nu T / (2\pi))^{0.5}`, such that :math:`\delta \ll R(z,t)`. The axial velocity :math:`u_z(r,z,t)` thus has the form

.. math::

  u_z(r,z,t) = \begin{cases}
  \bar{u}_z(z,t) & r \leq R(z,t)-\delta\\
  \bar{u}_z(z,t) (R(z,t)-r)/\delta & R(z,t)-\delta < r \leq R(z,t),
  \end{cases} \qquad (10)

where :math:`\bar{u}_z(z,t)` is the average axial velocity outside the boundary layer. This leads to a flat velocity profile outside the boundary layer and linearly increasing profile (from 0 to :math:`\bar{u}_z(z,t)`) inside the boundary layer. Note that a physiological cardiac cycle at rest has between 60 and 70 beats per minute (0.6 s :math:`\leq T \leq` 1.1 s), therefore the boundary layer is 0.07--0.09 cm in size. This is much smaller than the minimal inlet radius of arteries considered in this work, namely 0.14 cm, and therefore (10) is appropriate for the desired velocity profile. The first and second terms of (9) can then be expressed as a power series in :math:`\delta`

.. math::

  q = 2\pi \int_0^{R(z,t)} u_z(r,z,t) r dr = A \bar{u}_z(z,t) \left( 1 - \frac{\delta}{R(z,t)} + \mathcal{O}(\delta^2) \right),\\
  2\pi \frac{\partial}{\partial z} \int_0^{R(z,t)} u_z(r,z,t)^2 r dr = A \bar{u}_z(z,t) \left( 1 - \frac{4}{3} \frac{\delta}{R(z,t)} + \mathcal{O}(\delta^2) \right).

Using these solutions the second term of (9) becomes

.. math::

  2\pi \frac{\partial}{\partial z} \int_0^{R(z,t)} u_z(r,z,t)^2 r dr = \frac{q(z,t)^2}{A(z,t)} \left( 1 + \frac{2}{3} \frac{\delta}{R(z,t)} + \mathcal{O}(\delta^2) \right).

This leaves the term on the right-hand side (RHS) of (9) to be evaluated using the velocity profile

.. math::

  2 \pi \nu R(z,t) \frac{\partial u_z(r,z,t)}{\partial r} = - \frac{2 \pi \nu R(z,t)}{\delta} \frac{q(z,t)}{A(z,t)} + \mathcal{O}(\delta)

such that finally, keeping only leading order terms in :math:`\delta`, the momentum equation reads

.. math::

  \frac{\partial q(z,t)}{\partial t} + \frac{\partial}{\partial z} \left( \frac{q(z,t)^2}{A(z,t)} \right) + \frac{A(z,t)}{\rho} \frac{\partial p(z,t)}{\partial z} = - \frac{2 \pi \nu R(z,t)}{\delta} \frac{q(z,t)}{A(z,t)}. \qquad (11)

In order to solve the system of (6) and (11) they need to be written in conservation form

.. math::

  \frac{\partial \boldsymbol{U}}{\partial t} + \frac{\partial \boldsymbol{F}}{\partial z} = \boldsymbol{S}. \qquad (12)

The quantity :math:`B` is introduced and chosen to fulfill

.. math::
  B(r_0(z), p(z,t)) = \frac{1}{\rho} \int A(z,t) dp(z,t),

with :math:`r_0(z)` initial radius at rest such that

.. math::

  \frac{\partial B(r_0(z), p(z,t))}{\partial z} = \frac{A}{\rho} \frac{\partial p(z,t)}{\partial z} + \frac{\partial B(r_0(z), p(z,t))}{\partial r_0(z)} \frac{\partial r_0(z)}{\partial z}

Then, adding the term :math:`(\partial B / \partial r_0) (\partial r_0 / \partial z)` to both sides of (11), the system of equations can be written in conservation form

.. math::

  \begin{split}
  \dfrac{\partial}{\partial t} \begin{pmatrix} A(z,t) \\ q(z,t) \end{pmatrix} + \dfrac{\partial}{\partial z} & \begin{pmatrix} q(z,t)\\ \dfrac{q(z,t)^2}{A(z,t)} + B(r_0(z), p(z,t)) \end{pmatrix} =\\ & \qquad \begin{pmatrix} 0 \\ - \dfrac{2 \pi \nu R(z,t)}{\delta} \dfrac{q(z,t)}{A(z,t)} + \dfrac{\partial B(r_0(z), p(z,t))}{\partial r_0(z)} \dfrac{\partial r_0(z)}{\partial z} \end{pmatrix}. \qquad (13)
  \end{split}

Currently, (13) contains three unknowns (:math:`q, A, p`) for two equations, thus a third relation is needed to solve the system of equations. The aforementioned equation, referred to as the state equation, describes the relationship between :math:`A(z,t)` and :math:`p(z,t)`. One choice for the state equation is the linearly elastic relation

.. math::

  p(z,t) - p_0 = \frac{4}{3} \frac{Eh}{r_0(z)} \left( 1 - \sqrt{\frac{A_0(z)}{A(z,t)}} \right) \qquad (14),

where the constant :math:`p_0` is the diastolic pressure, :math:`E` is the Young's modulus of the vessel wall, :math:`h` is the wall width and :math:`A_0(z) = \pi r_0(z)^2`. The relationship :math:`Eh/r_0` is based on compliance estimates

.. math::

  \frac{Eh}{r_0(z)} = k_1 \exp (k_2 r_0(z)) + k_3, \qquad (15)

with :math:`k_1, k_2, k_3` as constants. Using (14) and defining :math:`f(r_0) = 4Eh/(3r_0)` the quantities :math:`B(r_0, p), (\partial B / \partial r_0) (\partial r_0 / \partial z)` can be evaluated from (13)

.. math::

  B(r_0(z), p(z,t)) = \frac{1}{\rho} \int \frac{f(r_0) A_0(r_0)}{p(z,t)^2/f(r_0) - 2p(z,t) + f(r_0)} dp = \frac{1}{\rho} \frac{f(r_0) A_0(r_0)}{(1 - p(z,t)/f(r_0))},\\
  \begin{split}
  \frac{\partial B(r_0(z), p(z,t))}{\partial r_0(z)} \frac{\partial r_0(z)}{\partial z} &=\\
  &\frac{1}{\rho} \left( 2 \sqrt{A(r_0)} \left( \sqrt{\pi} f(r_0) + \sqrt{A_0} \frac{df(r_0)}{dr_0 } \right) - A(r_0) \frac{df(r_0)}{dr_0} \right) \frac{dr_0}{dz},
  \end{split}

thus, (13) becomes

.. math::

  \begin{split}
  &\dfrac{\partial}{\partial t} \begin{pmatrix} A(z,t) \\ q(z,t) \end{pmatrix} + \dfrac{\partial}{\partial z} \begin{pmatrix} q(z,t)\\ \dfrac{q(z,t)^2}{A(z,t)} + \frac{f(r_0)}{\rho} \sqrt{A_0(z) A(z,t)} \end{pmatrix} =\\
  &\begin{pmatrix} 0 \\ -\dfrac{2 \pi \nu q(z,t) R(z,t)}{\delta A(z,t)} + \dfrac{1}{\rho} \left( 2 \sqrt{A(z,t)} \left( \sqrt{\pi} f(r_0) + \sqrt{A_0(z)} \frac{df(r_0)}{dr_0 } \right) - A(z,t) \dfrac{df(r_0)}{dr_0} \right) \dfrac{dr_0(z)}{dz} \end{pmatrix}. \qquad (16)
  \end{split}

To nondimensionalise we choose appropriate scaling parameters for the system variables:

================================  ====================
Variable                          Physical meaning
================================  ====================
:math:`z \sim R`                  length scale
:math:`r_0(z) \sim R`             radius at rest
:math:`q(z,t) \sim Q`             flow rate
:math:`t \sim R^3/Q`              time
:math:`A(z,t) \sim R^2`           cross-sectional area
:math:`p(z,t) \sim \rho Q^2/R^4`  pressure
================================  ====================

The resulting dimensionless system of equations is

.. math::

  \begin{split}
  &\dfrac{\partial}{\partial t} \begin{pmatrix} A(z,t) \\ q(z,t) \end{pmatrix} + \dfrac{\partial}{\partial z} \begin{pmatrix} q(z,t)\\ \dfrac{q(z,t)^2}{A(z,t)} + f(r_0) \sqrt{A_0(z) A(z,t)} \end{pmatrix} =\\
  &\begin{pmatrix} 0 \\ -\dfrac{2 \pi R(z,t)}{\delta \mathcal{Re}} \dfrac{q(z,t)}{A(z,t)} +\left( 2 \sqrt{A(z,t)} \left( \sqrt{\pi} f(r_0) + \sqrt{A_0(z)} \frac{df(r_0)}{dr_0 } \right) - A(z,t) \dfrac{df(r_0)}{dr_0} \right) \dfrac{dr_0(z)}{dz} \end{pmatrix}. \qquad (17)
  \end{split}

Boundary conditions
-------------------

Boundary conditions are applied at both ends of each vessel and are either an inlet, outlet or bifurcation condition.

Inlet
^^^^^

The inlet boundary condition only used at the inlet of the parent vessel. For a given :math:`q_0^{n+1}` :math:`A_0^{n+1}` is calculated as

.. math::

  A_0^{n+1} = A_0^n - \frac{\Delta t}{\Delta z} \left( q_{1/2}^{n+1/2} - q_{-1/2}^{n+1/2} \right), \qquad (18)

where :math:`q_{-1/2}^{n+1/2}` can be evaluated using

.. math::

  q_0^{n+1/2} = (q_{1/2}^{n+1/2} + q_{-1/2}^{n+1/2})/2 \qquad (19)

with :math:`q_0^{n+1/2}` evaluated directly from the inlet flux function and :math:`q_{1/2}^{n+1/2}`, evaluated from the Lax-Wendroff approximation

.. math::

  \boldsymbol{U}_j^{n+1/2} = \frac{\boldsymbol{U}_{j+1/2}^n + \boldsymbol{U}_{j-1/2}^n}{2} + \frac{\Delta t}{2} \left( - \frac{\boldsymbol{F}_{j+1/2}^n - \boldsymbol{F}_{j-1/2}^n}{\Delta z} + \frac{\boldsymbol{S}_{j+1/2}^n + \boldsymbol{S}_{j-1/2}^n}{2} \right) \qquad (20)

Outlet
^^^^^^

The outlet boundary condition is a three-element Windkessel (3WK), which is given by

.. math::

  \frac{\partial p(z,t)}{\partial t} = R_1 \frac{\partial q(z,t)}{\partial t} - \frac{p(z,t)}{R_2 C} + \frac{q(z,t) (R_1 + R_2)}{R_2 C}.

The 3WK model uses an electrical circuit analog representation of the downstream arterial tree, where electrical current represents :math:`q` and voltage represents :math:`p`, using resistance (:math:`R_i`) and compliance (:math:`C`) parameters . Discretisation yields

.. math::

  \frac{p_m^{n+1} - p_m^n}{\Delta t} = R_1 \frac{q_m^{n+1} - q_m^n}{\Delta t} - \frac{p_m^n}{R_2 C_T} + \frac{q_m^n (R_1 + R_2)}{R_2 C_T}, \qquad (21)

which is used as the outlet boundary condition. Solutions for :math:`A_m^{n+1}` and the discretised state equation

.. math::

  p_m^{n+1} = \frac{4}{3} \frac{E h}{(r_0)_m} \left( 1 - \sqrt{\frac{(A_0)_m}{A_m^{n+1}}} \right)

are found using an iterative scheme, starting with an initial guess for :math:`p_m^{n+1}`. Then, :math:`q_m^{n+1}` can be evaluated using (21). Using

.. math::

  A_m^{n+1} = A_m^n - \frac{\Delta t}{\Delta z} \left( q_{m+1/2}^{n+1/2} - q_{m-1/2}^{n+1/2} \right)

the next iteration of :math:`p_m^{n+1}` can then be calculated until the difference between two iterations has dropped below a threshold value.

Bifurcation
^^^^^^^^^^^

Lastly, bifurcation boundary conditions apply between a parent vessel p and two daughter vessels d1 and d2. Conservation of flow implies

.. math::

  \left( q^{(p)} \right)_M^n = \left( q^{(d1)} \right)_0^n + \left( q^{(d2)} \right)_0^n \qquad (22)

and continuity of pressure yields

.. math::

  \left( p^{(p)} \right)_M^n = \left( p^{(d1)} \right)_0^n = \left( p^{(d2)} \right)_0^n. \qquad (23)

Written in terms of A (23) becomes

.. math::

  \left( f^{(p)} \right)_M \left( 1 - \sqrt{\frac{\left( A_0^{(p)} \right)_M}{\left( A^{(p)} \right)_M^n}} \right) = \left( f^{(d1)} \right)_0 \left( 1 - \sqrt{\frac{\left( A_0^{(d1)} \right)_0}{\left( A^{(d1)} \right)_0^n}} \right), \qquad (24)\\
  \left( f^{(p)} \right)_M \left( 1 - \sqrt{\frac{\left( A_0^{(p)} \right)_M}{\left( A^{(p)} \right)_M^n}} \right) = \left( f^{(d2)} \right)_0 \left( 1 - \sqrt{\frac{\left( A_0^{(d2)} \right)_0}{\left( A^{(d2)} \right)_0^n}} \right). \qquad (25)

On both sides of the boundary q and A are calculated from the Lax-Wendroff discretisation

.. math::

  \left( A^{(i)} \right)_{\mathcal{M}}^{n+1} = \left( A^{(i)} \right)_{\mathcal{M}}^n - \frac{\Delta t}{\Delta z} \left(\left( F_1^{(i)} \right)_{\mathcal{M}+1/2}^{n+1/2} - \left( F_1^{(i)} \right)_{\mathcal{M}-1/2}^{n+1/2} \right) \qquad (26)\\
  \begin{split}
  \left( q^{(i)} \right)_{\mathcal{M}}^{n+1} = \left( q^{(i)} \right)_{\mathcal{M}}^n - \frac{\Delta t}{\Delta z} \left(\left( F_2^{(i)} \right)_{\mathcal{M}+1/2}^{n+1/2} - \right.&\left. \left( F_2^{(i)} \right)_{\mathcal{M}-1/2}^{n+1/2} \right) +\\
  &\frac{\Delta t}{2} \left(\left( S_2^{(i)} \right)_{\mathcal{M}+1/2}^{n+1/2} + \left( S_2^{(i)} \right)_{\mathcal{M}-1/2}^{n+1/2} \right), \qquad (27)
  \end{split}

where :math:`i = p, d1, d2` and :math:`\mathcal{M} = M` if :math:`i = p` and :math:`\mathcal{M} = 0` otherwise. We introduce the ghost points :math:`q_{M+1/2}^{n+1/2}` and :math:`A_{M+1/2}^{n+1/2}`, which are not part of the geometry of the parent vessel, but lie beyond the outlet point. These are, analogously to the inlet boundary condition, evaluated using

.. math:

  q_M^{n+1/2} = \frac{q_{M-1/2}^{n+1/2} + q_{M+1/2}^{n+1/2}}{2}, \qquad (28)\\
  A_M^{n+1/2} = \frac{A_{M-1/2}^{n+1/2} + A_{M+1/2}^{n+1/2}}{2}. \qquad (29)

(22)--(29) defines a system of eighteen equations for eighteen unknowns

.. math ::

  x_1 = \left( q^{(p)} \right)_M^{n+1} \qquad x_2 = \left( q^{(p)} \right)_M^{n+1/2} \qquad x_3 = \left( q^{(p)} \right)_{M+1/2}^{n+1/2}\\
  x_4 = \left( q^{(d1)} \right)_0^{n+1} \qquad x_5 = \left( q^{(d1)} \right)_0^{n+1/2} \qquad x_6 = \left( q^{(d1)} \right)_{-1/2}^{n+1/2}\\
  x_7 = \left( q^{(d1)} \right)_0^{n+1} \qquad x_8 = \left( q^{(d1)} \right)_0^{n+1/2} \qquad x_9 = \left( q^{(d1)} \right)_{-1/2}^{n+1/2}\\
  x_{10} = \left( A^{(p)} \right)_M^{n+1} \qquad x_{11} = \left( A^{(p)} \right)_M^{n+1/2} \qquad x_{12} = \left( A^{(p)} \right)_{M+1/2}^{n+1/2}\\
  x_{13} = \left( A^{(d1)} \right)_0^{n+1} \qquad x_{14} = \left( A^{(d1)} \right)_0^{n+1/2} \qquad x_{15} = \left( A^{(d1)} \right)_{-1/2}^{n+1/2}\\
  x_{16} = \left( A^{(d1)} \right)_0^{n+1} \qquad x_{17} = \left( A^{(d1)} \right)_0^{n+1/2} \qquad x_{18} = \left( A^{(d1)} \right)_{-1/2}^{n+1/2}.

The system of equations can be solved using Newton's method

.. math::

  \boldsymbol{x}_{k+1} = \boldsymbol{x}_k - \left( \boldsymbol{J}(\boldsymbol{x}_k) \right)^{-1} \boldsymbol{f_J}(\boldsymbol{x}_k) \text{ for } k = 0, 1, 2, \ldots,

where k indicates the current iteration, :math:`\boldsymbol{x} = (x_1, x_2, \ldots, x_{18})`, :math:`\boldsymbol{J}(\boldsymbol{x}_k)` is the Jacobian of the system of equations and :math:`\boldsymbol{f_J}(\boldsymbol{x})` are the residual equations.
