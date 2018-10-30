.. _theory:

Blood flow dynamics in 1D
=========================

artery.fe implements the 1D system of equations derived by [Olufsen:2000] that is commonly used to numerically model blood flow dynamics. Its derivation is provided here for the mathematically interested reader and for completeness.

Arteries are modelled as axisymmetric tubes in a cylindrical coordinate system with radial direction *r* and axial direction *z*. Therefore, the continuity equation is

.. math::

  \frac{\partial u_z(r,z,t)}{\partial z} + \frac{1}{r} \frac{\partial(ru_r(r,z,t))}{\partial r} = 0 \qquad (1),

where :math:`\boldsymbol{u} = (u_z(r,z,t), u_r(r,z,t))`. Integration of (1) over the cross-sectional area with :math:`A(z,t) = \pi R(z,t)^2` yields

.. math::
  2 \pi \int_0^{R(z,t)} \left( \frac{\partial u_z(r,z,t)}{\partial z} + \frac{1}{r} \frac{\partial(ru_r(r,z,t))}{\partial r} \right) r dr = 0,\\
  2 \pi \int_0^{R(z,t)} \frac{\partial u_z(r,z,t)}{\partial z} r dr + 2 \pi \left[ r u_r(r,z,t) \right]_0^{R(z,t)} = 0 \qquad (2)

where :math:`R(z,t)` describes the vessel radius. The first term of (2) can be evaluated using Leibniz' integral rule ("differentiation under the integral sign")

.. math::

  \frac{\partial}{\partial x} \int_{a(x)}^{b(x)} f(x, y) dy = \int_{a(x)}^{b(x)} \frac{\partial f(x,y)}{\partial x} dy + f(x,b(x)) \frac{\partial b(x)}{\partial x} - f(x,a(x)) \frac{\partial a(x)}{\partial x}

and results in

.. math::
  2 \pi \frac{\partial}{\partial z} \int_0^{R(z,t)} u_z(r,z,t) r dr - 2 \pi \frac{\partial R(z,t)}{\partial z}\left[ r u_z(r,z,t) \right]_{R(z,t)} + 2 \pi \left[ r u_r(r,z,t) \right]_0^{R(z,t)} = 0.

Defining flux through the vessel as

.. math::
  q(r,z,t) = 2\pi \int_0^{R(z,t)} u_z(r,z,t) r dr


yields

.. math::
  \frac{\partial q(r,z,t)}{\partial z} - 2 \pi \frac{\partial R(z,t)}{\partial z}\left[ r u_z(r,z,t) \right]_{R(z,t)} + 2 \pi \left[ r u_r(r,z,t) \right]_0^{R(z,t)} = 0. \qquad (3)

Due to no-slip :math:`r` vanishes at the boundary :math:`R(z,t)`

.. math::
  \left. u_z(r,z,t) \right|_{R(z,t)} = 0 \qquad (4)

and longitudinal tether of the vessel leads to

.. math::

  \left. u_r(r,z,t) \right|_{R(z,t)} = \frac{\partial R(z,t)}{\partial t}. \qquad(5)

Thus, (3) reduces to

.. math::

  \frac{\partial q(r,z,t)}{\partial z} + \frac{\partial A(z,t)}{\partial t} = 0. \qquad (6)

The momentum equation for axisymmetric flow with no swirl is

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

To solve the remaining terms it is necessary to make assumptions about the velocity profile of blood flow through an artery. Blood flow is considered pulsatile laminar and vessels are considered slightly tapered, therefore the velocity profile is assumed to be mostly flat with a thin boundary layer with cardiac cycle length :math:`T` and width :math:`\delta_b = (\nu T / (2\pi))^{0.5}`, such that :math:`\delta_b \ll R(z,t)`. The axial velocity :math:`u_z(r,z,t)` thus has the form

.. math::

  u_z(r,z,t) = \begin{cases}
  \bar{u}_z(z,t) & r \leq R(z,t)-\delta_b\\
  \bar{u}_z(z,t) (R(z,t)-r)/\delta_b & R(z,t)-\delta_b < r \leq R(z,t),
  \end{cases} \qquad (10)

where :math:`\bar{u}_z(z,t)` is the average axial velocity outside the boundary layer. This leads to a flat velocity profile outside the boundary layer and linearly increasing profile (from 0 to :math:`\bar{u}_z(z,t)`) inside the boundary layer. Note that a physiological cardiac cycle at rest has between 40 and 70 beats per minute (0.6 s :math:`\leq T \leq` 1.1 s), therefore the boundary layer is 0.07--0.09 cm in size. The minimal inlet radius of arteries considered in this work is 0.14 cm, therefore (10) is appropriate for the desired velocity profile. The first and second terms of (9) can then be expressed as a power series in :math:`\delta_b`

.. math::

  q = 2\pi \int_0^{R(z,t)} u_z(r,z,t) r dr = A \bar{u}_z(z,t) \left( 1 - \frac{\delta_b}{R(z,t)} + \mathcal{O}(\delta_b^2) \right),\\
  2\pi \frac{\partial}{\partial z} \int_0^{R(z,t)} u_z(r,z,t)^2 r dr = A \bar{u}_z(z,t) \left( 1 - \frac{4}{3} \frac{\delta_b}{R(z,t)} + \mathcal{O}(\delta_b^2) \right).

Using these solutions the second term of (9) becomes

.. math::

  2\pi \frac{\partial}{\partial z} \int_0^{R(z,t)} u_z(r,z,t)^2 r dr = \frac{q(z,t)^2}{A(z,t)} \left( 1 + \frac{2}{3} \frac{\delta_b}{R(z,t)} + \mathcal{O}(\delta_b^2) \right).

This leaves the term on the right-hand side (RHS) of (9) to be evaluated using the velocity profile

.. math::

  2 \pi \nu R(z,t) \frac{\partial u_z(r,z,t)}{\partial r} = - \frac{2 \pi \nu R(z,t)}{\delta_b} \frac{q(z,t)}{A(z,t)} + \mathcal{O}(\delta_b)

such that finally, keeping only leading order terms in :math:`\delta_b`, the momentum equation reads

.. math::

  \frac{\partial q(z,t)}{\partial t} + \frac{\partial}{\partial z} \left( \frac{q(z,t)^2}{A(z,t)} \right) + \frac{A(z,t)}{\rho} \frac{\partial p(z,t)}{\partial z} = - \frac{2 \pi \nu R(z,t)}{\delta_b} \frac{q(z,t)}{A(z,t)}. \qquad (11)

In order to solve the system of (6) and (11) they need to be written in conservation form

.. math::

  \frac{\partial \boldsymbol{U}}{\partial t} + \frac{\partial \boldsymbol{F}}{\partial z} = \boldsymbol{S}. \qquad (12)

The quantity :math:`B` is introduced and chosen to fulfill

.. math::
  B(r_0(z), p(z,t)) = \frac{1}{\rho} \int A(z,t) dp(z,t),

with :math:`r_0(z)` initial radius at rest such that

.. math:

  \frac{\partial B(r_0(z), p(z,t))}{\partial z} = \frac{A}{\rho} \frac{\partial p(z,t)}{\partial z} + \frac{\partial B(r_0(z), p(z,t))}{\partial r_0(z)} \frac{\partial r_0(z)}{\partial z}

Then, adding the term :math:`(\partial B / \partial r_0) (\partial r_0 / \partial z)` to both sides of (11), the system of equations can be written in conservation form

.. math::

  \begin{split}
  \dfrac{\partial}{\partial t} \begin{pmatrix} A(z,t) \\ q(z,t) \end{pmatrix} + \dfrac{\partial}{\partial z} & \begin{pmatrix} q(z,t)\\ \dfrac{q(z,t)^2}{A(z,t)} + B(r_0(z), p(z,t)) \end{pmatrix} =\\ & \qquad \begin{pmatrix} 0 \\ - \dfrac{2 \pi \nu R(z,t)}{\delta_b} \dfrac{q(z,t)}{A(z,t)} + \dfrac{\partial B(r_0(z), p(z,t))}{\partial r_0(z)} \dfrac{\partial r_0(z)}{\partial z} \end{pmatrix}. \qquad (13)
  \end{split}

Currently, (13) contains three unknowns (:math:`q, A, p`) for two equations, thus a third relation is needed to solve the system of equations. The aforementioned equation, referred to as the state equation, describes the relationship between :math:`A(z,t)` and :math:`p(z,t)`. One choice for the state equation is the linearly elastic relation

.. math:

  p(z,t) - p_0 = \frac{4}{3} \frac{Eh}{r_0(z)} \left( 1 - \sqrt{\frac{A_0(z)}{A(z,t)}} \right) \qquad (14),

where the constant $p_0$ is the diastolic pressure, $E$ is the Young's modulus of the vessel wall, $h$ is the wall width and $A_0(z) = \pi r_0(z)^2$. The relationship $Eh/r_0$ is based on compliance estimates

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
  &\begin{pmatrix} 0 \\ -\dfrac{2 \pi \nu q(z,t) R(z,t)}{\delta_b A(z,t)} + \dfrac{1}{\rho} \left( 2 \sqrt{A(z,t)} \left( \sqrt{\pi} f(r_0) + \sqrt{A_0(z)} \frac{df(r_0)}{dr_0 } \right) - A(z,t) \dfrac{df(r_0)}{dr_0} \right) \dfrac{dr_0(z)}{dz} \end{pmatrix}. \qquad (16)
  \end{split}

To nondimensionalise we define some characteristic parameters

================================  ================
Parameter                         Physical meaning
================================  ================
R                                 radius
Q                                 flow rate
:math:`\rho`                      blood density
:math:`\nu`                       blood viscosity
:math:`\mathcal{Re} = Q/(\nu R)`  Reynold's number
================================  ================

and rescale variables accordingly

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
  &\begin{pmatrix} 0 \\ -\dfrac{2 \pi R(z,t)}{\delta_b \mathcal{Re}} \dfrac{q(z,t)}{A(z,t)} +\left( 2 \sqrt{A(z,t)} \left( \sqrt{\pi} f(r_0) + \sqrt{A_0(z)} \frac{df(r_0)}{dr_0 } \right) - A(z,t) \dfrac{df(r_0)}{dr_0} \right) \dfrac{dr_0(z)}{dz} \end{pmatrix}. \qquad (17)
  \end{split}
