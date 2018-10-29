Tutorial
=========

The base folder of the artery.fe repository contains the demo file `run_from_config.py`, which we will disect here to demonstrate how to run a simulation. Start by importing artery.fe and applying the short name 'af' for convenience.

``import arteryfe as af``

Parameters are loaded from a .cfg file using the class `ParamParser`, where the .cfg file is stored in `config_location`

``param = af.ParamParser(config_location)``

Parameters can then be explicitly loaded using the syntax

``order = param.param['order']``

Inlet flow rates for the parent artery should be prescribed using a .csv file whose location is provided in the parameter `inlet_flow_location`. It should not contain any header lines. The first column should correspond to time points, while the second column should correspond to the flow rate values.

Parameters are nondimensionalised and the Reynold's number is calculated using

.. code-block:: python
Ru, Rd, L, k1, k2, k3, Re, nu, p0, R1, R2, CT, q_ins, T =\
        af.nondimensionalise_parameters(rc, qc, Ru, Rd, L, k1, k2, k3,
                                   rho, nu, p0, R1, R2, CT, q_ins, T)

We are now ready to create an `ArteryNetwork` object. This is done in three steps, where the first step calls the `ArteryNetwork` constructor, the second step calls a function to set up the geometry of the artery, and the third step calls a function to set up the solver for the problem.

.. code-block:: python
an = af.ArteryNetwork(order, rc, qc, Ru, Rd, L, k1, k2,
                      k3,	rho, Re, nu, p0, R1, R2, CT)
an.define_geometry(Nx, Nt, T, N_cycles)
an.define_solution(output_location, q_ins[0], theta)

The solution is then calculated using

``an.solve(q_ins, Nt_store, N_cycles_store, store_area, store_pressure)``

The solution is stored in `output` in a separate folder defined in the parameter `output_location`.
