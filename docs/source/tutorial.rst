Tutorial
=========

This tutorial provides a description of the demo files provided with artery.fe.

Running a simulation
--------------------

The base folder of the artery.fe repository contains the demo file `run_from_config.py`, which we will disect here to demonstrate how to run a simulation. Start by importing artery.fe and applying the short name 'af' for convenience::

  import arteryfe as af

Parameters are loaded from a .cfg file using the class `ParamParser`, where the .cfg file is stored in `config_location`::

  param = af.ParamParser(config_location)

Parameters can then be explicitly loaded using the syntax::

  order = param.param['order']

Inlet flow rates for the parent artery should be prescribed using a .csv file whose location is provided in the parameter `inlet_flow_location`. It should not contain any header lines. The first column should correspond to time points, while the second column should correspond to the flow rate values.

Parameters are nondimensionalised and the Reynold's number is calculated using::

  Ru, Rd, L, k1, k2, k3, Re, nu, p0, R1, R2, CT, q_ins, T =\
        af.nondimensionalise_parameters(rc, qc, Ru, Rd, L, k1, k2, k3,
                                   rho, nu, p0, R1, R2, CT, q_ins, T)

We are now ready to create an `ArteryNetwork` object. This is done in three steps, where the first step calls the `ArteryNetwork` constructor, the second step calls a function to set up the geometry of the artery, and the third step calls a function to set up the solver for the problem::

  an = af.ArteryNetwork(order, rc, qc, Ru, Rd, L, k1, k2,
                      k3,	rho, Re, nu, p0, R1, R2, CT)
  an.define_geometry(Nx, Nt, T, N_cycles)
  an.define_solution(output_location, q_ins[0], theta)

The solution is then calculated using::

  an.solve(q_ins, Nt_store, N_cycles_store, store_area, store_pressure)

The solution is stored in `output` in a separate folder defined in the parameter `output_location`.


Visualising the output
----------------------

The base folder of artery.fe additionally includes a file `postprocess.py`, which handles the postprocessing and visualisation of results. As before we import artery.fe, and here also numpy::

  import arteryfe as af
  import numpy as np

To read parameters back from the output use::

  order, Nx, Nt, T0, T, L, rc, qc, rho, mesh_locations, names, locations =
        read_output(data_location)

In the example case provided in the demo the simulation runs for four cardiac cycles, but we are only interested in plotting the solution for the last cardiac cycle. We use a Numpy array to define the time variable on the redimensionalised time parameters::

  T0 = redimensionalise(rc, qc, rho, T0, 'time')
  T = redimensionalise(rc, qc, rho, T, 'time')
  t = np.linspace(T0, T, Nt)

The variable `names` contains the names of the unknowns computed during the artery.fe simulation, which are flow rate ('flow'), cross-sectional area ('area') and arterial pressure ('pressure'). Thus, the first for loop::

  for i, name in enumerate(names):

creates the same visualisation for each unknown, while the second for loop::

  for j in range(2**order-1):

iterates over all arteries in the geometry. FEniCS writes data to XDMF and HDF5 files, and these can be converted to Numpy matrices using::

  M = XDMF_to_matrix(Nx, Nt, mesh_locations[j],
            '%s/%s_%i.xdmf' % (locations[i], name, j), name)
  M = redimensionalise(rc, qc, rho, M, name)

Redimensionalised pressure is given in units of Pa, but for clinical purposes the use of mmHg (millimetres of mercury) is much more common::

  if name == 'pressure':
    M = unit_to_mmHg(M)

The spatial variable for the plots for each artery is given by::

  x = np.linspace(0, L[j], Nx+1)

And lastly, the plots are created using::

  plot_matrix(t, x, M, name, '%s/%s_%i.png' % (locations[i], name, j))
