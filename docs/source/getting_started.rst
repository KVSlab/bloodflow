Getting started
================

artery.fe is a Python package that simulates blood flow through arterial trees in 1D using the finite element (FE) method implemented in FEniCS_.

Installation and dependencies
-----------------------------

Clone the repository::

  $ git clone https://github.com/KVSlab/bloodflow.git

We recommend installing artery.fe using the provided Dockerfile. This
ensures all dependencies are correctly installed. Build the Docker image
by running::

  $ docker build --no-cache -t arteryfe:2017.2.0 .

To create and enter a Docker container run::

  $ docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared "arteryfe:2017.2.0"

Alternatively, artery.fe can be installed using the provided `setup.py` file by running::

  $ python3 setup.py install

artery.fe requires FEniCS_ version 2017.2.0 or higher and Python 3.5 or higher

.. _FEniCS: https://fenicsproject.org/


Running artery.fe
-----------------------------

The file ``run_from_config.py`` provides an example for running a simulation using artery.fe and reproduces the results presented in reference [Kolachalama:2007]::

  $ python3 demo_arterybranch.py config/demo_arterybranch.cfg

to run a simulation over four cardiac cycles, storing the output for the final cardiac cycle only. This automatically creates a directory inside the output directory, which has the same name as the .cfg file used in the simulation. All output is stored in this directory. Additionally, the simulation creates a file 'data.cfg' inside the output directory, which can be used to configure postprocessing and create figures from the output. To produce figures from the output use::

  $ python3 postprocess.py output/4cycles_last/data.cfg

This creates three directories area, flow, and pressure inside the output directory, which contain the corresponding figures.


License
-------

artery.fe is free software made available under the BSD 3-clause
License. For details see the LICENSE file.
