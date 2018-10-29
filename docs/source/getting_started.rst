Getting started
=========

Implementation of the 1D blood flow equations in FEniCS.

Installation and dependencies
-----------------------------

We recommend installing artery.fe using the provided Dockerfile. This
ensures all dependencies are correctly installed. Build the Docker image
by running

``docker build --no-cache -t arteryfe:2017.2.0 .``

To create and enter a Docker container run

``docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared "arteryfe:2017.2.0"``

Alternatively, artery.fe can be installed using the provided
``setup.py`` file by running

``python setup.py install``

This requires FEniCS version 2017.2.0 or higher.


Running artery.fe
-----------------------------

The file ``run_from_config.py`` provides an example for running a simulation using artery.fe and reproduces the results presented in reference [Kolachalama:2007]. Use

``python run_from_config.py config/4cycles_last.cfg``

to run a simulation over four cardiac cycles, storing the output for the final cardiac cycle only. This automatically creates a directory inside the output directory, which has the same name as the .cfg file used in the simulation. All output is stored in this directory. Additionally, the simulation creates a file 'data.cfg' inside the output directory, which can be used to configure postprocessing and create figures from the output. To produce figures from the output use

``python postprocess.py output/4cycles_last/data.cfg``

This creates three directories area, flow, and pressure inside the output directory, which contain the corresponding figures.


License
-------

artery.fe is free software made available under the BSD 3-clause
License. For details see the LICENSE file.

.. raw:: html

   <!--
   ## Other

   To use the package, the Artery_Network file has to be imported. All interaction with the solver goes throught the Artery_Network class. The utils file helps handling data.

   Parameters should be without dimension before the package takes them into use. The utils-file provides adimensionalisation methods. For the package to work correctly, an Artery_Network object should be created. Define_geometry should be called next, with spatial and temporal discretisation, and then Define_solution may be called. Solve should be called lastly. This will generate an output folder, containing a file called data.cfg, mesh-files, and folders for area, flow or pressure containing the solution in xdmf-format, according to the specified storage options. All files are enumerated from 0 to the number of arteries in the same way as in the package.

   A post processing file is preconfigured to make plots of the data. The only parameter needed is the location of the data.cfg file in the output folder.

   A run_from_config file is preconfigured to read parameters from a cfg-file and run the necessary functions in the right order. The structure of the config files may be found in the example-config-files in the config folder. In a FEniCS-enabled terminal window, an example command is:

   > python3 run_from_config.py 'config/4cycles.cfg'

   The above command will create an output folder containing the solution on four cardiac cycles.

   Unit test are provided in the test folder, along with associated configuration files. To run unit tests, one can either run a test file directly, passing the config-file-location as a (string) parameter, or import the file to run the tests individually.
   -->
