# artery.fe

Implementation of the 1D blood flow equations in FEniCS.

## Documentation

The documentation was build using Sphinx autobuild and is hosted on Readthedocs https://bloodflow.readthedocs.io

## Installation and dependencies

We recommend installing artery.fe using the provided Dockerfile. This ensures all dependencies are correctly installed. Build the Docker image by running

`docker build .`

To create and enter a Docker container run

`docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared "artery.fe:2017.2.0"`

Alternatively, artery.fe can be installed using the provided ``setup.py`` file by running

`python setup.py install`

This requires FEniCS version 2017.2.0 or higher.

<!---
## Attribution

Will be added after JOSS publication
-->

## License

artery.fe is free software made available under the BSD 3-clause License. For details see the LICENSE file.

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
