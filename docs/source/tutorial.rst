Tutorial
=========

To use the package, the Artery_Network file has to be imported. All interaction with the solver goes throught the Artery_Network class. The utils file helps handling data.

Parameters should be without dimension before the package takes them into use. The utils-file provides adimensionalisation methods. For the package to work correctly, an Artery_Network object should be created. Define_geometry should be called next, with spatial and temporal discretisation, and then Define_solution may be called. Solve should be called lastly. This will generate an output folder, containing a file called data.cfg, mesh-files, and folders for area, flow or pressure containing the solution in xdmf-format, according to the specified storage options. All files are enumerated from 0 to the number of arteries in the same way as in the package.

A post processing file is preconfigured to make plots of the data. The only parameter needed is the location of the data.cfg file in the output folder.

A run_from_config file is preconfigured to read parameters from a cfg-file and run the necessary functions in the right order. The structure of the config files may be found in the example-config-files in the config folder. In a FEniCS-enabled terminal window, an example command is:

> python3 run_from_config.py 'config/4cycles.cfg'

The above command will create an output folder containing the solution on four cardiac cycles.

Unit test are provided in the test folder, along with associated configuration files. To run unit tests, one can either run a test file directly, passing the config-file-location as a (string) parameter, or import the file to run the tests individually.
