# artery.fe

Implementation of the 1D blood flow equations in FEniCS.

## Documentation

The documentation was build using Sphinx autobuild and is hosted on Readthedocs https://bloodflow.readthedocs.io

## Installation and dependencies

We recommend installing artery.fe using the provided Dockerfile. This ensures all dependencies are correctly installed. Build the Docker image by running

`docker build --no-cache -t arteryfe:2017.2.0 .`

To create and enter a Docker container run

`docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared "artery.fe:2017.2.0"`

Alternatively, artery.fe can be installed using the provided ``setup.py`` file by running

`python3 setup.py install`

This requires FEniCS version 2017.2.0 or higher.

## Usage

The demo can be run using

`python3 demo_arterybranch.py config/demo_arterybranch.cfg`

## Support

Active contributions and bug reports from users are welcome and encouraged. Should you experience any issues, please do not hesitate to contact @akdiem for advice.

## License

artery.fe is free software made available under the BSD 3-clause License. For details see the LICENSE file.
