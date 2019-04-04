import sys

import numpy as np
import configparser

import arteryfe as af


def main(config_location):
    """Read config-file.
    Run the necessary functions to compute the solution.
    :param string config_location: Location of config file
    """
    param = af.ParamParser()

    # Create artery network
    an = af.ArteryNetwork(param)

    # Solve problem and store data
    an.solve()


if __name__ == '__main__':
    main(sys.argv[1])
