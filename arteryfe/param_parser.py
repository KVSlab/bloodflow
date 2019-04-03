__author__ = "Alexandra Diem <alexandra@simula.no>"

from configparser import SafeConfigParser
import argparse
import numpy as np


class ParamParser(object):

    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--cfg")
        args = parser.parse_args()
        self.fparams = args.cfg

        try:
            f = open(self.fparams, 'r')
            data = f.read()
            f.close()
        except Exception as e:
            print(e)
            import sys; sys.exit(1)
        self.param, self.geo, self.solution = self.get_params()


    def get_params(self):
        """
        Reads config file provided in self.fparams.
        Parameters are stored as dictionary params['name'] = value
        """
        config = SafeConfigParser()
        config.optionxform = str
        config.read(self.fparams)

        # Read parameters
        param = ParamParser.get_param_section(config)

        # Read geometry
        geom = ParamParser.get_section(config, "Geometry")

        # Read solution parameters
        sol = ParamParser.get_section(config, "Solution")

        return param, geom, sol


    @staticmethod
    def get_param_section(config):
        """
        Get config file options from section containing strings.

        :param config: ConfigParser object.
        :param section: Name of the section to be read.
        """
        section = 'Parameters'
        options = config.items(section)
        section_dict = {}
        for key, value in options:
            if "," in value:
                value = np.array([float(val) for val in value.split(',')])
            else:
                value = eval(value)
            section_dict[key] = value
        return section_dict


    @staticmethod
    def get_section(config, section):
        """
        Get config file options from section containing strings.

        :param config: ConfigParser object.
        """
        options = config.items(section)
        section_dict = {}
        for key, value in options:
            try:
                section_dict[key] = eval(value)
            except:
                section_dict[key] = value
        return section_dict
