__author__ = "Alexandra Diem <alexandra@simula.no>"

from configparser import SafeConfigParser
import numpy as np


class ParamParser(object):

    def __init__(self, fparams):
        self.fparams = fparams
        try:
            f = open(fparams, 'r')
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

        # Read geometry tags section
        geo = ParamParser.get_geo_section(config)

        # Read solution tags
        solution = ParamParser.get_solution_section(config)

        return param, geo, solution


    @staticmethod
    def get_param_section(config):
        """
        Get config file options from section.

        :param config: ConfigParser object.
        """
        section = 'Parameters'
        options = config.items(section)
        section_dict = {}
        for key, value in options:
            if "," in value:
                value = np.array([float(v) for v in value.split(',')])
            elif key == 'order':
                value = int(value)
            else:
                value = float(value)
            section_dict[key] = value
        return section_dict


    @staticmethod
    def get_geo_section(config):
        """
        Get config file options from section.

        :param config: ConfigParser object.
        """
        section = 'Geometry'
        options = config.items(section)
        section_dict = {}
        for key, value in options:
            section_dict[key] = int(value)
        return section_dict


    @staticmethod
    def get_solution_section(config):
        """
        Get config file options from section.

        :param config: ConfigParser object.
        """
        section = 'Solution'
        options = config.items(section)
        section_dict = {}
        for key, value in options:
            if 'location' in key:
                section_dict[key] = value
            elif key == 'theta':
                section_dict[key] = float(value)
            else:
                section_dict[key] = int(value)
        return section_dict
