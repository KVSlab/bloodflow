import sys
import numpy as np
import matplotlib.pyplot as plt

from fenics import *
import configparser

sys.path.insert(0, '../src')
from artery_network import Artery

config = configparser.ConfigParser()
config.read('config.cfg')
par = config['Parameter']
sol = config['Solution']

a = Artery(config.getfloat('Parameter', 'Ru'),
		   config.getfloat('Parameter', 'Rd'),
		   config.getfloat('Parameter', 'L'),
		   config.getfloat('Parameter', 'k1'),
		   config.getfloat('Parameter', 'k2'),
		   config.getfloat('Parameter', 'k3'),
		   config.getfloat('Parameter', 'nu'),
		   config.getfloat('Parameter', 'p0'),
		   config.getfloat('Parameter', 'R1'),
		   config.getfloat('Parameter', 'R2'),
		   config.getfloat('Parameter', 'CT'))

a.solve(config.getint('Solution', 'Nx'),
		config.getint('Solution', 'Nt'),
		config.getfloat('Solution', 'T'),
		config.getint('Solution', 'N_cycles'),
		np.zeros(config.getint('Solution', 'Nt')))
		
plot(a.solution[0].split()[1])
plt.savefig()
