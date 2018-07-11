import sys
import numpy as np
from scipy.interpolate import interp1d

from fenics import *
import configparser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sys.path.insert(0, '../src')

from utils import *
from artery import Artery
import conservation_solver as cs


config = configparser.ConfigParser()
config.read('config.cfg')

a = Artery(config.getfloat('Parameters', 'Ru'),
		   config.getfloat('Parameters', 'Rd'),
		   config.getfloat('Parameters', 'L'),
		   config.getfloat('Parameters', 'k1'),
		   config.getfloat('Parameters', 'k2'),
		   config.getfloat('Parameters', 'k3'),
		   config.getfloat('Parameters', 'nu'),
		   config.getfloat('Parameters', 'p0'),
		   config.getfloat('Parameters', 'R1'),
		   config.getfloat('Parameters', 'R2'),
		   config.getfloat('Parameters', 'CT'))

Nt = config.getint('Geometry', 'Nt')

# Import the inlet flow data
data_q = np.genfromtxt('../data/example_inlet.csv', delimiter = ',')

tt = data_q[:, 0]
qq = data_q[:, 1]

T = data_q[-1, 0]

q = interp1d(tt, qq)
t = np.linspace(0, T, Nt)
q_ins = q(t)
#q_ins = np.zeros(Nt)

a.define_geometry(config.getint('Geometry', 'Nx'), Nt, T,
		config.getint('Geometry', 'N_cycles'))

cs.solve_artery(a, q_ins)

area = np.zeros([a.Nx+1, a.Nt])
flow = np.zeros([a.Nx+1, a.Nt])
pressure = np.zeros([a.Nx+1, a.Nt])

t += (a.N_cycles-1)*a.T


f = interpolate(a.f, a.V).vector().array()[::-1]
A0 = interpolate(a.A0, a.V).vector().array()[::-1]


for n in range(Nt):
	area[:, n] = (a.solution[n].vector().array()[::2])[::-1]
	flow[:, n] = (a.solution[n].vector().array()[1::2])[::-1]
	pressure[:, n] = a.pressure(f, A0, area[:, n])

x = np.linspace(0, a.L, a.Nx+1)

# Area plot
plot_matrix(t, x, area, 'area', '../output/r0/area.png')

# Flow plot
plot_matrix(t, x, flow, 'flow', '../output/r0/flow.png')

# Pressure plot
plot_matrix(t, x, unit_to_mmHg(pressure), 'pressure', '../output/r0/pressure.png')

