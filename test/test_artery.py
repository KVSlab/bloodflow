import sys
import numpy as np
from scipy.interpolate import interp1d

from fenics import *
import configparser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sys.path.insert(0, '../src')
from artery_network import Artery
import utils


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

Nt = config.getint('Solution', 'Nt')

# Import the inlet flow data
data_q = np.genfromtxt('../data/example_inlet.csv', delimiter = ',')

tt = data_q[:, 0]
qq = data_q[:, 1]

T = data_q[-1, 0]	

q = interp1d(tt, qq)
t = np.linspace(0, T, Nt)
q_ins = q(t)
#q_ins = np.zeros(Nt)

a.solve(config.getint('Solution', 'Nx'), Nt, T,
		config.getint('Solution', 'N_cycles'), q_ins)

area = np.zeros([a.Nx+1, a.Nt])
flow = np.zeros([a.Nx+1, a.Nt])
pressure = np.zeros([a.Nx+1, a.Nt])

t += (a.N_cycles-1)*a.T

f = interpolate(a.f, a.V).vector().array()[::-1]
A0 = interpolate(a.A0, a.V).vector().array()[::-1]


for n in range(Nt):
	area[:,n] = (a.solution[n].vector().array()[::2])[::-1]
	flow[:,n] = (a.solution[n].vector().array()[1::2])[::-1]
	pressure[:, n] = utils.unit_to_mmHg(a.p0 + (f*(1-np.sqrt(A0)/area[:, n])))

x = np.linspace(0, a.L, a.Nx+1)

# Area plot
utils.plot_matrix(t, x, area, 'area', '../output/r0/area.png')

# Flow plot
utils.plot_matrix(t, x, flow, 'flow', '../output/r0/flow.png')

# Pressure plot
utils.plot_matrix(t, x, pressure, 'pressure', '../output/r0/pressure.png')

