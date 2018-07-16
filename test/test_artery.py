import sys
import numpy as np
from scipy.interpolate import interp1d

from fenics import *
import configparser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sys.path.insert(0, 'src/')

from utils import *
from artery_network import Artery_Network

config = configparser.ConfigParser()
config.read('test/config.cfg')

order = config.getint('Parameters', 'order')
Ru = [float(f) for f in config.get('Parameters', 'Ru').split(',')]
Rd = [float(f) for f in config.get('Parameters', 'Rd').split(',')]
L = [float(f) for f in config.get('Parameters', 'L').split(',')]


k1 = config.getfloat('Parameters', 'k1')
k2 = config.getfloat('Parameters', 'k2')
k3 = config.getfloat('Parameters', 'k3')
nu = config.getfloat('Parameters', 'nu')
p0 = config.getfloat('Parameters', 'p0')
R1 = config.getfloat('Parameters', 'R1')
R2 = config.getfloat('Parameters', 'R2')
CT = config.getfloat('Parameters', 'CT')

Nt = config.getint('Geometry', 'Nt')
Nx = config.getint('Geometry', 'Nx')
N_cycles = config.getint('Geometry', 'N_cycles')

# Import the inlet flow data
data_q = np.genfromtxt('data/example_inlet.csv', delimiter = ',')
tt = data_q[:, 0]
qq = data_q[:, 1]
T = data_q[-1, 0]
q = interp1d(tt, qq)
t = np.linspace(0, T, Nt)
q_ins = q(t)
#q_ins = np.zeros(Nt)

# Create artery network
an = Artery_Network(order, Ru, Rd, L, k1, k2, k3, nu, p0, R1, R2, CT)
an.define_geometry(Nx, Nt, T, N_cycles)
an.define_solution(q_ins[0])
an.solve(q_ins)
"""
a = an.arteries[0]

area = np.zeros([a.Nx+1, a.Nt])
flow = np.zeros([a.Nx+1, a.Nt])
pressure = np.zeros([a.Nx+1, a.Nt])

t += (a.N_cycles-1)*a.T


f = interpolate(a.f, a.V).vector().array()[::-1]
A0 = interpolate(a.A0, a.V).vector().array()[::-1]


for n in range(Nt):
	area[:, n] = (a.solution[n].vector().array()[::2])[::-1]  # "U[0]"
	flow[:, n] = (a.solution[n].vector().array()[1::2])[::-1]  # "U[1]"
	pressure[:, n] = a.pressure(f, A0, area[:, n])

x = np.linspace(0, a.L, a.Nx+1)

# Area plot
plot_matrix(t, x, area, 'area', '../output/r0/area.png')

# Flow plot
plot_matrix(t, x, flow, 'flow', '../output/r0/flow.png')

# Pressure plot
plot_matrix(t, x, unit_to_mmHg(pressure), 'pressure', '../output/r0/pressure.png')
"""
