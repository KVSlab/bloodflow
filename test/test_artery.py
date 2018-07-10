import sys
import numpy as np
from scipy.interpolate import interp1d

from fenics import *
import configparser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

Nt = config.getint('Solution', 'Nt')

# Import the inlet flow data
data_q = np.genfromtxt('../data/example_inlet.csv', delimiter = ',')

tt = data_q[:, 0]
qq = data_q[:, 1]

T = data_q[-1, 0]	

q = interp1d(tt, qq)
t = np.linspace(0, T, Nt)
q_ins = q(t)

a.solve(config.getint('Solution', 'Nx'), Nt, T,
		config.getint('Solution', 'N_cycles'), q_ins)

area = np.zeros([a.Nx+1, a.Nt])
flow = np.zeros([a.Nx+1, a.Nt])
pressure = np.zeros([a.Nx+1, a.Nt])

t += (a.N_cycles-1)*a.T

f = interpolate(a.f, a.V)
A0 = interpolate(a.A0, a.V)

for n in range(Nt):
	area[:,n] = a.solution[n].vector().array()[::2]
	flow[:,n] = a.solution[n].vector().array()[1::2]
	pressure[:, n] = 76/101325*(a.p0 + f.vector().array()\
					*(1-np.sqrt(A0.vector().array()/area[:, n])))

x = np.linspace(0, a.L, a.Nx+1)

ttt, xxx = np.meshgrid(t, x)

# Area plot
fig = plt.figure(figsize=(8, 6))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(ttt, xxx, area, rstride=1, cstride=1,  cmap='viridis', linewidth=0, antialiased=False)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('A')
ax.set_ylim(min(x), max(x))
ax.set_xlim(min(t), max(t))
plt.savefig('../output/r0/area.png')


# Flow plot
fig = plt.figure(figsize=(8, 6))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(ttt, xxx, flow, rstride=1, cstride=1,  cmap='viridis', linewidth=0, antialiased=False)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('A')
ax.set_ylim(min(x), max(x))
ax.set_xlim(min(t), max(t))
plt.savefig('../output/r0/flow.png')


# Pressure plot
fig = plt.figure(figsize=(8, 6))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(ttt, xxx, pressure, rstride=1, cstride=1,  cmap='viridis', linewidth=0, antialiased=False)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('A')
ax.set_ylim(min(x), max(x))
ax.set_xlim(min(t), max(t))
plt.savefig('../output/r0/pressure.png')

