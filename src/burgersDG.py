from fenics import *
from mshr import *
import numpy as np

T = 0.5             # final time
num_steps = 250      # number of time steps
dt = T / num_steps   # time step size

Nx, Ny = 100, 100


#mesh = RectangleMesh(Point(0.0,0.0), Point(1.0,0.2), Nx, Ny)
mesh = UnitSquareMesh(Nx,Ny)
n = FacetNormal(mesh)
h = mesh.hmin()
h_avg = avg(h)
V = FunctionSpace(mesh, 'DG', 1)

u_D = Expression('sin(x[0]*pi)*sin(x[1]*pi)', degree=2)
#u_D = Expression('x[0]*(1-x[0])', degree=3)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, u_D, boundary)

u_n = interpolate(u_D, V)

u = Function(V)
v = TestFunction(V)

# Transport direction
b = Constant((1.0,0))

# Penalty weight
alpha = Constant(10)

F = u*v*dx - dt/2*pow(u,2)*dot(b, grad(v))*dx\
  + dt/2*jump(b*pow(u,2), n)*avg(v)*dS\
  + dt/2*dot(jump(v, n), avg(b*pow(u,2)))*dS\
  + dt/2*alpha/h_avg*jump(u)*jump(v)*dS\
  + dt/2*dot(n, b)*pow(u,2)*v*ds\
  + dt/2*alpha/h_avg*u*v*ds\
  - u_n*v*dx


# Create XDMF files for visualization output
xdmffile_u = XDMFFile('10.xdmf')

# Create progress bar
progress = Progress('Time-stepping')
set_log_level(PROGRESS)

t = 0

for n in range(num_steps):

	t += dt
	solve(F == 0, u, bc)
	u_n.assign(u)
	
	# Save solution to file (XDMF/HDF5)
	xdmffile_u.write(u, t)

	progress.update(t / T)
	


