import numpy as np
from fenics import *

def solve(mesh, V, V2, Nx, Nt, L, T, q_ins, A0, U, F, S, compute_A_out):
	""" Return the solution of dU/dt + dF/dx = S.
	:param mesh: 1D mesh of the interval [0, L] with Nx points
	:param V: Scalar function space on mesh
	:param V2: Vector function space on mesh
	:param Nx: Number of spatial steps
	:param Nt: Number of temporal steps
	:param q_in: Inlet boundary condition for the second component of U
	:param A0: Initial value of the first component of U
	:param F: Flux function
	:param S: Source function
	:param compute_A_out: Return next outlet boundary value of A
	"""
	# Define trial functions
	A, q = split(U)

	# Define test functions
	v1, v2 = TestFunctions(V2)
	v = (v1, v2)

	# Inlet flow at a given time t_n (initially t_0)
	q_in = Function(V)
	q_in.assign(Constant(q_ins[0]))

	# Outlet area
	A_out = Function(V)
	A_out.assign(Constant(A0(L)))

	# Initial value
	U_n = Function(V2)
	U_n.assign((A0, q_in))

	# Spatial boundary conditions
	tol = 1.e-14
	def inlet_bdry(x, on_boundary):
		return on_boundary and near(x[0], 0, tol)	
	def outlet_bdry(x, on_boundary):
		return on_boundary and near(x[0], L, tol)
	bc_outlet = DirichletBC(V2.sub(0), A_out, outlet_bdry)
	bc_inlet = DirichletBC(V2.sub(1), q_in, inlet_bdry)
	bcs = [bc_inlet, bc_outlet]


	# Variational form
	FF = dot(U, v)*dx\
	   + dt*dot(F, v)*ds\
	   - dt*dot(F, grad(v))*dx\
	   + dt*dot(S, v)*dx\
	   - dot(U_n, v)*dx

	# Progress bar
	progress = Progress('Time-stepping')
	set_log_level(PROGRESS)

	t = 0
	
	# Time-stepping
	for n in range(0, Nt-1):
		
		print('Iteration '+str(n))
		
		t += dt
		
		# Update inlet boundary condition
		q_in.assign(Constant(qq[n+1]))
		
		# Update outlet boundary condition
		A_out_value = compute_A_out(U_n)
		A_out.assign(Constant(A_out_value))

		# U is solution of FF == 0
		solve(FF == 0, U, bcs)
		
		# Update previous solution
		U_n.assign(U)
		
		# Update progress bar
		progress.update(t)
