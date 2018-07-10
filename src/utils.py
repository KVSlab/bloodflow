import numpy as np

import matplotlib.pyplot as plt


# Pressure unit converting functions ('unit' = g cm-1 s-2)
def unit_to_mmHg(p):
	return 76/101325*p
	
def mmHg_to_unit(p):
	return 101325/76*p
	


def plot_matrix(t, x, M, label, output):
	"""Create a plot of a matrix.
	Store it in a given location.
	:param t: Vector containing time values
	:param x: Vector containing space values
	:param M: Matrix representing function to be plotted, with dimension = t*x
	:param string label: Name of the function
	:param string output: Location (with filename) of output file
	"""
	
	T, X = np.meshgrid(t, x)
	fig = plt.figure(figsize=(8, 6))
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(T, X, M, rstride=1, cstride=1,  cmap='viridis',
						   linewidth=0, antialiased=False)
	ax.set_xlabel('t')
	ax.set_ylabel('x')
	ax.set_zlabel(label)
	ax.set_ylim(min(x), max(x))
	ax.set_xlim(min(t), max(t))
	plt.savefig(output)
