import numpy as np

from fenics import *
from mshr import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from configparser import ConfigParser


def unit_to_mmHg(p):
    """Convert pressure units.
    :param p: Pressure in g cm-1 s-2
    :return: Pressure in mmHg
    """
    return 76/101325*p


def mmHg_to_unit(p):
    """Convert pressure units.
    :param p: Pressure in mmHg
    :return: Pressure in g cm-1 s-2
    """
    return 101325/76*p


def adimensionalise_parameters(rc, qc, Ru, Rd, L, k1, k2, k3, rho,
                               nu, p0, R1, R2, CT, q_ins, T):
    """Make quantities independent of units.
    :param rc: Characteristic length (radius)
    :param qc: Characteristic flow
    :param Ru: Radius
    :param Rd: Radius
    :param L: Length
    :param k1: Pressure
    :param k2: Inverse length
    :param k3: Pressure
    :param rho: Density
    :param nu: Viscosity
    :param p0: Pressure
    :param q_ins: Flow
    :param T: Time
    :param R1: Windkessel model resistance
    :param R2: Windkessel model resistance
    :param CT: Windkessel model compliance
    :return: Adimensionalised quantities, including Reynolds' number
    """
    Ru = Ru/rc
    Rd = Rd/rc
    L = L/rc
    k1 = k1*rc**4/rho/qc**2
    k2 = k2*rc
    k3 = k3*rc**4/rho/qc**2
    Re = qc/nu/rc
    nu = nu*rc/qc
    p0 = p0*rc**4/rho/qc**2
    R1 = R1*rc**4/rho/qc
    R2 = R2*rc**4/rho/qc
    CT = CT*rho*qc**2/rc**7
    q_ins = q_ins/qc
    T = T*qc/rc**3
    return Ru, Rd, L, k1, k2, k3, Re, nu, p0, R1, R2, CT, q_ins, T


def adimensionalise(rc, qc, rho, x, nature):
    """Adimensionalise a quantity.
    :param rc: Characteristic radius
    :param qc: Characteristic flow
    :param rho: Density
    :param x: Quantity to be adimensionalised
    :param string nature: Nature of the quantity to be adimensionalised
    :return: adimensionalised quantity
    """
    if nature == 'time':
        x = x*qc/rc**3
    elif nature == 'area':
        x = x/rc**2
    elif nature == 'flow':
        x = x/qc
    elif nature == 'pressure':
        x = x*rc**4/rho/qc**2
    return x


def redimensionalise(rc, qc, rho, x, nature):
    """Give a quantity its units back.
    :param rc: Characteristic radius
    :param qc: Characteristic flow
    :param rho: Density
    :param x: Quantity to be redimensionalised
    :param string nature: Nature of the quantity to be redimensionalised
    :return: Redimensionalised quantity
    """
    if nature == 'time':
        x = x*rc**3/qc
    elif nature == 'area':
        x = x*rc**2
    elif nature == 'flow':
        x = x*qc
    elif nature == 'pressure':
        x = x*rho*qc**2/rc**4
    return x


def read_output(filename):
    """Read data file.
    The cfg-file should have the format in artery_network.dump_metadata.
    :param filename: Location of file to be read
    :return: All variables stored in the file
    """
    config = ConfigParser()
    config.read(filename)

    order = config.getint('data', 'order')
    Nx = config.getint('data', 'Nx')
    Nt = config.getint('data', 'Nt')
    T0 = config.getfloat('data', 'T0')
    T = config.getfloat('data', 'T')
    L = [float(f) for f in config.get('data', 'L').split(',')]
    rc = config.getfloat('data', 'rc')
    qc = config.getfloat('data', 'qc')
    rho = config.getfloat('data', 'rho')
    mesh_locations = config.get('data', 'mesh_locations').split(',')
    names = config.get('data', 'names').split(',')
    locations = config.get('data', 'locations').split(',')

    return order, Nx, Nt, T0, T, L, rc, qc, rho, mesh_locations,\
           names, locations


def XDMF_to_matrix(Nx, Nt, mesh_location, location, name):
    """Read XDMF-file and store it in a matrix.
    :param Nx: Number of spatial steps (number of nodes: Nx+1)
    :param Nt: Number of time steps
    :param mesh_location: Location of mesh to be loaded
    :param location: Location of XDMF-file
    :param name: Name of function within XDMF-file
    :return: Matrix containing values from XDMFFile
    """
    print('Loading %s to matrix.' % (location))
    M = np.zeros([Nx+1, Nt])
    f = XDMFFile(location)
    mesh = Mesh(mesh_location)
    V = FunctionSpace(mesh, 'CG', 1)
    u = Function(V)
    for n in range(Nt):
        f.read_checkpoint(u, name, n)
        M[:, n] = u.vector().get_local()[::-1]
    f.close()
    return M


def plot_matrix(t, x, M, label, output):
    """Create plot of a matrix.
    Store plot in a given location.
    :param t: Vector containing time values
    :param x: Vector containing space values
    :param M: Matrix representing function to be plotted, with dimension t*x
    :param string label: Name of function
    :param string output: Location of output file
    """
    print('Making plot of %s-matrix.' % (label))
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
    print('Saving matrix to %s.' % (output))
    plt.savefig(output)


def is_near(a, b, tol=1.e-14, reltol=1.e-10):
    """Check equality between two floats to a certain tolerance.
    Name contains 'is_' to differentiate it from FEniCS near-function.
    :param a: First number
    :param b: Second number
    :param tol: Tolerance for equality
    :param reltol: Relative tolerance for equality
    :return: True if the two numbers are to be considered equal
    """
    # Neglect relative error if numbers are close to zero
    if np.abs(b) > 1.e-10:
        return np.abs(a-b) < tol or np.abs(a/b-1) < reltol
    else:
        return np.abs(a-b) < tol


def write_file(f, u, label, t):
    set_log_level(40)
    f.write_checkpoint(u, label, t)
    set_log_level(30)


def read_file(f, u, label, i):
    set_log_level(40)
    f.read_checkpoint(u, label, i)
    set_log_level(30)
    return u


def print_progress(n_cycle, n, dt):
    print('Current cycle: %i, Cycle iteration: %i, Time-step %i'\
          % (n_cycle, n, dt), end='\r')
