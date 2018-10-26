import numpy as np

from dolfin import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from configparser import SafeConfigParser
from scipy.interpolate import interp1d


def unit_to_mmHg(p):
    """
    Converts pressure value in g cm-1 s-2 to mmHg.

    Arguments
    ---------
    p : float
        Pressure value in g cm-1 s-2

    Returns
    -------
    return : float
        Pressure value in mmHg
    """
    return 76/101325*p


def mmHg_to_unit(p):
    """
    Converts pressure value in mmHg to g cm-1 s-2.

    Arguments
    ---------
    p : float
        Pressure value in mmHg

    Returns
    -------
    return : float
        Pressure value in g cm-1 s-2
    """
    return 101325/76*p


def nondimensionalise_parameters(rc, qc, Ru, Rd, L, k1, k2, k3, rho,
       nu, p0, R1, R2, CT, q_ins, T):
    """
    Nondimensionalise parameters.

    Arguments
    ---------
    rc : float
        Characteristic radius (length)
    qc : float
        Characteristic flow
    Ru : float
        Upstream radius
    Rd : float
        Downstream radius
    L : float
        Vessel length
    k1 : float
     	First constant from the relation Eh/r0
    k2 : float
        Second constant from the relation Eh/r0
    k3 : float
        Third constant from the relation Eh/R0
    rho : float
        Density of blood
    Re : float
        Reynolds' number
    nu : float
        Viscosity of blood
    p0 : float
        Diastolic pressure

    Returns
    -------
    return : tuple
        Tuple of dimensionless quantities, including Reynold's number
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


def nondimensionalise(rc, qc, rho, x, nature):
    """
    Nondimensionalise a parameter using the characteristic parameters.

    Arguments
    ---------
    rc : float
        Characteristic radius (length)
    qc : float
        Characteristic flow
    rho : float
        Density of blood
    x : float
        Parameter to redimensionalise
    nature : string
        Nature of parameter to be redimensionalised

    Returns
    -------
    return : float
     	Dimensionless quantity
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
    """
    Redimensionalise a parameter using the characteristic parameters.

    Arguments
    ---------
    rc : float
        Characteristic radius (length)
    qc : float
        Characteristic flow
    rho : float
        Density of blood
    x : float
        Parameter to redimensionalise
    nature : string
        Nature of parameter to be redimensionalised

    Returns
    -------
    return : float
        Redimensionalised quantity
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


def read_inlet(data_location, Nt):
    """
    Read inlet flow rate data from file and returns

    Arguments
    ---------
    data_location: string
        Location of inlet flow data file
    Nt : int
        Number of time steps

    Returns
    -------
    return : tuple
        Length of a cardiac cycle, inlet flow rate data
    """
    # Import inlet flow data
    data_q = np.genfromtxt(data_location, delimiter=',')
    tt = data_q[:, 0]
    qq = data_q[:, 1]
    T = data_q[-1, 0]

    # Interpolate inlet flow
    q = interp1d(tt, qq)
    t = np.linspace(0, T, Nt)
    q_ins = q(t)

    return T, q_ins


def read_output(filename):
    """
    Read data file generated in the output folder.

    Arguments
    ---------
    rc : string
        Data file name

    Returns
    -------
    return : tuple
        Tuple of all parameters stored in the file
    """
    config = SafeConfigParser()
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
    """
    Read XDMF-file and store it in a Numpy array.

    Arguments
    ---------
    Nx : int
        Number of spatial points per artery
    Nt : int
        Number of time steps per cardiac cycle
    mesh_location : string
        File name of the mesh to be loaded
    location : string
        File name of the XDMF file
    name : string
        Name of the function in the XDMF file

    Returns
    -------
    return : numpy.array
        Array containing the values loaded from the XDMF file
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
    """
    Creates and stores a plot of a numpy array.

    Arguments
    ---------
    t : numpy.array
        Array containing time point values
    x : numpy.array
        Array containing spatial point values
    M : numpy.array
        Array to be plotted, with dimension (t, x)
    output : string
        File name to store plot
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


def is_near(a, b, tol=1.e-11, reltol=1.e-10):
    """
    Check near-equality between two floats to a certain tolerance. Name
    contains 'is' to differentiate it from DOLFIN near()-function.

    Arguments
    ---------
    a : float
        First number
    b : float
        Second number
    tol : float
        Tolerance for near-equality
    reltol : float
        Relative tolerance for near-equality

    Returns
    -------
    return : boolean
        True if a and b are near-equal
    """
    # Neglect relative error if numbers are close to zero
    if np.abs(b) > 1.e-10:
        return np.abs(a-b) < tol or np.abs(a/b-1) < reltol
    else:
        return np.abs(a-b) < tol


def write_file(f, u, label, t):
    """
    Wrapper function for DOLFIN's write_checkpoint() to avoid print of warnings.

    Arguments
    ---------
    f : XDMFFile
        XDMFFile to write to
    u : Function
        Function to write to f
    label : string
        Label for the Function
    t : float
        Time point
    """
    set_log_level(40)
    f.write_checkpoint(u, label, t)
    set_log_level(30)


def read_file(f, u, label, i):
    """
    Wrapper function for DOLFIN's read_checkpoint() to avoid print of warnings.

    Arguments
    ---------
    f : XDMFFile
        XDMFFile to write to
    u : Function
        Function to write to f
    label : string
        Label for the Function
    i : int
        Index of time point
    """
    set_log_level(40)
    f.read_checkpoint(u, label, i)
    set_log_level(30)


def print_progress(n_cycle, n, dt):
    """
    Print progress of simulation.

    Arguments
    ---------
    n_cycle : int
        Current cardiac cycle
    n : int
        Current iteration within cardiac cycle
    dt : float
        Current time step
    """
    print('Current cycle: %i, Cycle iteration: %i, Time-step %i'\
    % (n_cycle, n, dt), end='\r')
