"""
Implementation of the 1D blood flow equations.
Author : Syver Døving Agdestein
"""

from fenics import *
import numpy as np

T = 1.0
L = 20.0
Nt = 100
Nx = 100

r0 = 1.0
E = 1.0
H = 0.01
Re = 1.0
Deltab = 0.1

f = 4*E*H/3/r0

mesh = Rectangle(Point(0,0), Point(T,L), Nt Nx)
V = FunctionSpace(mesh, 'P', degree = 1)

A = Function(V)
q = Function(V)



U = (A,q)

F = (q, pow(q,2)/A + sqrt(A0*A))

S = (0, 2*pi*R/deltab/Re*q/A
