---
title: 'ArtEniCS: Implementing the 1D blood flow equations in FEniCS'
tags:
- Python
- FEniCS
- Blood flow
authors:
- name: Syver Døving Agdestein
  orcid: 0000-0002-1589-2916
- name: A. K. Diem
date: 24 July 2018
bibliography: paper.bib
---

# Summary

This package implements the 1D blood flow equations using FEniCS. The package provides tools for modeling a network of arteries and solving the equations governing blood flow on the arterial network.

FEniCS is a finite element framework that is being increasingly used for cardiovascular modeling. At the present moment, there does not seem to be any implementations of the 1D blood flow equations in FEniCS. In addition to provide tools for modeling artery networks and computing blood flow, this package may serve as a reference for others who intend to use FEniCS for vascular modeling. In a more general framework, the package may also be used as an example of how to solve equations on conservation form on multiple domains, with shared boundary conditions.

In order for the package to work correctly, the user has to provide appropriate input. The input consists of parameters for the structure of the artery network, as well as intrinsic characteristics for each artery. The user also has to specify parameters for the solver, and output parameters. The only part of the solution that the user has to provide before solving is the inlet flow for the first artery.

The main part of the package is the class Artery_Network, which provides methods for describing a network of arteries, meshing the network, computing boundary conditions and solving the 1D blood flow equations using FEniCS. Each artery is defined according to the Artery class, containing the characteristics and the necessary FEniCS-objects for representing an artery. While the structure of the boundary conditions is taken care of within the Artery class, it is up to the Artery_Network class to compute the right values for the given boundary conditions.

The module provides tools for running, interpreting output and visualising data. An example run-file is preconfigured to read parameters from a config file and to run the necessary functions from Artery_Network in the right order, and thus generating an output folder containing the solution to the 1D blood flow equations. The utils-file provides methods for interpreting output. A post processing file is configured to make plots.

## Approach

To solve the 1D blood flow equations, a Crank-Nicolson finite difference scheme has been implemented in time, and a FEniCS-solver has been used for the spatial variable.

The adimensionalised problem may be written on conservation form:

$$
\frac{\partial U}{\partial t} + \frac{\partial F(U)}{\partial x} = S(U)
$$

where $U = (A, q)$, $F(U) = (q, \frac{q^2}{A} + f(r_0) \sqrt{A_0 A})$, $S(U) = (0, -\frac{2\sqrt{\pi}}{\delta_b Re} \frac{q}{\sqrt{A}} + (2 \sqrt{A} (\sqrt{\pi} f(r_0) + \sqrt{A_0} \frac{df}{dr_0}(r_0)) - A \frac{df}{dr_0}(r_0))\frac{dr_0}{dx})$. A is the cross-sectional area, q is the flow.

The above equation governs the blood flow on one artery. The associated boundary conditions are as follows: The initial radii respect the relation $r_0(x) = Ru (\frac{Rd}{Ru})^{x/L}$. For the root artery, the inlet flow is provided as input. For each artery, the initial flow is computed according to the artery's share of the cross-sectional area at the bifurcation. At each bifurcation point, the area and the flow is computed using Richtmyer's two-step Lax-Wendroff method, by solving a system of equations with 18 unknowns. At the end arteries, the pressure is computed using a Windkessel model. When the pressure is known, the area can be computed from the state equation: $p - p_0 = f(r_0) (1 - \sqrt{\frac{A_0}{A}})$. A more detailed explanation is elaborated in Olufssen:2000, Kolachalama:2007 and Diem:2016.


# Acknowlegdments

The author would like to thank the Simula Research Laboratory for its contributions to this package.


# References
