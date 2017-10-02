"""
Library of pair potentials.

`PairPotential` instances are built on top of plain functions that
return the energy and its derivatives. In this module we collect some
common pair potentials.

Example:
-------

The Lennard-Jones potential calculated at a reduced squared distance
equal to 1.0:

    lj = lennard_jones(1.0, epsilon=1.0, sigma=1.0)
"""

from math import sqrt

__all__ = ['constant', 'lennard_jones', 'harmonic_sphere']

def constant(rsq, epsilon):
    """Constant potential."""
    return epsilon, 0.0, 0.0

def lennard_jones(rsq, epsilon, sigma):
    """
    Lennard-Jones potential.

    u(r) = 4 * epsilon * [(r/sigma)^12 - (r/sigma)^6]
    """
    sigsq = sigma**2
    u = 4 * epsilon * ((sigsq/rsq)**6 - (sigsq/rsq)**3)
    w = 24 * epsilon * (2*(sigsq/rsq)**6 - (sigsq/rsq)**3) / rsq
    h = 0.0
    return u, w, h

def harmonic_sphere(rsq, epsilon, sigma):
    """
    Harmonic sphere potential.

    u(r) = 0.5 * epsilon * [1-(r/sigma)]**2 if r<=sigma else 0
    """
    r = sqrt(rsq)
    return 0.5 * epsilon * (1.0 - r/sigma)**2, \
        epsilon * (1.0 - r/sigma) / (sigma*r), \
        0.0
