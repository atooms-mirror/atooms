"""
Library of pair potentials.

`PairPotential` instances are built on top of plain functions that
return the energy and its derivatives. In this module we collect some
common pair potentials.

Aliases are defined for convenience. These are all equivalent calls

    lennard_jones(1.0, epsilon=1.0, sigma=1.0)
    lj(1.0, epsilon=1.0, sigma=1.0)
    LJ(1.0, epsilon=1.0, sigma=1.0)

"""

from math import sqrt

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

# Aliases

lj = lennard_jones
LJ = lennard_jones

def square_well(rsq, epsilon, sigma):
    """Square well potential."""
    if rsq > sigma**2:
        return 0.0, 0.0, 0.0
    else:
        return epsilon, 0.0, 0.0

SW = square_well
sw = square_well

def harmonic_sphere(self, rsq, epsilon, sigma):
    """Harmonic sphere potential.

    u(r) = 0.5 * epsilon * [1-(r/sigma)**2]
    """
    r = sqrt(rsq)
    return 0.5 * epsilon * (1.0 - r/sigma)**2, \
        epsilon * (1.0 - r/sigma) / (sigma*r), \
        0.0

harm = harmonic_sphere
