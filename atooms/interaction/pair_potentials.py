"""Pair potentials."""

from math import sqrt

def lennard_jones(rsq, epsilon, sigma):
    sigsq = sigma**2
    u = 4 * epsilon * ((sigsq/rsq)**6 - (sigsq/rsq)**3)
    w = 24 * epsilon * (2*(sigsq/rsq)**6 - (sigsq/rsq)**3) / rsq
    h = 0.0
    return u, w, h

# Aliases

lj = lennard_jones
LJ = lennard_jones

def square_well(rsq, epsilon, sigma):
    if rsq > sigma**2:
        return 0.0, 0.0, 0.0
    else:
        return epsilon, 0.0, 0.0

SW = square_well
sw = square_well

def harmonic_sphere(self, rsq, epsilon, sigma):
    r = sqrt(rsq)
    return 0.5 * epsilon * (1.0 - r/sigma)**2, \
        epsilon * (1.0 - r/sigma) / (sigma*r), \
        0.0

harm = harmonic_sphere
