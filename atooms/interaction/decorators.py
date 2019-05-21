"""
Potential decorators to produce hard pair potentials.

A purely hard potential has hardcore > rcut. When this condition is
met, client code can perform some optimizations, ex. tabulation is
skipped.
"""

from math import sqrt

__all__ = ['hard_sphere', 'square_well']

from .library import constant
from .cutoff import CutOff


def hard_sphere(potential):
    """
    Hard sphere potential

    Parameters: sigma
    """
    potential.func = constant
    potential.params['epsilon'] = 0.0
    potential.hard_core = potential.params['sigma']
    potential.cutoff = CutOff('c', 0.0)
    return potential


def square_well(potential):
    """
    Square well potential

    Parameters: epsilon, sigma, delta
    """
    potential.func = constant
    potential.hard_core = potential.params['sigma']
    potential.cutoff = CutOff('c', potential.params['sigma'] + potential.params['delta'])
    return potential
