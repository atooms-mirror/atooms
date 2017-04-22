# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Individual particles interact via a `Potential`. A simple example is a
two-body potential that depends only on the scalar distance between
two particles. The potentials can be cut off and smoothed by a
`CutOff`. `Interaction` accounts for the total interaction of all the
particles in a system (order N^2 calculation).
"""

class InteractionBase(object):

    def __init__(self, name, potential):
        self.name = name
        self.potential = potential
        self.total_energy = 0.0
        self.total_virial = 0.0
        self.total_stress = None  # this will be ndim,ndim dimensional numpy array

# InteractionTerm is just an interface and does not need any python implementation.
# System will simply return a list : system.interaction_terms

# factory method

__factory_map = {'base': InteractionBase}

def Interaction(name, *args, **kwargs):

    """ Factory class shortcut """

    if name not in __factory_map:
        return InteractionBase(name, *args, **kwargs)
        # raise ValueError('unknown class %s', name)
    else:
        return __factory_map.get(name)(name, *args, **kwargs)
