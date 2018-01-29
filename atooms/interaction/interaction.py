# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Base interaction class.

Actual backends should implement this interface.
"""

import numpy

class Interaction(object):

    def __init__(self, potential, name=''):
        """
        The interaction is calculated given a set of potentials.

        - `potential` is a list of `Potential` instances.
        - `name` is a string tag that can be used to distinguish
        different interaction instances.
        """
        self.potential = potential
        self.name = name
        self.forces = None
        self.energy = None
        self.virial = None
        self.stress = None  # this will be (ndim,ndim) numpy array
        self.hessian = None

    def compute(self, observable, particle, cell):
        """
        Compute interaction between `particle` instances in a `cell`.

        At a minimum, `observable` can take the following values:
        `energy`, `forces`, `stress`. Note that each of these
        observables imply the ones preceeding it, e.g. computing
        forces implies energy calculation. The following observables
        are set to `None`.
        """
        if observable == 'energy':
            self.energy = 0.0
            self.virial = None
            self.stress = None
            self.forces = None
        elif observable == 'forces':
            self.energy = 0.0
            self.virial = 0.0
            self.stress = None
            self.forces = numpy.zeros((len(particle), len(cell.side)))
        elif observable == 'stress':
            self.energy = 0.0
            self.virial = 0.0
            self.stress = numpy.zeros((len(cell.side), len(cell.side)))
            self.forces = numpy.zeros((len(particle), len(cell.side)))
        else:
            raise ValueError('unsupported observable %s' % observable)

    def report(self):
        txt = ''
        for p in self.potential:
            try:
                txt += p.report()
            except AttributeError:
                pass
        return txt

