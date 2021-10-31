# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Base interaction class.

Actual backends should implement this interface.
"""

import numpy


class Interaction(object):

    def __init__(self):
        self.variables = {'position': 'particle.position',
                          'species': 'particle.species:int32',
                          'side': 'cell.side'}
        """
        A list of variables needed to compute the interaction. The keys
        must match the interfaces of compute(), the fields are
        variables accepted by System.dump(). It is possible to specify
        the required data type using the optionl colon syntax
        <property>[:<dtype>]. The dtype must a valid identifier for
        numpy array creation.
        """
        self.order = 'F'        
        self.forces = None
        self.energy = None
        self.virial = None
        self.stress = None  # this will be (ndim,ndim) numpy array
        self.hessian = None
        
    def compute(self, observable, position=None, species=None, side=None):
        """
        Compute interaction between `particle` instances in a `cell`.

        At a minimum, `observable` can take the following values:
        `energy`, `forces`, `stress`. Note that each of these
        observables imply the ones preceeding it, e.g. computing
        forces implies energy calculation. The following observables
        are set to `None`.
        """
        assert position is not None
        assert species is not None
        assert side is not None
        if observable == 'energy':
            self.energy = 0.0
            self.virial = None
            self.stress = None
            self.forces = None
        elif observable == 'forces':
            self.energy = 0.0
            self.virial = 0.0
            self.stress = None
            self.forces = numpy.zeros_like(position)
        elif observable == 'stress':
            self.energy = 0.0
            self.virial = 0.0
            self.stress = numpy.zeros((len(side), len(side)))
            self.forces = numpy.zeros_like(position)
        else:
            raise ValueError('unsupported observable %s' % observable)
