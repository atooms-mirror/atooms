# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""Simulation cell."""

import numpy
import warnings
from atooms.core import ndim as _ndim


class Cell(object):

    def __init__(self, side=None, center=None, periodic=None):
        if side is None:
            self.side = numpy.zeros(_ndim)
        else:
            self.side = numpy.asarray(side, dtype=numpy.float64)
        # By default, the center is at the origin of the reference frame
        if center is None:
            self.center = numpy.zeros_like(self.side)
        else:
            self.center = numpy.asarray(center, dtype=numpy.float64)
        self.shape = 'cubic'

        # Periodic boundary conditions apply separately on each axis
        if periodic is None:
            self.periodic = numpy.ndarray(self.side.size, dtype=bool)
            self.periodic[:] = True
        else:
            self.periodic = numpy.asarray(periodic, dtype=bool)

    @property
    def volume(self):
        return numpy.prod(self.side)

    @property
    def origin(self):
        warnings.warn('cell origin is deprecated', DeprecationWarning)
        return numpy.zeros(_ndim)
