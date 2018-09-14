# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""Simulation cell."""

import numpy
from atooms.core import ndim as _ndim


class Cell(object):

    def __init__(self, side=None, origin=None):
        if side is None:
            self.side = numpy.zeros(_ndim)
        else:
            self.side = numpy.asarray(side, dtype=numpy.float64)
        if origin is None:
            self.origin = numpy.zeros(_ndim)
        else:
            self.origin = numpy.asarray(side, dtype=numpy.float64)          
        self.shape = 'cubic'

    @property
    def volume(self):
        return numpy.prod(self.side)
