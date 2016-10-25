# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import numpy
from atooms.core import ndim

class Cell(object):

    """Cell class"""
    
    def __init__(self, side, origin = numpy.zeros(ndim)):
        self.side   = numpy.array(side)
        self.origin = origin

    @property
    def volume(self):
        return numpy.prod(self.side)
