# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import numpy

from cutoff import CutOff
from pair_potentials import *

# Potentials factory

_factory = {}

def update(module, factory):
    import sys
    import inspect
    for name, func in inspect.getmembers(sys.modules[module], inspect.isfunction):
        if name is not 'update':
            factory[name] = func

update(__name__, _factory)


class PairPotential(object):

    interacting_bodies = 2

    def __init__(self, func, params, species, cutoff=None, hard_core=0.0, npoints=20000):
        self.func = func
        self.params = params
        self.species = species
        self.cutoff = cutoff
        self.hard_core = hard_core
        self.npoints = npoints
        self._adjusted = False
        if self.func in _factory:
            self.func = _factory[func]

    def __str__(self):
        if type(self.func) is str:
            return self.func
        else:
            return self.func.__name__

    def _adjust(self):
        """Adjust the cutoff to the potential."""
        self._adjusted = True
        if self.cutoff is not None:
            u = self.func(self.cutoff.radius**2, **self.params)
            self.cutoff.tailor(self.cutoff.radius**2, u)

    def tabulate(self, npoints=None, rmin=0.01, rmax=None):
        # We tabulate one point more than the cutoff, so that the
        # value for discontinuous potential is not smoothed.
        if npoints is None:
            npoints = self.npoints
        if self.cutoff is not None:
            rmax = self.cutoff.radius + 0.05
        else:
            if rmax is None:
                raise ValueError('rmax is needed to tabulate a cutoff-less potential')
            
        rsq = numpy.ndarray(npoints)
        u0 = numpy.ndarray(npoints)
        u1 = numpy.ndarray(npoints)
        drsq = rmax**2 / (npoints-1)

        rsq[0], rsq[1] = 0.0, drsq
        u0[0], u1[0], _ = self.compute(rsq[1])
        u0[1], u1[1], _ = self.compute(rsq[1])
        for i in range(2,npoints):
            rsq[i] = i*drsq
            if self.is_zero(rsq[i-2]):
                u0[i], u1[i], _ = 0, 0
            else:
                u0[i], u1[i], _ = self.compute(rsq[i])
        return rsq, u0, u1

    def compute(self, rsquare):
        if not self._adjusted:
            self._adjust()
        # Compute the potential and smooth it
        u = self.func(rsquare, **self.params)
        if self.cutoff is not None:
            u = self.cutoff.smooth(rsquare, u)
        # if rsquare < self.hard_core**2:
        #     u0, u1 = float("inf"), float("inf")
        return u

    def is_zero(self, rsquare):
        if self.cutoff is not None:
            return self.cutoff.is_zero(rsquare)
        else:
            return False
