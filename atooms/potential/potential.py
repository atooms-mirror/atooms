# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import numpy

class PairPotentialBase(object):

    interacting_bodies = 2

    def __init__(self, name, params, species, cutoff=None, hard_core=0.0, npoints=20000):
        self.name = name
        self.params = params
        self.species = species
        self.cutoff = cutoff
        self.hard_core = hard_core
        self._tailor()

    def _tailor(self):
        u0, u1 = self._compute(self.cutoff.radius**2)
        self.cutoff.tailor(self.cutoff.radius**2, u0, u1)

    def is_zero(self, rsquare):
        return self.cutoff.is_zero(rsquare)

    def tabulate(self, npoints=None):
        if npoints is not None:
            self.npoints = npoints
        rcut = self.cutoff.radius
        rmin = 0.01
        rmax = rcut+0.05
        rsq = numpy.ndarray(npoints)
        u0 = numpy.ndarray(npoints)
        u1 = numpy.ndarray(npoints)
        drsq = rmax**2 / (npoints-1)
        #for i in range(self.npoints-1,0,-1):
        for i in range(1,npoints):
            rsq[i] = i*drsq
            if not self.is_zero(rsq[i]):
                u0[i], u1[i] = self.compute(rsq[i])
            else:
                u0[i], u1[i] = 0, 0
        rsq[0] = 0.0
        u0[0] = max(u0)
        u1[0] = max(u1)
        # for i in range(self.npoints):
        #     if math.isnan(u0[i]):
        #         u0[i] = max(u0)
        #     if math.isnan(u1[i]):
        #         u1[i] = max(u1)
        return rsq, u0, u1

    def _compute(self, rsquare):
        raise NotImplementedError()

    def compute(self, rsquare):
        if rsquare < self.hard_core**2:
            u0, u1 = float("inf"), float("inf")
        u0, u1 = self._compute(rsquare)
        u0, u1 = self.cutoff.smooth(rsquare, u0, u1)
        return u0, u1

def NullPotential(PairPotentialBase):

    def _compute(rsquare):
        return 0.0, 0.0

# factory method

__factory_map = {'null': NullPotential}

def PairPotential(name, *args, **kwargs):

    """ Factory class shortcut """

    if not name in __factory_map:
        return PairPotentialBase(name, *args, **kwargs)
        #raise ValueError('unknown class %s', name)
    else:
        return __factory_map.get(name)(name, *args, **kwargs)
