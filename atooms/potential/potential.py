# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import numpy

class PairPotentialBase(object):

    interacting_bodies = 2

    def __init__(self, name, params, species, cutoff=None, npoints=20000):
        self.name = name
        self.params = params
        self.species = species
        self.cutoff = cutoff
        self.npoints = npoints

        self._lookup = None

    def _tailor(self):
        pass

    def is_zero(self, rsquare):
        self.cutoff.is_zero(rsquare)

    def tabulate(self):
        rcut = 2.5
        rmin = 0.01
        rmax = rcut+0.05
        rsq = numpy.ndarray(self.npoints)
        u0 = numpy.ndarray(self.npoints)
        u1 = numpy.ndarray(self.npoints)
        drsq = rmax**2 / (self.npoints-1)
        #for i in range(self.npoints-1,0,-1):
        for i in range(1,self.npoints):
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

    def compute(self, rsquare):
        raise NotImplementedError()

    def compute_analytical(self, rsquare):
        raise NotImplementedError()

    def write_potential(self, fh, rmin, rmax, nr=100):
        pass


# factory method

__factory_map = {'base' : PairPotentialBase}

def PairPotential(name, *args, **kwargs):

    """ Factory class shortcut """

    if not name in __factory_map:
        return PairPotentialBase(name, *args, **kwargs)
        #raise ValueError('unknown class %s', name)
    else:
        return __factory_map.get(name)(name, *args, **kwargs)
