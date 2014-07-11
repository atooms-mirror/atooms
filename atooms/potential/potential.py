# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

# TODO: dangling fix needed for import of potentials in trajectory

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

    def _tabulate(self):
        pass

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
