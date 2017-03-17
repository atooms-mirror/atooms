import numpy

class PairPotentialBase(object):

    interacting_bodies = 2

    def __init__(self, name, params, species, cutoff=None, hard_core=0.0, npoints=20000):
        self.name = name
        self.params = params
        self.species = species
        self.cutoff = cutoff
        self.hard_core = hard_core
        self.npoints = npoints
        self._adjusted = False

    def _adjust(self):
        """Adjust the cutoff to the potential."""
        self._adjusted = True
        if self.cutoff is not None:
            try:
                u = self._compute(self.cutoff.radius**2)
                self.cutoff.tailor(self.cutoff.radius**2, u)
            except NotImplementedError:
                # If _compute() is not implemented, we ignore the error.
                # This way we can use the base potential as a placeholder
                # to move parameters around, e.g. in writing trajectories
                # or forcefields.
                pass              

    def tabulate(self, npoints=None, rmax=None):
        # We tabulate one point more than the cutoff, so that the
        # value for discontinuous potential is not smoothed.
        if npoints is None:
            npoints = self.npoints
        if self.cutoff is not None:
            rmax = self.cutoff.radius + 0.05
        else:
            if rmax is None:
                raise ValueError('rmax is needed to tabulate a cutoff-less potential')
            
        rmin = 0.01
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

    def _compute(self, rsquare):
        raise NotImplementedError()

    def compute(self, rsquare):
        # if rsquare < self.hard_core**2:
        #     u0, u1 = float("inf"), float("inf")
        # else:
        if not self._adjusted:
            self._adjust()
        u = self._compute(rsquare)
        if self.cutoff is not None:
            u = self.cutoff.smooth(rsquare, u)
        return u

    def is_zero(self, rsquare):
        if self.cutoff is not None:
            return self.cutoff.is_zero(rsquare)
        else:
            return False


def lennard_jones(rsq, epsilon, sigma):
    sigsq = sigma**2
    u = 4 * epsilon * ((sigsq/rsq)**6 - (sigsq/rsq)**3)
    w = 24 * epsilon * (2*(sigsq/rsq)**6 - (sigsq/rsq)**3) / rsq
    h = 0.0
    return u, w, h

def NullPotential(PairPotentialBase):

    def _compute(rsq):
        return 0.0, 0.0, 0.0

class LennardJones(PairPotentialBase):

    def _compute(self, rsq):
        return lennard_jones(rsq, **self.params)


# factory method

__factory_map = {'null': NullPotential}

def PairPotential(name, *args, **kwargs):

    """ Factory class shortcut """

    if not name in __factory_map:
        return PairPotentialBase(name, *args, **kwargs)
        #raise ValueError('unknown class %s', name)
    else:
        return __factory_map.get(name)(name, *args, **kwargs)
