# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""Pair potential classes and factory."""

import warnings
import numpy

from atooms.interaction import library
from atooms.interaction import decorators

def tabulate(potential, parameters, cutoff='c', rc=2.5, npoints=10000,
             rmin=0.5, fmt='lammps', metadata='', fileout=None, precision=14):

    """Tabulate a potential."""

    from atooms.core.utils import tipify
    from atooms.interaction.potential import PairPotential
    from atooms.interaction.cutoff import CutOff

    if isinstance(parameters, dict):
        param_dict = parameters
    else:
        param_dict = {}
        for param in parameters.split(','):
            key, value = param.split('=')
            param_dict[key] = tipify(value)

    potential = PairPotential(potential, param_dict, (1, 1))
    if cutoff is not None:
        potential.cutoff = CutOff(cutoff, rc)
    rsq, u0, u1, u2 = potential.tabulate(npoints, rmin=rmin, what='uwh')
    r = rsq**0.5
    if fmt == 'lammps':
        u1 *= r
        txt = """

POTENTIAL
N {}

""".format(len(rsq))
        i = 1
        for x, y, z in zip(r, u0, u1):
            txt += '{} {} {} {}\n'.format(i, x, y, z)
            i += 1

    elif fmt == 'uwh':
        txt = '# {} columns: r, u, w, h\n'.format(metadata)
        for x, y, z, w in zip(r, u0, u1, u2):
            txt += '{:.14g} {:.14g} {:.14g} {:.14g}\n'.format(x, y, z, w)

    else:
        u1 *= r
        txt = '# {} columns: r, u, f\n'.format(metadata)
        for x, y, z in zip(r, u0, u1):
            txt += '{} {} {}\n'.format(x, y, z)

    if fileout is None:
        return txt
    else:
        with open(fileout, 'w') as fh:
            fh.write(txt)


class PairPotential(object):

    """Pair potential between two particles."""

    interacting_bodies = 2

    def __init__(self, func, params, species, cutoff=None,
                 hard_core=0.0, npoints=20000):
        """
        If `func` is a function, it will be used it to compute the
        potential u(r) and its derivatives. `func` takes the squared
        distance between the particles as a first argument plus an
        arbitrary number of keyword arguments. The `params` dict will
        be passed to the function. It must return the tuple (u, w, h)
        where

        - u = u(r)
        - w = - (du/dr)/r (w=-W/r^2 where W is the virial)
        - h = - (dw/dr)/r (term for the Hessian matrix)

        At present only u, w are used.

        If `func` is a string, it will be looked up into a database of
        pair potentials (the `_factory` dict). The database is loaded
        with the functions found in the `pair_potentials` module.

        `species` is a list or tuple of size 2 containing the
        particles' id associated to the potential.

        Examples:
        --------
        The Lennard-Jones potential:

        `PairPotential('lennard_jones', {'epsilon': 1.0, 'sigma': 1.0}, [1, 1])`
        """
        self.func = func
        self.params = params
        self.species = species
        self.cutoff = cutoff
        self.hard_core = hard_core
        self.npoints = npoints
        self._adjusted = False

        if not hasattr(self.func, '__call__'):
            # If func is not callable, look up the potential in the
            # potential library
            if self.func in library.__all__:
                self.func = getattr(library, self.func)
            elif self.func in decorators.__all__:
                # Decorate the potential
                # Used to implement hard potentials
                getattr(decorators, self.func)(self)
            else:
                raise ValueError('unknown potential %s' % self.func)

    def __str__(self):
        if type(self.func) is str:
            return self.func
        else:
            return self.func.__name__

    def report(self):
        txt = """\
potential {0.species}: {0.func.__name__}
parameters: {0.params}
cutoff: {0.cutoff} at {0.cutoff.radius}
""".format(self)
        if self.hard_core > 0:
            txt += "hardcore: {0.hard_core}\n".format(self)
        return txt

    def _adjust(self):
        """Adjust the cutoff to the potential."""
        self._adjusted = True
        if self.cutoff is not None and self.cutoff.radius > 0:
            u = self.func(self.cutoff.radius**2, **self.params)
            self.cutoff.tailor(self.cutoff.radius**2, u)

    def tabulate(self, npoints=None, rmax=None, rmin=0.0, what='uw'):
        """
        Tabulate the potential from 0 to `rmax`.

        The potential cutoff is only used to determine `rmax` if this
        is not given. The full potential is tabulated, it is up to the
        calling code to truncate it. We slightly overshoot the
        tabulation, to avoid boundary effects at the cutoff or at
        discontinuities.

        The `what` parameters can be 'u', 'uw', 'uwh'.
        """
        if not self._adjusted:
            self._adjust()

        if npoints is None:
            npoints = self.npoints
        if self.cutoff is None:
            if rmax is None:
                raise ValueError('rmax is needed to tabulate a cutoff-less potential')
        else:
            rmax = self.cutoff.radius

        rsq = numpy.ndarray(npoints)
        u0 = numpy.ndarray(npoints)
        u1 = numpy.ndarray(npoints)
        u2 = numpy.ndarray(npoints)
        # We overshoot 2 points beyond rmax (cutoff) to avoid
        # smoothing discontinuous potentials.
        # This is necessary also for the Allen Tildesley lookup table,
        # which for any distance within the cutoff will look up two
        # points forward in the table.
        # Note that the cutoff is applied to the function only to smooth it
        # not to cut it.
        drsq = (rmax**2 - rmin**2) / (npoints - 3)
        warnings.filterwarnings("ignore")
        for i in range(npoints):
            rsq[i] = rmin**2 + i * drsq
            try:
                u0[i], u1[i], u2[i] = self.compute(rsq[i])
            except ZeroDivisionError:
                u0[i], u1[i], u2[i] = float('nan'), float('nan'), float('nan')
        warnings.resetwarnings()

        # For potentials that diverge at zero, we remove the singularity by hand
        import math
        if math.isnan(u0[0]):
            u0[0], u1[0], u2[0] = u0[1], u1[1], u2[0]
        if 'h' in what:
            return rsq, u0, u1, u2
        else:
            return rsq, u0, u1

    def compute(self, rsquare):
        """Compute the potential and its derivatives."""
        if not self._adjusted:
            self._adjust()
        # Compute the potential and smooth it via the cutoff
        u = self.func(rsquare, **self.params)
        if self.cutoff is not None:
            u = self.cutoff.smooth(rsquare, u)
        # if rsquare < self.hard_core**2:
        #     u0, u1 = float("inf"), float("inf")
        return u

    def is_zero(self, rsquare):
        """
        Returns `True` if `rsquare` is beyond the squared cutoff distance.
        """
        if self.cutoff is not None:
            return self.cutoff.is_zero(rsquare)
        else:
            return False
