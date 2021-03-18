import os
import logging
import numpy
import f2py_jit


_log = logging.getLogger(__name__)


class VerletList(object):

    def __init__(self, skin=0.2, update="largest", update_period=10, source=os.path.join(os.path.dirname(__file__),'neighbor_list.f90')):
        self.skin = skin
        self.update = update
        self.update_period = update_period
        self.displacement = None
        self.rcut = None
        self.number_of_neighbors = None
        self.neighbors = None
        self.is_adjusted = False
        self.calls = 0
        self._last_call = 0        
        self._uid = f2py_jit.build_module(source)
        self._f90 = f2py_jit.import_module(self._uid)

    def adjust(self, box, pos, rcut, boost=1.0):
        """
        Adjust neighbor list to positions of particles in the box and cut off distance.
        """
        from math import gamma, pi
        
        if self.is_adjusted:
            return
        self.rcut = rcut + self.skin
        ndim, N = pos.shape
        volume = pi**(ndim/2.) / gamma(ndim/2. + 1) * numpy.max(self.rcut)**ndim
        rho = N / numpy.product(box)
        nmax = int(rho * volume * 1.5)
        self.displacement = numpy.zeros_like(pos, order='F')
        self.number_of_neighbors = numpy.zeros(N, dtype=numpy.int32, order='F')
        self.neighbors = numpy.zeros((nmax, N), dtype=numpy.int32, order='F')
        self._pos_old = pos.copy(order='F')        
        self.is_adjusted = True

    def add_displacement(self, pos, box):
        """
        Accumulate displacements assuming that PBC has been applied, hence we need to refold them.
        """
        self._f90.neighbor_list.add_displacement(pos, self._pos_old, box, self.displacement)

    def need_update(self):
        """Return True if neighbor list must be completely updated."""
        if self.update == 'periodic':
            return self.calls % self.update_period == 0
        elif self.update == 'largest':
            return self._f90.neighbor_list.need_update_largest(self.displacement, self.skin)
        else:
            raise ValueError('unknown update method {}',format(self.update))

    def compute(self, box, pos, ids):
        """Compute Verlet lists"""        
        self.add_displacement(pos, box)
        if self.calls == 0 or self.need_update():
            self._f90.neighbor_list.compute(box, pos, ids, self.rcut, self.neighbors, self.number_of_neighbors)
            self.displacement[:, :] = 0.
            _log.debug('updated list after {} callls'.format(self.calls - self._last_call))
            self._last_call = self.calls
        self.calls += 1
