# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Trajectory decorators."""

import numpy

class Centered(object):

    """ Center positions in the box on the fly """

    def __new__(cls, component):
        cls = type('Centered', (Centered, component.__class__), component.__dict__)
        return object.__new__(cls)

    def __init__(self, component):
        pass

    def read_sample(self, sample):
        # TODO: check that we have not subtracted it yet (use a list)
        s = super(Centered, self).read_sample(sample)
        for p in s.particle:
            p.position -= s.cell.side / 2.0
        return s

class Sliced(object):

    """ Only return a slice of a trajectory """

    # This is still necessary. slicing via __getitem__ has a large memory fingerprint
    # since we couldnt write it as a generator (maybe it is possible?)
    # TODO: adjust uslice to pick up blocks without truncating them!

    def __new__(cls, component, uslice):
        cls = type('Sliced', (Sliced, component.__class__), component.__dict__)
        return object.__new__(cls)

    def __init__(self, component, uslice):
        self._sliced_samples = range(len(self.steps))[uslice]
        self.steps = self.steps[uslice]

    def read_sample(self, sample):
        i = self._sliced_samples[sample]
        return super(Sliced, self).read_sample(i)
        

class Unfolded(object):

    """ Decorate Trajectory to unfold particles positions on the fly """

    def __new__(cls, component):
        cls = type('Unfolded', (Unfolded, component.__class__), component.__dict__)
        return object.__new__(cls)

    def __init__(self, component):
        pass
    
    def read_init(self):
        s = super(Unfolded, self).read_init()
        # Cache the initial sample and cell
        s = super(Unfolded, self).read_sample(0)
        self._old = numpy.array([p.position for p in s.particle])
        self._last_read = 0

    def read_sample(self, sample):
        # Return here if first sample
        if sample == 0:
            return super(Unfolded, self).read_sample(sample)

        # Compare requested sample with last read
        delta = sample - self._last_read
        if delta < 0:
            raise ValueError('cannot unfold jumping backwards (delta=%d)' % delta)
        if delta > 1:
            # Allow to skip some samples by reading them internally
            # We read delta-1 samples, then delta is 1
            for i in range(delta-1):
                self.read_sample(self._last_read+1)

        s = super(Unfolded, self).read_sample(sample)
        self._last_read = sample
      
        # Unfold positions
        # Note that since L can be variable we get it at each step
        # TODO: I am not entirely sure this is correct with NPT.
        # The best thing in this case is to get unfolded positions
        # from the simulation.
        L = s.cell.side
        pos = numpy.array([p.position for p in s.particle])
        dif = pos - self._old
        dif = dif - numpy.rint(dif/L) * L
        self._old += dif

        # Return unfolded system
        for i in xrange(len(pos)):
            s.particle[i].position = self._old[i][:]
        return s

# TODO: see if we can avoid reading anything on construction
# TODO: how to better handle conversions between subclasses?
# We cannot use _convert() because we rely on close() method being called

# To properly implement decorators in python see 
# http://stackoverflow.com/questions/3118929/implementing-the-decorator-pattern-in-python
# asnwer by Alec Thomas. if we don't subclass at runtime we won't be able to use the decorated
# mathod in other non-subclassed methods.

class MatrixFix(object):

    # Subclass component at runtime
    def __new__(cls, component, matrix_species):
        cls = type('MatrixFix', (MatrixFix, component.__class__), component.__dict__)
        return object.__new__(cls)

    # Object initialization
    def __init__(self, component, matrix_species):
        self._component = component # actually unnecessary
        self.matrix_species = matrix_species

    # def _init_read(self):
    #     s = super(MatrixFix, self)._init_read()
    #     matrix = []
    #     fluid = []
    #     for p in s.particle:
    #         if p.id in self.matrix_species:
    #             matrix.append(p)
    #         else:
    #             fluid.append(p)
    #     s.particle = fluid
    #     s.matrix = matrix
    #     return s

    def read_sample(self, *args, **kwargs):
        s = super(MatrixFix, self).read_sample(*args, **kwargs)
        # Get rid of matrix particles in trajectory
        s.particle = [p for p in s.particle if not p.id in self.matrix_species]
        return s

# TODO: can we write a base Trajectory Decorator like this?
class TrajectoryDecorator(object):

    # Subclass component at runtime
    def __new__(cls, component, cb_read_initial_state):
        cls = type('TrajectoryDecorator', (TrajectoryDecorator, component.__class__), component.__dict__)
        return object.__new__(cls)

    # Object initialization
    def __init__(self, component, cb_read_initial_state):
        self._component = component # actually unnecessary
        self._cb_read_initial_state = cb_read_initial_state
        #self._cb_read = cb_read

    def _init_read(self):
        s = super(TrajectoryDecorator, self)._init_read()
        self._cb_read_initial_state(self, s)
        return s

    # def read(self, sample):
    #     s = super(TrajectoryDecorator, self).read(sample)
    #     self._cb_read(s, sample)
    #     return s


class NormalizeId(object):

    """Make sure all chemical ids start from 1"""

    def __new__(cls, component):
        cls = type('NormalizeId', (NormalizeId, component.__class__), component.__dict__)
        return object.__new__(cls)

    def __init__(self, component): 
        pass

    def _normalize(self, pl):
        pid = [p.id for p in pl]
        id_min = numpy.min(pid)
        if id_min == 0:
            for p in pl:
                p.id += 1
        return pl

    def read_sample(self, *args, **kwargs):
        s = super(NormalizeId, self).read_sample(*args, **kwargs)
        s.particle = self._normalize(s.particle)
        return s


class MatrixFlat(object):

    # Subclass component at runtime
    def __new__(cls, component):
        cls = type('MatrixFlat', (MatrixFlat, component.__class__), component.__dict__)
        return object.__new__(cls)

    # Object initialization
    def __init__(self, component):
        self._component = component # actually unnecessary
        self._matrix = None
        
    def __setup_matrix(self, s):

        if self._matrix is not None:
            return

        infinite_mass = 1e20
        max_isp = max([p.id for p in s.particle])

        self._matrix = []
        for m in s.matrix:
            self._matrix.append(m)
            self._matrix[-1].mass = infinite_mass
            self._matrix[-1].id += max_isp

        # Sort matrix particles by index
        self._matrix.sort(key = lambda a : a.id)

    def read_sample(self, *args, **kwargs):
        s = super(MatrixFlat, self).read_sample(*args, **kwargs)
        self.__setup_matrix(s)
        for m in self._matrix:
            s.particle.append(m)
        return s


class Filter(object):
    
    """Apply a filter that transforms each read sample in a trajectory"""

    def __new__(cls, component, filt, *args):
        cls = type('Filter', (Filter, component.__class__), component.__dict__)
        return object.__new__(cls)

    def __init__(self, component, filt, *args):
        """filt is a function that receives a System and returns a modified version of it"""
        import copy
        self.filt = filt
        self._args = args

    def read_sample(self, sample):
        sy = super(Filter, self).read_sample(sample)
        # Apply filter to the system, that's all        
        # HACK!: when further decorating the class, the referenced function becomes a bound method of the decorated class and therefore passes the first argument (self)
        # Workaround is to always expect a trajectory and forcibly add self if no further decorators are added
        try:
            return self.filt(sy, *self._args)
        except:
            return self.filt(self, sy, *self._args)
            
