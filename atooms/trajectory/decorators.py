# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Trajectory callbacks and class decorators.

- "Callbacks" are simple functions that modify the `system` instance
  returned by `read_sample()`. They can be registered to `trajectory`
  instance via `add_callback()`.
- "Class decorators" can be used for
  more complex modifications of trajectory behavior.
"""

import numpy

__all__ = ['center', 'normalize_id', 'sort', 'filter_id',
           'set_density', 'set_temperature', 'fix_cm', 'fold',
           'Sliced', 'Unfolded']

# Callbacks

def center(system):
    """
    Center particles in the simulation cell.

    It won't check if that is done multiple times.
    """
    for p in system.particle:
        p.position -= system.cell.side / 2.0
    return system

def normalize_id(system, alphabetic=False):
    """
    Change particle species id's to start from 1 (fortran style).

    Species names, given by the variable `particle.name`, can be
    reassigned alphabetically (ex. 'A', 'B' ...)  using the
    `alphabetic` flag.
    """
    map_ids = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E'}
    pid = [p.id for p in system.particle]
    id_min = numpy.min(pid)
    if id_min == 0:
        for p in system.particle:
            p.id += 1
    if alphabetic:
        for p in system.particle:
            p.name = map_ids[p.id]
    return system

def sort(system):
    """Sort particles by species id."""
    return sorted(system.particle, key=lambda a: a.id)

def filter_id(system, species):
    """Return particles of a given `species` id."""
    system.particle = [p for p in system.particle if p.id == species]
    return system

def set_density(system, rho):
    """Set density of system to `rho` by rescaling the cell."""
    rho_old = system.density
    x = (rho_old / rho)**(1./3)
    system.cell.side *= x
    for p in system.particle:
        p.position *= x
    return system

def set_temperature(system, T):
    """Set system temperature to `T` by reassigning velocities."""
    from atooms.system.particle import cm_velocity
    for p in system.particle:
        p.maxwellian(T)
    v_cm = cm_velocity(system.particle)
    for p in system.particle:
        p.velocity -= v_cm
    return system

def fix_cm(s):
    # Get current position of CM from unfolded positions
    cm = s.cm_position
    for p in s.particle:
        p.position -= cm
    return s

def fold(s):
    """Center and fold positions into central cell."""
    for p in s.particle:
        p.position -= s.cell.side / 2
        p.fold(s.cell)
    return s


# Class decorators

# To properly implement decorators in python see
# http://stackoverflow.com/questions/3118929/implementing-the-decorator-pattern-in-python
# asnwer by Alec Thomas. if we don't subclass at runtime we won't be able to use the decorated
# mathod in other non-subclassed methods.

class Sliced(object):

    """Only return a slice of a trajectory."""

    # This is still necessary. slicing via __getitem__ has a large memory fingerprint
    # since we couldnt write it as a generator (maybe it is possible?)
    # TODO: adjust uslice to pick up blocks without truncating them

    def __new__(cls, component, uslice):
        cls = type('Sliced', (Sliced, component.__class__), component.__dict__)
        return object.__new__(cls)

    def __init__(self, component, uslice):
        self._sliced_frames = range(len(self.steps))[uslice]
        self.steps = self.steps[uslice]

    def read_sample(self, frame):
        i = self._sliced_frames[frame]
        return super(Sliced, self).read_sample(i)


class Unfolded(object):

    """Decorate Trajectory to unfold particles positions on the fly."""

    def __new__(cls, component, fixed_cm=False):
        cls = type('Unfolded', (Unfolded, component.__class__), component.__dict__)
        return object.__new__(cls)

    def __init__(self, component, fixed_cm=False):
        self._initialized_read = False
        self.fixed_cm = fixed_cm

    def read_init(self):
        s = super(Unfolded, self).read_init()
        # Cache the initial sample and cell
        s = super(Unfolded, self).read_sample(0)
        self._old = numpy.array([p.position for p in s.particle])
        self._old_cm = s.cm_position
        self._last_read = 0

    def read_sample(self, frame):
        # Return here if first frame
        if frame == 0:
            return super(Unfolded, self).read_sample(frame)

        # Compare requested frame with last read
        delta = frame - self._last_read
        if delta < 0:
            raise ValueError('cannot unfold jumping backwards (delta=%d)' % delta)
        if delta > 1:
            # Allow to skip some frames by reading them internally
            # We read delta-1 frames, then delta is 1
            for i in range(delta-1):
                self.read_sample(self._last_read+1)

        s = super(Unfolded, self).read_sample(frame)
        self._last_read = frame

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
