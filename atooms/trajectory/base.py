# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import warnings

from .utils import get_period

class TrajectoryBase(object):

    """Trajectory base class.

    __init__ is supposed to deal with file existence, creating
    handles, setup steps list.

    Read and write implement the following template.

    1. read_init() and write_init() are called only once to initialize
    data structures (ex. counts samples and steps) or grab metadata
    (stuff that doesn't change)

    2. read_sample() and write_sample() are used to actually
    read/write a system

    Additionally, write_sample append the step to the step list

    In future implementation, we might pass a list of objects to be
    written, to store for instance, integrator data and so on.
    """

    # TODO: there is a problem with putting metatdata reading in read_init. It means that steps and timestep are not known before calling read(). These should then be properties that get initialized by calling read_init, rather than their specific read_timestep, read_steps methods

    # We might consider renaming, although it is a bad idea.
    # What init methods are supposed to be is parsing / writing metadata.
    # read_init -> read_metadata
    # write_init -> write_metadata

    # metadata is:
    # read: dt, steps, cell (if invariant). Subclasses may have additional simulation info: integration algorithm etc
    # write, dt, cell (if invariant)

    # steps wants to become a property then.

    # in xyz, rename read_metadata -> read_header

    suffix = None

    def __init__(self, filename, mode='r'):
        """When mode is 'r', it must set the list of available steps."""
        self.filename = filename
        self.mode = mode
        # fmt is a list of strings describing data to be written by
        # write_sample(). Subclasses may use it to filter out some
        # data from their format or can even ignore it entirely.
        self.fmt = []
        self.precision = 6
        self.steps = []
        # These are cached properties
        self._grandcanonical = None
        self._timestep = None
        self._block_period = None
        # Internal state
        self._initialized_write = False
        self._initialized_read = False
        # Sanity checks
        if self.mode == 'r' and not os.path.exists(self.filename):
            raise IOError('trajectory file %s does not exist' % self.filename)

    def __len__(self):
        return len(self.steps)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __iter__(self):
        for i in xrange(len(self.steps)):
            yield self.read(i)

    def __getitem__(self, key):
        if isinstance(key, slice):
            # This works but it loads the whole trajectory in ram.
            # The Sliced decorator doesn't have this issue.
            # If we make this a generator, then access a single sample
            # wont work. Unless we put it in separate functions?
            samples = range(len(self.steps))
            return [self.read(i) for i in samples[key]]

        elif isinstance(key, int):
            if key < 0:
                key += len(self)
            if key >= len(self):
                raise IndexError("Index (%d) is out of range." % key)
            return self.read(key)

        else:
            raise TypeError("Invalid argument type [%s]" % type(key))

    def close(self):
        pass

    def read(self, index):
        if not self._initialized_read:
            self.read_init()
            self._initialized_read = True
        return self.read_sample(index)

    def write(self, system, step):
        if self.mode == 'r':
            raise IOError('trajectory file not open for writing')
        if not self._initialized_write:
            self.write_init(system)
            self._initialized_write = True
        self.write_sample(system, step)
        # Step is added last, sample index starts from 0 by default
        self.steps.append(step)

    def read_init(self):
        """It may setup data structures needed by the trajectory. Need not be implemented."""
        pass

    def write_init(self, system):
        """Subclass should use it to open files for writing."""
        pass

    # These methods must be implemented by subclasses
    def read_sample(self, index):
        """It must return the sample (system) with the given index"""
        raise NotImplementedError()

    def write_sample(self, system, step):
        """It must write a sample (system) to disk. Noting to return."""
        raise NotImplementedError()

    # To read/write timestep and block period sublcasses may implement
    # these methods. The default is dt=1 and blockperiod determined dynamically.
    def read_timestep(self):
        return 1.0

    def write_timestep(self, value):
        pass

    def read_blockperiod(self):
        return None

    def write_blockperiod(self, value):
        pass

    @property
    def timestep(self):
        if self._timestep is None:
            self._timestep = self.read_timestep()
        return self._timestep

    @timestep.setter
    def timestep(self, value):
        self.write_timestep(value)
        self._timestep = value

    @property
    def block_period(self):
        if self._block_period is None:
            self._block_period = self.read_blockperiod()
        if self._block_period is None:
            # If period is still None (read_blockperiod is not
            # implemented) we determine it dynamically
            self._block_period = get_period(self.steps)
        return self._block_period

    @block_period.setter
    def block_period(self, value):
        self._block_period = value
        self.write_blockperiod(value)

    # Some additional useful properties

    @property
    def grandcanonical(self):
        # In subclasses, cache it for efficiency, since we might have to discover it
        if self._grandcanonical is None:
            self._grandcanonical = False
        return self._grandcanonical

    @property
    def times(self):
        """All available times."""
        return [s*self.timestep for s in self.steps]

    @property
    def time_total(self):
        """Total simulation time."""
        return self.steps[-1] * self.timestep

    def time_when_msd_is(self, msd_target, sigma=1.0):
        """Estimate the time when the MSD reaches target_msd in units of sigma^2.
        Bounded by the actual maximum time of trajectory tmax.
        """
        # TODO: when using decorators, this will introduce a diamond, let's move it to a generic function in post processing
        raise NotImplementedError('time_when_msd_is is broken')
        # self._unfold()
        # msd_total = numpy.sum((self._pos_unf[-1] - self._pos_unf[0])**2) / self._pos_unf[0].shape[0]
        # return min(1.0, msd_target * sigma**2 / msd_total) * self.steps[-1] * self.timestep

    def timeseries(self, callback, *args, **kwargs):
        """Returns a timeseries of a callback results"""
        for i, s in enumerate(self):
            yield self.steps[i], callback(s, *args, **kwargs)


class SuperTrajectory2(TrajectoryBase):

    """Collection of subtrajectories"""

    # Optimized version

    def __init__(self, files, trajectoryclass, mode='r'):
        """Group list of files into a single trajectory"""
        self.files = files
        if len(self.files) == 0:
            raise ValueError('no files found in %s' % self.files)
        f = os.path.dirname(self.files[0])
        super(SuperTrajectory2, self).__init__(f, mode)
        self.trajectoryclass = trajectoryclass
        self._timestep = 1.0

        # Make sure subtrajectories are sorted by increasing step
        self.files.sort()
        self.steps = []
        for i, f in enumerate(self.files):
            # This is slow, just to get the step index.
            # If we accept not to have the steps list updated at this stage
            # we can optimize this by about 10% on xyz files (16.12.2016)
            with self.trajectoryclass(f) as t:
                self.steps.append(t.steps[0])

    def read_sample(self, sample):
        f = self.files[sample]
        with self.trajectoryclass(f) as t:
            return t[0]

