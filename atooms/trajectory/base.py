# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

import os
import copy

from .utils import get_block_size


# Fields thesaurus: common synonims for fields, such as position -> pos
# It can be used to match fields read by different trajectory classes.
# The first element of each entry is the official one.
FIELDS_DICTIONARY = [('position', 'pos'),
                     ('velocity', 'vel'),
                     ('species', 'spe', 'id')]

def canonicalize_fields(fields):
    for i, field in enumerate(fields):
        for entry in FIELDS_DICTIONARY:
            if field in entry:
                fields[i] = entry[0]
                break
    return fields


class TrajectoryBase(object):
    """
    Trajectory abstract base class.

    A trajectory is composed by one or several frames, each frame
    being a sample of a `System` taken at a given `step` during a
    simulation. `Trajectory` instances are iterable and can be opened
    and closed using the `with` syntax..

        #!python
        with Trajectory(inpfile) as th:
            for system in th:
                pass

    To be fully functional, concrete classes must implement
    `read_sample()` and `write_sample()` methods.

    `read()` is a template composed of the two following steps:

    - `read_init()`: called only once to initialize internal data
    structures, grab metadata, etc. Need *not* be implemented by
    subclasses.

    - `read_sample(n)`: actually return the system at frame n. It must
      be implemented by subclasses.

    Similarly, `write()` is a template composed of `write_init()` and
    `write_sample()`. Only the latter method must be implemented by
    subclasses.

    The `cache` variable reduces acess time when reading the same
    trajectory multiple times. We use shallow copies to cut down the
    overhead. Cache is disabled by default as there is no control on
    its size yet.
    """

    suffix = None

    # TODO: add class callbacks

    def __init__(self, filename, mode='r', cache=False):
        """
        The `mode` can be 'r' (read) or 'w' (write).
        """
        self.filename = filename
        self.mode = mode
        self.callbacks = []
        self.fields = []
        """
        A list of strings describing data to be written by `write_sample`
        and/or read by `read_sample`. Subclasses may use it to filter
        out some data from their format. They can ignore it entirely.
        """
        self.precision = 6
        self.metadata = {}
        """
        Dictionary of metadata about the trajectory. It can be used by
        subclasses to hold trajectory format info or even dynamically
        on a per sample basis,
        """
        # These are cached properties
        self._steps = None
        self._timestep = None
        self._grandcanonical = None
        self._block_size = None
        # Internal state
        self._initialized_write = False
        self._initialized_read = False
        # Sanity checks
        if self.mode == 'r' and not os.path.exists(self.filename):
            raise IOError('trajectory file %s does not exist' % self.filename)
        # Cache frames to optimize reading the same trajectory multiple times
        # We use shallow copies to cut down the overhead
        self.cache = cache
        self._cache = None

    # Trajectory is iterable and supports with syntax

    def __len__(self):
        # We try first with read_len() which returns None by default
        frames = self.read_len()
        if frames is None:
            # We get the steps, which might take a bit longer
            return len(self.steps)
        else:
            return frames

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __iter__(self):
        for i in range(len(self.steps)):
            yield self.read(i)

    def __getitem__(self, key):
        if isinstance(key, slice):
            # This works but it loads the whole trajectory in ram.
            # If we make this a generator, then access a single sample
            # wont work.
            # The Sliced decorator doesn't have this issue.
            # Or just slice range(self.steps) instead.
            frames = range(len(self.steps))
            return [self.read(i) for i in frames[key]]

        elif isinstance(key, int):
            if key < 0:
                key += len(self)
            if key >= len(self):
                raise IndexError("Index (%d) is out of range (%d)." % (key, len(self)))
            return self.read(key)

        else:
            raise TypeError("Invalid argument type [%s]" % type(key))

    def close(self):
        pass

    def read(self, index):
        """Read and return system at frame `index`."""
        if not self._initialized_read:
            self.read_init()
            self._initialized_read = True

        if self.cache and self._cache and index in self._cache:
            # We get the system from the cache
            system = self._cache[index]
        else:
            system = self.read_sample(index)
            if self.cache:
                # Store the system in cache
                if self._cache is None:
                    self._cache = {}
                self._cache[index] = copy.copy(system)

        # TODO: add some means to access the current frame / step in a callback? 11.09.2017
        for cbk, args, kwargs in self.callbacks:
            system = cbk(system, *args, **kwargs)
        return system

    def write(self, system, step):
        """Write `system` at given `step`."""
        if self.mode == 'r':
            raise IOError('trajectory file not open for writing')
        if not self._initialized_write:
            self.write_init(system)
            self._initialized_write = True
        self.write_sample(system, step)
        # Step is added last, frame index starts from 0 by default
        # If step is already there we overwrite (do not append)
        # TODO: just check last step
        if step not in self.steps:
            self.steps.append(step)

    def read_init(self):
        """
        Read metadata and/or set up data structures. Need not be
        implemented.
        """
        pass

    def write_init(self, system):
        """Subclass should use it to open files for writing."""
        pass

    # These methods must be implemented by subclasses

    def read_sample(self, index):
        """Return the system at the given frame `index`."""
        raise NotImplementedError()

    def write_sample(self, system, step):
        """Write a `system` to file. Noting to return."""
        raise NotImplementedError()

    # Callbacks will be applied to the output of read_sample()

    def register_callback(self, cbk, *args, **kwargs):
        if cbk not in self.callbacks:
            self.callbacks.append([cbk, args, kwargs])

    def add_callback(self, cbk, *args, **kwargs):
        """Same as register_callback."""
        self.register_callback(cbk, *args, **kwargs)

    # To read/write timestep and block size sublcasses may implement
    # these methods. The default is dt=1 and block determined dynamically.

    def read_len(self):
        """Return the number of frames. Optional."""
        return None

    def read_steps(self):
        """Return a list of steps."""
        return []

    def read_timestep(self):
        return 1.0

    def write_timestep(self, value):
        self._timestep = value

    def read_block_size(self):
        return None

    def write_block_size(self, value):
        pass

    @property
    def steps(self):
        if self._steps is None:
            if 'r' in self.mode:
                self._steps = self.read_steps()
            else:
                self._steps = []
        return self._steps

    @steps.setter
    def steps(self, value):
        self._steps = value

    @property
    def timestep(self):
        if self._timestep is None:
            if self.mode == 'r':
                self._timestep = self.read_timestep()
            else:
                self._timestep = 1.0
        return self._timestep

    @timestep.setter
    def timestep(self, value):
        self.write_timestep(value)
        self._timestep = value

    @property
    def block_size(self):
        if self._block_size is None:
            self._block_size = self.read_block_size()
        if self._block_size is None:
            # If size is still None (read_block_size is not
            # implemented) we determine it dynamically
            self._block_size = get_block_size(self.steps)
        return self._block_size

    @block_size.setter
    def block_size(self, value):
        self._block_size = value
        self.write_block_size(value)

    # Some additional useful properties

    @property
    def grandcanonical(self):
        """
        True if the trajectory is grandcanonical, i.e. the number of
        particles changes.
        """
        # In subclasses, cache it for efficiency, since we might have to discover it
        if self._grandcanonical is None:
            self._grandcanonical = False
        return self._grandcanonical

    @property
    def times(self):
        """All available times."""
        return [s*self.timestep for s in self.steps]

    @property
    def total_time(self):
        """Total simulation time."""
        return (self.steps[-1] - self.steps[0]) * self.timestep


class SuperTrajectory(TrajectoryBase):

    """Collection of subtrajectories."""

    # Optimized version

    def __init__(self, files, trajectoryclass, mode='r'):
        """
        Group a list of `files` into a single trajectory of class
        `trajectoryclass`.
        """
        self.files = files
        if len(self.files) == 0:
            raise ValueError('no files found in %s' % self.files)
        f = os.path.dirname(self.files[0])
        super(SuperTrajectory, self).__init__(f, mode)
        self.trajectoryclass = trajectoryclass

        # Make sure subtrajectories are sorted by increasing step
        self.files.sort()
        # This list holds the file containing a given step
        self._steps_file = []
        self._steps_frame = []
        # This caches the last trajectory used to minimize __init__() overhead
        self._last_trajectory = None
        self.steps = []
        for f in self.files:
            # This is slow, just to get the step index.
            # If we accept not to have the steps list updated at this stage
            # we can optimize this by about 10% on xyz files (16.12.2016)
            with self.trajectoryclass(f) as t:
                for j, step in enumerate(t.steps):
                    if len(self.steps) == 0 or step != self.steps[-1]:
                        self.steps.append(step)
                        self._steps_file.append(f)
                        self._steps_frame.append(j)

    def read_sample(self, frame):
        f = self._steps_file[frame]
        j = self._steps_frame[frame]
        # Optimization: use the last trajectory in cache (it works
        # well if frames are read sequentially)
        if self._last_trajectory is None:
            self._last_trajectory = self.trajectoryclass(f)
        elif self._last_trajectory.filename != f or \
             self._last_trajectory.trajectory.closed:
            # Careful: we must check if the file object has not been closed in the meantime.
            # This can happen with class decorators.
            self._last_trajectory.close()
            self._last_trajectory = self.trajectoryclass(f)
        else:
            # In cache
            pass
        t = self._last_trajectory
        return t[j]

    def read_timestep(self):
        with self.trajectoryclass(self.files[0]) as t:
            return t.timestep

    def close(self):
        if self._last_trajectory is not None:
            self._last_trajectory.close()
            self._last_trajectory = None
        super(SuperTrajectory, self).close()
