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
    simulation. Trajectory instances are iterable and behave as file
    objects: they should be opened and closed using the `with` syntax

        #!python
        with Trajectory(inpfile) as th:
            for system in th:
                pass

    alternatively

        #!python
        th = Trajectory(inpfile)
        for system in th:
            pass
        th.close()

    To write the state of a `System` to a trajectory, we must open the
    trajectory in write mode.

        #!python
        with Trajectory(outfile, 'w') as th:
            th.write(system, step=0)

    Trajectories can use the `fields` variable to define which system
    properties will be written to disk. Concrete classes may thus
    provide means to write arbitrary fields to disk via particle
    properties. Example:

        #!python
        for particle in system:
            particle.some_custom_property = 1.0
        with Trajectory(outfile, 'w') as th:
            th.fields = ['position', 'some_custom_property']
            th.write(system, step=0)

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

    The `cache` variable reduces access time when reading the same
    trajectory multiple times. We use shallow copies to cut down the
    overhead. Cache is disabled by default as there is no control on
    its size yet.
    """

    suffix = None

    # Mutable class variables are shared by subclasses https://stackoverflow.com/a/13404513

    # To avoid that subclasses modify the parent class variable, we
    # use the trick to define class_callbacks as None, then assign it
    # to a list in the class method
    class_callbacks = None

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
        self._overwrite = False
        # These are cached properties
        self._steps = None
        self._timestep = None
        self._grandcanonical = None
        self._block_size = None
        # Internal state
        self._initialized_write = False
        self._initialized_read = False
        # Sanity checks
        # Subclasses may define mode variants, such as r+, therefore the
        # autoritative mode is the first letter
        if self.mode[0] == 'r' and self.filename is not None and \
           not os.path.exists(self.filename):
            raise IOError('trajectory file %s does not exist' % self.filename)
        # Cache frames to optimize reading the same trajectory multiple times
        # We use shallow copies to cut down the overhead
        self.cache = cache
        self._cache = None

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

    def append(self, system):
        self.write(system)

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

        # Apply callbacks, class level and instance level
        #
        # Pass the frame index to the callback via system
        # by storing a "temporary" frame instance variable
        # (deleted after this loop)
        system.frame = index
        callbacks = copy.copy(self.callbacks)
        if self.class_callbacks is not None:
            callbacks += self.class_callbacks
        for cbk, args, kwargs in callbacks:
            system = cbk(system, *args, **kwargs)
        if hasattr(system, 'frame'):
            del(system.frame)
        return system

    def write(self, system, step=None):
        """
        Write `system` at a given integer `step`.

        If `step` is not provided, it is defined internally by
        incrementing by one the last added step, staring from zero.
        """
        if self.mode[0] == 'r':
            raise IOError('trajectory file not open for writing')
        if not self._initialized_write:
            self.write_init(system)
            self._initialized_write = True

        # If we do not provide a step, we incrementally add 1 to the
        # last step, starting from 0.
        if step is not None:
            current_step = step
        else:
            if len(self.steps) == 0:
                current_step = 0
            else:
                current_step = self.steps[-1] + 1

        # If overwriting is not allowed (default), we check that we
        # are adding a step larger than the last added step.
        if not self._overwrite:
            if len(self.steps) > 0 and current_step <= self.steps[-1]:
                raise ValueError('cannot add step {} when overwrite is False'.format(current_step))

        # Write the sample.
        self.write_sample(system, current_step)
        # Step is added last, frame index starts from 0 by default.
        if step is None:
            self.steps.append(current_step)
        else:
            # If overwriting is allowed, we append the step only if it
            # not already there. This enables a small optimization by
            # avoiding this check if overwriting is off.
            if not self._overwrite:
                self.steps.append(current_step)
            else:
                if current_step not in self.steps:
                    self.steps.append(current_step)

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
        """
        Register a callback `cbk` to be applied when reading a frame.

        The callback must receive a `System` instance as input an
        return the modified `System` instance as output.
        """
        if (cbk, args, kwargs) not in self.callbacks:
            self.callbacks.append((cbk, args, kwargs))

    def add_callback(self, cbk, *args, **kwargs):
        """Same as `register_callback()`."""
        self.register_callback(cbk, *args, **kwargs)

    @classmethod
    def register_class_callback(cls, cbk, *args, **kwargs):
        """
        Register a class level callback `cbk` to be applied when reading a frame.

        The callback must receive a `System` instance as input an
        return the modified `System` instance as output.

        Class level callbacks are applied to all instances of a class.
        """
        if cls.class_callbacks is None:
            cls.class_callbacks = []
        if (cbk, args, kwargs) not in cls.class_callbacks:
            cls.class_callbacks.append((cbk, args, kwargs))

    @classmethod
    def add_class_callback(self, cbk, *args, **kwargs):
        """Same as `register_class_callback()`."""
        self.register_class_callback(cbk, *args, **kwargs)

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
            if self.mode[0] == 'r':
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
            if self.mode[0] == 'r':
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
