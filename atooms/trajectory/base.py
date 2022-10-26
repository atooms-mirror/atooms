# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Base trajectory classes.
"""

import os
import copy
import warnings


from .utils import get_block_size


# FutureWarning must show a traceback so that the user can localize the deprecation
# This should be moved up the chain
def _warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    import traceback
    import sys
    if category is FutureWarning:
        traceback.print_stack(file=sys.stderr, limit=6)
    sys.stderr.write(warnings.formatwarning(message, category, filename, lineno, line))


warnings.showwarning = _warn_with_traceback


def canonicalize_fields(fields):
    from atooms.core.utils import canonicalize

    warnings.warn('canonicalize_fields() is deprecated, use canonicalize() instead', FutureWarning)
    th = TrajectoryBase(None)
    return canonicalize(fields, th.thesaurus)


class TrajectoryBase(object):
    """
    Trajectory abstract base class.

    A trajectory is composed by one or several frames, each frame
    being a sample of a `System` taken at a given `step` during a
    simulation. Trajectory instances are iterable and behave as file
    objects: they can be opened and closed using the `with` syntax

    ```python
    with Trajectory(inpfile) as th:
        for system in th:
            pass
    ```

    alternatively

    ```python
    th = Trajectory(inpfile)
    for system in th:
        pass
    th.close()
    ```

    To write the state of a `System` to a trajectory, we must open the
    trajectory in write mode

    ```python
    with Trajectory(outfile, 'w') as th:
        th.write(system, step=0)
    ```

    Trajectories can use the `variables` attribute to define which system
    properties will be written to disk. Concrete classes may thus
    provide means to write arbitrary variables to disk via particle
    properties. Example:

    ```python
    for particle in system:
        particle.some_custom_property = 1.0
    with Trajectory(outfile, 'w') as th:
        th.variables = ['position', 'some_custom_property']
        th.write(system, step=0)
    ```

    To be fully functional, concrete classes must implement
    `read_system()` and `write_system()` methods.

    `read()` is a template composed of the two following steps:

    - `read_init()`: called only once to initialize internal data
    structures, grab metadata, etc. Need *not* be implemented by
    subclasses.

    - `read_system(n)`: actually return the system at frame n. It must
      be implemented by subclasses.

    Similarly, `write()` is a template composed of `write_init()` and
    `write_system()`. Only the latter method must be implemented by
    subclasses.

    The `cache` variable reduces access time when reading the same
    trajectory multiple times. We use shallow copies to cut down the
    overhead. Cache is disabled by default as there is no control on
    its size yet.
    """

    suffix = None

    def __init__(self, filename, mode='r', cache=False):
        self.filename = filename
        self.mode = mode
        """Can be 'r' (read) or 'w' (write)."""
        self.callbacks = []
        """List of callbacks to be applied to the system when reading a frame, see `add_callback()`"""
        self.precision = 6
        """Number of digits used to write trajectory files."""
        self.metadata = {}
        """
        Dictionary of metadata about the trajectory. It can be used by
        subclasses to hold trajectory format info or dynamically
        on a per sample basis.
        """
        self.thesaurus = {
            'position': 'particle.position',
            'velocity': 'particle.velocity',
            'species': 'particle.species',
            'radius': 'particle.radius',
            'mass': 'particle.mass',
            'pos': 'particle.position',
            'x': 'particle.position[0]',
            'y': 'particle.position[1]',
            'z': 'particle.position[2]',
            'vel': 'particle.velocity',
            'vx': 'particle.velocity[0]',
            'vy': 'particle.velocity[1]',
            'vz': 'particle.velocity[2]',
            'id': 'particle.species',
            'type': 'particle.species',
            'name': 'particle.species',
        }
        """
        Dictionary of common shortcuts and synonims for system
        attributes. Can be updated by subclasses.
        """

        # These are cached properties
        self._variables = ()
        """
        Tuple of system attributes to be written by `write_system` and/or
        read by `read_system`. Its entries are canonicalized using
        `self.thesaurus` everytime the attribute is set. Subclasses
        may use it to allow the user to modify the trajectory layout
        or they can ignore it entirely.
        """
        self._steps = None
        self._timestep = None
        self._grandcanonical = None
        self._block_size = None

        # Internal state
        self._overwrite = False
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
        """If `True`, use a cache when reading the same frame multiple times"""
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
        """Equivalent to `write(system)`."""
        self.write(system)

    def close(self):
        """Can be implemented by subclasses to close file handles."""
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
            system = self.read_system(index)
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
        # Copying callbacks is an old hack that allowed client
        # code to safely pass variables in the callback across calls
        # Deprecating this feature would require checking if the
        # callback has variables not starting with __
        callbacks = copy.copy(self.callbacks)
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
        self.write_system(system, current_step)
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

    def read_system(self, frame):
        """Return the system at the given `frame`."""
        if hasattr(self, 'read_sample'):
            warnings.warn('read_sample() is deprecated, use read_system()', DeprecationWarning)
            return self.read_sample(frame)
        else:
            raise NotImplementedError()

    def write_system(self, system, step):
        """Write a `system` to file. Noting to return."""
        if hasattr(self, 'write_sample'):
            warnings.warn('write_sample() is deprecated, use write_system()', DeprecationWarning)
            return self.write_sample(system, step)
        else:
            raise NotImplementedError()

    # Callbacks will be applied to the output of read_system()

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

    # To read/write timestep and block size sublcasses may implement
    # these methods. The default is dt=1 and block determined dynamically.

    def read_len(self):
        """Return the number of frames. Optional."""
        return None

    def read_steps(self):
        """Return a list of steps."""
        return []

    def read_timestep(self):
        """Can be subclassed to parse the timestep."""
        return 1.0

    def write_timestep(self, value):
        """Set the timestep."""
        self._timestep = value

    def read_block_size(self):
        """
        Can be subclassed to parse the block size for trajectory with
        non-uniform intervals between frames (ex. exponential sampling).
        """
        return None

    def write_block_size(self, value):
        """Can be subclassed to write the block size."""
        pass

    # Properties

    @property
    def variables(self):
        return self._variables

    @variables.setter
    def variables(self, value):
        from atooms.core.utils import canonicalize
        self._variables = canonicalize(value, self.thesaurus)

    def copy(self, cls=None, fout=None, only=None, include=None, exclude=None, steps=None):
        """
        Return a copy of the trajectory
        """
        from atooms.core.progress import progress
        from atooms.core.utils import canonicalize
        from atooms.core.utils import mkdir

        # Output class
        if cls is None:
            cls = self.__class__
        else:
            from atooms.trajectory import Trajectory
            if isinstance(cls, str):
                cls = Trajectory.formats[cls]

        # Make sure parent folder exists
        if fout is not None:
            mkdir(os.path.dirname(fout))

        # Copy trajectory
        conv = cls(fout, 'w')
        # In a previous version of trajectory conversion, we only
        # included the variables of the original trajectory, to allow
        # the output trajectory to include variables not present in
        # the original one. I do not think this makes sense. Suppose
        # the original trajectory did not store velocities, but the
        # output one writes them by default. The results will be wrong.
        variables = list(self.variables)
        if only is not None:
            variables = only
        if include is not None:
            for pattern in canonicalize(include, self.thesaurus):
                if pattern not in variables:
                    variables.append(pattern)
        if exclude is not None:
            for pattern in canonicalize(exclude, self.thesaurus):
                if pattern not in variables:
                    variables.append(pattern)
        conv.variables = variables
        conv.precision = self.precision
        conv.timestep = self.timestep
        conv.block_size = self.block_size
        # In python 3, zip returns a generator so this is ok
        #
        # for system, step in zip(inp, inp.steps):
        #     conv.write(system, step)
        #
        # In python 2, zipping t and t.steps will load everything
        # in RAM. In this case, it is better to use enumerate()
        if steps is None:
            for i, system in progress(enumerate(self), total=len(self)):
                conv.write(system, self.steps[i])
        else:
            # Only include requested steps (useful to prune
            # non-periodic trajectories)
            for step in steps:
                idx = self.steps.index(step)
                conv.write(self[idx], step)
        return conv

    @property
    def fields(self):
        warnings.warn('fields is deprecated, use variables instead', FutureWarning)
        return self.variables

    @fields.setter
    def fields(self, value):
        warnings.warn('fields is deprecated, use variables instead', FutureWarning)
        self._variables = value

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
        """
        Return the block size in frames of the trajectory.

        A block is a sequence of steps that repeats periodically in
        the trajectory. The block size is 1 for trajectories with
        constant time intervals between frames. A block size larger
        than 1 occurs with exponential sampling, ex. for steps
        sequences like 1,2,4,11,12,14,...
        """
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
        # We look up the first step in each subtrajectory
        _step = []
        for f in self.files:
            with self.trajectoryclass(f) as t:
                _step.append(t.steps[0])
        _, self.files = zip(*sorted(zip(_step, self.files)))
        self.files = list(self.files)

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
                    # The second test is to avoid storing twice the same frame
                    # if it appears at the end of a file and the beginning
                    # of the next file
                    if len(self.steps) == 0 or step != self.steps[-1]:
                        self.steps.append(step)
                        self._steps_file.append(f)
                        self._steps_frame.append(j)

        # Propagate variables
        with self.trajectoryclass(self.files[0]) as t:
            self.variables = t.variables

    def read_system(self, frame):
        # Optimization: use the last trajectory in cache (it works
        # well if frames are read sequentially)
        f = self._steps_file[frame]
        j = self._steps_frame[frame]

        # Note: if the file object has been closed in the meantime,
        # this approach will fail. In atooms 2, there was a check
        # whether the file object (trajectory attribute) was
        # closed. We will assume that decorators have no side
        # effects. If this needs to be fixed again, we can find a way.
        if self._last_trajectory is None:
            # First trajectory
            self._last_trajectory = self.trajectoryclass(f)
        elif self._last_trajectory.filename != f:
            # This is a new trajectory file
            self._last_trajectory.close()
            self._last_trajectory = self.trajectoryclass(f)
        else:
            # Already in cache, we just pass for clarity
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
