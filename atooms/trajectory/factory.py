# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Abstract factory class for trajectories.

The callable `TrajectoryFactory` returns the appropriate trajectory
class based on either explicit rules or suffix euristics. The list of
available trajectory classes can be updated at run-time using the
`update` method.

Examples:
--------

- Reading an hdf5-format trajectory via the factory

    Trajectory = TrajectoryFactory()
    t = Trajectory('input.h5')

is equivalent to

    t = TrajectoryHDF5('input.h5')

- Ask for a specific trajectory format

    Trajectory = TrajectoryFactory()
    t = Trajectory('input.xyz', fmt='rumd')

This is equivalent to

    t = TrajectoryRUMD('input.xyz')
"""

import os
import sys
import re
import inspect


# Note: trajectories should implement a method to check if a file is
# of their own format or not, to avoid relying on suffix check out
# http://stackoverflow.com/questions/456672/class-factory-in-python.
# However, one would have to deal with conflicts anyway, so the
# settings in update() are fine imo.

class TrajectoryFactory(object):

    def __init__(self):
        self.formats = {}
        self.suffixes = {}
        self.callbacks = []

    def _add(self, trj_class, trj_name, overwrite=True):
        """Atomically add one trajectory class to the factory"""
        # We extract the name of the trajectory and lowercase it
        # We expect names of the format [Tag1]Trajectory<Tag2>, where
        # <Tag2> is obligatory and Tag1 is optional. The trajectory key
        # is then Tag1+Tag2 lowercased.
        s = re.search(r'([a-zA-Z0-9]*)Trajectory([a-zA-Z0-9]*)', trj_name)
        if len(s.group(2)) == 0:
            return
        fmt = (s.group(1) + s.group(2)).lower()
        # Always add the new class to the dict of trajectory formats
        # but only add its suffix if not already present (unless we overwrite)
        self.formats[fmt] = trj_class
        if trj_class.suffix not in self.suffixes or overwrite:
            self.suffixes[trj_class.suffix] = trj_class

    def add(self, trajectory_class):
        """Add a trajectory class to the factory"""
        trajectory_name = str(trajectory_class)
        self._add(trajectory_class, trajectory_name)

    def update(self, module, overwrite=True):
        """Update factory parsing all trajectory classes from module"""
        classes = inspect.getmembers(sys.modules[module], inspect.isclass)
        trajectories = [c for c in classes if 'Trajectory' in c[0]]
        for trj_name, trj_class in trajectories:
            if trj_name == 'TrajectoryFactory':
                continue
            self._add(trj_class, trj_name, overwrite)

    def register_callback(self, cbk, *args, **kwargs):
        """
        Register a callback `cbk` to be applied when reading a frame.

        The callback is registered in the Trajectory instance when
        calling the factory, see __call__().
        """
        if (cbk, args, kwargs) not in self.callbacks:
            self.callbacks.append((cbk, args, kwargs))

    def __call__(self, filename, mode='r', fmt=None):
        if fmt is not None:
            return self.formats[fmt](filename, mode)

        th = None
        suffix = os.path.splitext(filename)[-1].replace('.', '')
        if suffix in self.suffixes:
            th = self.suffixes[suffix](filename, mode)
        else:
            # Fallback to hdf5
            if 'hdf5' in self.formats:
                try:
                    th = self.formats['hdf5'](filename, mode)
                except OSError:
                    pass
        if th is not None:
            # Register the callbacks in the actual trajectory instance
            for cbk, args, kwargs in self.callbacks:
                th.register_callback(cbk, *args, **kwargs)
            return th
        else:
            raise ValueError('unknown file format for %s' % filename)
