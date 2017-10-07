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
from .xyz import TrajectoryXYZ
try:
    from .hdf5 import TrajectoryHDF5
except ImportError:
    pass


# Note: trajectories should implement a method to check if a file is
# of their own format or not, to avoid relying on suffix check out
# http://stackoverflow.com/questions/456672/class-factory-in-python.
# However, one would have to deal with conflicts anyway, so the
# settings in update() are fine imo.

class TrajectoryFactory(object):

    def __init__(self):
        self.formats = {}
        self.suffixes = {}

    def update(self, module):
        classes = inspect.getmembers(sys.modules[module], inspect.isclass)
        trajectories = [c for c in classes if 'Trajectory' in c[0]]
        for trj_name, trj_class in trajectories:
            # We extract the name of the trajectory and lowercase it
            # We expect names of the format [Tag1]Trajectory<Tag2>, where
            # <Tag2> is obligatory and Tag1 is optional. The trajectory key
            # is then Tag1+Tag2 lowercased.
            s = re.search(r'([a-zA-Z0-9]*)Trajectory([a-zA-Z0-9]*)', trj_name)
            if len(s.group(2)) == 0:
                continue
            fmt = (s.group(1) + s.group(2)).lower()
            try:
                self.suffixes[trj_class.suffix] = trj_class
                self.formats[fmt] = trj_class
            except AttributeError:
                pass
        # Lock some common suffixes to specific formats
        self.suffixes['xyz'] = TrajectoryXYZ

    def __call__(self, filename, mode='r', fmt=None):
        if fmt is not None:
            return self.formats[fmt](filename, mode)

        suffix = os.path.splitext(filename)[-1].replace('.', '')
        if suffix in self.suffixes:
            return self.suffixes[suffix](filename, mode)
        else:
            # Fallback to hdf5
            try:
                return TrajectoryHDF5(filename, mode)
            except ImportError:
                raise ValueError('hdf5 library not installed')
            except:
                raise ValueError('unknown file format for %s' % filename)
