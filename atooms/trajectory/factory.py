# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Abstract factory class for trajectories.

Construct a callable that returns the appropriate trajectory class
based on explicit or suffix-based rules. The list of trajectory
classes can be updated at run time using update().

Examples: 

- Read an hdf5-format trajectory file

>>> Trajectory = TrajectoryFactory()
>>> t = Trajectory('input.h5')

This is equivalent to 

>>> t = TrajectoryHDF5('input.h5')

Force reading an xyz trajectory file using rumd format

>>> Trajectory = TrajectoryFactory()
>>> t = Trajectory('input.xyz', fmt='rumd')

This is equivalent to

>>> t = TrajectoryRUMD('input.xyz')

"""

import os
import sys
import inspect
from .xyz import TrajectoryXYZ
from .hdf5 import TrajectoryHDF5

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
        trajectories = [c for c in classes if c[0].startswith('Trajectory')] 
        for trj_name, trj_class in trajectories:
            # We extract the name of the trajectory and lowercase it
            fmt = trj_name[len('Trajectory'):].lower()
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
            except:
                raise ValueError('unknown file format for %s' % filename)
