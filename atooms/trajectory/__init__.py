# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import sys
import inspect
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())

from .utils import convert, split
from .base import SuperTrajectory
from .decorators import *

# Import all trajectories
from xyz import TrajectorySimpleXYZ, TrajectoryXYZ, TrajectoryNeighbors
from pdb import TrajectoryPDB
from hoomd import TrajectoryHOOMD
from rumd import TrajectoryRUMD, SuperTrajectoryRUMD
from lammps import TrajectoryLAMMPS
try:
    from hdf5 import TrajectoryHDF5
except:
    pass

# Load third party trajectories (if any).
# Perhaps use some dynamic import http://stackoverflow.com/questions/301134/dynamic-module-import-in-python ?
# try:
#     from atooms.trajectory.plugins import *
# except:
#     pass

def factory(extra_module=None):
    classes = inspect.getmembers(sys.modules[__name__], inspect.isclass)
    trajectories = [c for c in classes if c[0].startswith('Trajectory')] 
    if extra_module is not None:
        classes = inspect.getmembers(sys.modules[extra_module], inspect.isclass)
        trajectories += [c for c in classes if c[0].startswith('Trajectory')] 
    __factory_suf = {}
    __factory_fmt = {}
    for trj_name, trj_class in trajectories:
        # We extract the name of the trajectory and lowercase it
        fmt = trj_name[len('Trajectory'):].lower()
        try:
            __factory_suf[trj_class.suffix] = trj_class
            __factory_fmt[fmt] = trj_class
        except AttributeError:
            pass

    # Lock some common suffixes to specific formats
    __factory_suf['xyz'] = TrajectoryXYZ
    return __factory_fmt, __factory_suf

available_formats = factory()[0]

# # List all imported trajectory classes, store their names and suffixes
# # to setup the factory.
# classes = inspect.getmembers(sys.modules[__name__], inspect.isclass)
# trajectories = [c for c in classes if c[0].startswith('Trajectory')] 

# # Factory method which mimics an abstract factory class.
# __factory_suf = {}
# __factory_fmt = {}
# for trj_name, trj_class in trajectories:
#     # We extract the name of the trajectory and lowercase it
#     fmt = trj_name[len('Trajectory'):].lower()
#     try:
#         __factory_suf[trj_class.suffix] = trj_class
#         __factory_fmt[fmt] = trj_class
#     except AttributeError:
#         pass

# # This variable can be used to inspect available trajectory classes
# available_formats = __factory_fmt

# # Lock some common suffixes to specific formats
# __factory_suf['xyz'] = TrajectoryXYZ

# Trajectories should implement a method to check if a file is of
# their own format or not, to avoid relying on suffix check out
# http://stackoverflow.com/questions/456672/class-factory-in-python.
# However, one would have to deal with conflicts manually, so the
# setting above is fine imo.

def Trajectory(filename, mode='r', fmt=None):
    """Factory class shortcut."""
    # If we are passed an explicit format, we use that
    if fmt is not None:
        return __factory_fmt[fmt](filename, mode)
    
    suffix = os.path.splitext(filename)[-1].replace('.', '')
    if suffix in __factory_suf:
        return __factory_suf[suffix](filename, mode)
    else:
        raise ValueError('unknown file suffix %s' % suffix)
