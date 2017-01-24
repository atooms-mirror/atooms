# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import sys
import logging

from .utils import convert, split
from .base import SuperTrajectory
from .decorators import *

logging.getLogger(__name__).addHandler(logging.NullHandler())

# Factory method which mimics an abstract factory class
__factory_map = {}

from xyz import *
from pdb import TrajectoryPDB
__factory_map['xyz'] = TrajectoryXYZ
__factory_map['pdb'] = TrajectoryPDB

from hoomd import TrajectoryHOOMD
__factory_map['tgz'] = TrajectoryHOOMD

from rumd import TrajectoryRUMD, SuperTrajectoryRUMD
from lammps import TrajectoryLAMMPS

try:
    from hdf5 import TrajectoryHDF5
    __factory_map['h5'] = TrajectoryHDF5
    __factory_map['dat'] = TrajectoryHDF5
except:
    pass

# Load plugins (if plugins are found)
try:
    from atooms.plugins.trajectory import *
except:
    # No plugins found
    pass

# TODO: trajectories should implement a method to check if a file
# is of their own format or not, to avoid relying on suffix
# check out http://stackoverflow.com/questions/456672/class-factory-in-python

# def Trajectory(filename, fmt=None, **kwargs):
# Make interface consistent with base, fmt is specified thorugh extension only
# It should be called Trajectory(fname, mode='w', fmt='h5') or
# Trajectory(fname, 'w', fmt='h5')
def Trajectory(filename, mode='r', fmt='' ): #, *args, **kwargs):
    """ Factory class shortcut """
    suffix = os.path.splitext(filename)[-1].replace('.', '')
    if not suffix in __factory_map:
        # always try hdf5 to accomodate non standard suffixes
        # we could even search within the path... :-)
        suffix = 'h5'
        #raise ValueError('unknown file type %s' % suffix)
        
    return __factory_map.get(suffix)(filename, mode)
#    return __factory_map.get(suffix)(filename, **kwargs)
