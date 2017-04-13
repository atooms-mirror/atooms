# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Read and write trajectory files in various formats."""

import os
import sys
import logging

from atooms.utils import NullHandler
logging.getLogger(__name__).addHandler(NullHandler())

from .utils import convert, split
from .base import SuperTrajectory
from .decorators import *

# Import all trajectories
from .xyz import TrajectorySimpleXYZ, TrajectoryXYZ, TrajectoryNeighbors
from .pdb import TrajectoryPDB
from .hoomd import TrajectoryHOOMD
from .rumd import TrajectoryRUMD, SuperTrajectoryRUMD
from .lammps import TrajectoryLAMMPS
try:
    from .hdf5 import TrajectoryHDF5
except ImportError:
    pass
from .ram import TrajectoryRam
from .field import TrajectoryField

# We build an instance of trajectory factory here and update the
# trajectory classes. Client code can update the factory with new
# classes at run time by calling update(__name__) again.
#
# The list of available formats is stored as a dict in
# Trajectory.formats, the suffixes as Trajectory.suffixes.

from .factory import TrajectoryFactory
Trajectory = TrajectoryFactory()
"""An instance of a Trajectory factory, see `TrajectoryFactory`."""
Trajectory.update(__name__)

# Update factory with plugins modules
import pkgutil
import atooms.plugins
for _, mod_name, _ in pkgutil.iter_modules(atooms.plugins.__path__, prefix='atooms.plugins.'):
    m = __import__(mod_name)
    Trajectory.update(mod_name)

# Additional plugins can be put in the atooms_plugins module
try:
    import atooms_plugins
    for _, mod_name, _ in pkgutil.iter_modules(atooms_plugins.__path__, prefix='atooms_plugins.'):
        m = __import__(mod_name)
        Trajectory.update(mod_name)
except ImportError:
    pass
