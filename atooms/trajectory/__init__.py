# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""Read and write trajectory files in various formats."""

import logging
import pkgutil

from .utils import split
from .base import SuperTrajectory
from .decorators import *

# Import all trajectories
from .xyz import TrajectorySimpleXYZ, TrajectoryXYZ, TrajectoryNeighbors
from .exyz import TrajectoryEXYZ
from .pdb import TrajectoryPDB
from .hoomd import TrajectoryHOOMD
from .rumd import TrajectoryRUMD, SuperTrajectoryRUMD
from .lammps import TrajectoryLAMMPS, TrajectoryFolderLAMMPS
try:
    from .hdf5 import TrajectoryHDF5
except ImportError:
    pass
try:
    from .gsd import TrajectoryGSD
except ImportError:
    pass
from .ram import TrajectoryRam

from atooms.core.utils import NullHandler
logging.getLogger(__name__).addHandler(NullHandler())

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
import atooms.plugins
for _, _mod_name, _ in pkgutil.iter_modules(atooms.plugins.__path__, prefix='atooms.plugins.'):
    m = __import__(_mod_name)
    Trajectory.update(_mod_name, overwrite=False)

# Make sure the current directory is in the path. This is needed by
# the plugin mechanism. Not sure why '' should not be there btw.
import sys
__added = False
if '' not in sys.path:
    __added = True
    sys.path.append('')

# Additional plugins can be put in the atooms_plugins module
try:
    import atooms_plugins
except ImportError:
    pass
else:
    for _, _mod_name, _ in pkgutil.iter_modules(atooms_plugins.__path__,
                                                prefix='atooms_plugins.'):
        try:
            m = __import__(_mod_name)
            Trajectory.update(_mod_name, overwrite=False)
        except ImportError as err:
            # Usually it is an error in the plugin module so we show
            print(err)
            # Could not import this trajectory
            pass

# Clean up
if __added:
    sys.path.remove('')
del(__added)
