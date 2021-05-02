import tempfile
import json

import numpy
import warnings

from atooms.core.utils import rmf
from atooms.system import System as _System
from atooms.system import Particle, Cell
from .interaction import Interaction
from .verlet_list import VerletList
from .helpers import _merge_source, _normalize_path

import f2py_jit


__all__ = ['Interaction', 'System', 'NeighborList', 'VerletList',
           'Trajectory', 'Particle', 'Cell']


class NeighborList(object):

    def __init__(self, rcut, neighbors='neighbor_list.f90',
                 helpers='helpers.f90', inline=True):
        self.rcut = numpy.array(rcut)
        self.neighbors = None
        self.number_neighbors = None
        self._module_path = None

        # Gather f90 sources into a single one
        source = _merge_source(helpers, neighbors)

        # Inline subroutine
        if inline:
            from f2py_jit.finline import inline_source
            source = inline_source(source, ignore='compute')  # avoid reinlining forces!

        # Compile and bundle the module with f2py
        # Build a unique module.
        # Every model with its own parameter combination corresponds to a unique module
        # and can be safely reused (up to changes in interaction / helpers)
        extra_args = '--opt="-O3 -ffast-math"'
        uid = f2py_jit.build_module(source,
                                    metadata={"neighbors": neighbors,
                                              "helpers": helpers},
                                    extra_args=extra_args,
                                    db_file='.atooms_jit.json')

        # Store module name (better not store the module itself, else we cannot deepcopy)
        self._uid = uid

    def _setup(self, npart, nneigh):
        """Allocate or reallocate arrays for neighbor list"""
        if self.neighbors is None or self.neighbors.shape[1] != npart or self.neighbors.shape[0] < nneigh:
            self.neighbors = numpy.ndarray(shape=(nneigh, npart), order='F', dtype=numpy.int32)
        if self.number_neighbors is None or len(self.number_neighbors) != npart:
            self.number_neighbors = numpy.ndarray(npart, order='F', dtype=numpy.int32)

    def compute(self, box, pos, ids):
        # Setup
        f90 = f2py_jit.import_module(self._uid)

        # Setup arrays
        # Estimate max number of neighbors based on average density
        # We take the largest cut off distance
        npart = pos.shape[1]
        rho = npart / box.prod()
        nneigh = int(4.0 / 3.0 * 3.1415 * rho * numpy.max(self.rcut)**3 * 1.50)
        self._setup(npart, nneigh)
        # Compute neighbors list
        #
        # If the f90 code returns an error, the arrays are reallocated
        # based on the largest number of neighbors returned by the f90
        # routine
        error = f90.neighbor_list.compute(box, pos, ids, self.rcut, self.neighbors, self.number_neighbors)
        if error:
            self._setup(npart, max(self.number_neighbors))
            error = f90.neighbor_list.compute(box, pos, ids, self.rcut, self.neighbors, self.number_neighbors)
            assert not error, "something wrong with neighbor_list"


class System(_System):

    def compute_interaction(self, what='forces'):
        if not hasattr(self, 'interaction'):
            raise ValueError('system has no interaction')
        box = self.dump('cell.side', view=True)
        ids = self.dump('particle.species', view=True, dtype='int32')
        pos = self.dump('particle.position', order='F', view=True)
        self.interaction.compute(what, box, pos, ids)

    def _compute_neighbor_list(self):
        assert hasattr(self, 'neighbor_list'), 'system has no neighbor_list'
        box = self.dump('cell.side', view=True)
        ids = self.dump('particle.species', view=True, dtype='int32')
        pos = self.dump('particle.position', order='F', view=True)
        self.neighbor_list.compute(box, pos, ids)
        # Assign neighbors to particles
        nl = self.neighbor_list
        for i, p in enumerate(self.particle):
            p.neighbors = nl.neighbors[0: nl.number_neighbors[i], i]


def _wrap_system(system):
    from atooms.trajectory.decorators import change_species
    new_system = System()
    new_system.update(system)
    new_system = change_species(new_system, 'F')
    return new_system


def _add_interaction(trajectory, system):
    """
    Add interaction to trajectory from model metadata information or
    accompanying json file.
    """
    try:
        import atooms.models
    except ImportError:
        return system

    try:
        # Lookup in atooms models database
        # TODO: accept serialized json potential/cutoff in metadata
        name = trajectory.metadata['model']
        db = atooms.models.load()
        model = db[name]
    except KeyError:
        # If that fails, look for an accompanying json file
        try:
            model = atooms.models.read_json(trajectory.filename + '.json')
        except IOError:
            return system
            #raise ValueError('cannot set interaction')

    system.interaction = Interaction(model,
                                     interaction='interaction.f90',
                                     helpers='helpers_3d.f90',
                                     inline=True)
    return system

# Boost all trajectory classes with specific jit callbacks
# and create a new Trajectory factory that only loads these classes
# They keep their original names but live in this module namespace
# The original classes are untouched.
from atooms.trajectory import Trajectory as __Trajectory
from atooms.trajectory.factory import TrajectoryFactory
Trajectory = TrajectoryFactory()
for key in __Trajectory.formats:
    new_name = __Trajectory.formats[key].__name__
    old_cls = __Trajectory.formats[key]
    cls = type(new_name, (old_cls, ), dict(old_cls.__dict__))
    cls.add_class_callback(_wrap_system)
    cls.add_self_callback(_add_interaction)
    Trajectory.add(cls)

# Lock the xyz format
Trajectory.suffixes['xyz'] = Trajectory.formats['xyz']
