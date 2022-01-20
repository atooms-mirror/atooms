from __future__ import absolute_import

import warnings
import numpy as np

# With absolute import (default in python 3) this does not clash
import gsd
import gsd.hoomd

from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System
from atooms.trajectory.base import TrajectoryBase


class TrajectoryGSD(TrajectoryBase):

    """
    Glotzer group's binary GSD format for HOOMD (https://glotzerlab.engin.umich.edu/hoomd-blue/)
    """
    suffix = 'gsd'

    def __init__(self, filename, mode='r', fields=None):
        super(TrajectoryGSD, self).__init__(filename, mode)
        self.variables = ['particle.species', 'particle.position']
        if fields is not None:
            self.variables = fields
            warnings.warn('fields is deprecated, use variables instead', FutureWarning)

        # self.mode can be 'w' or 'r', but gsd is a binary format, so it only accepts 'wb' or 'rb'.
        file_mode = self.mode + "b"
        # Trajectory file handle
        self._file = gsd.hoomd.open(name=self.filename, mode=file_mode)
        # When reading, we must define the steps.
        if self.mode == 'r':
            self.steps = [snap.configuration.step for snap in self._file]

    def read_system(self, frame):
        """ returns System instance. """
        snap = self._file[frame]
        ndim = snap.configuration.dimensions

        # Convert typeid from [0, 0, 1, ...] to ['A', 'A', 'B', ...] when snap.particles.types = ['A', 'B']
        distinct_species = snap.particles.types
        distinct_typeids = list(range(len(distinct_species)))
        typeid_to_species = {}
        for i in distinct_typeids:
            typeid_to_species[i] = distinct_species[i]

        box = snap.configuration.box[:ndim]    # atooms does not handle sheared boxes.
        cell = Cell(side=box)

        N = snap.particles.position.shape[0]
        particles = []
        for i in range(N):
            p = Particle(
                mass=snap.particles.mass[i],
                species=typeid_to_species[snap.particles.typeid[i]],
                position=snap.particles.position[i, :ndim],
                velocity=snap.particles.velocity[i, :ndim],
                radius=snap.particles.diameter[i] / 2
            )
            particles.append(p)

        return System(particle=particles, cell=cell)

    def write_system(self, system, step):
        """ Writes to the file handle self.trajectory."""
        variables = self.variables
        data = {what: system.dump(what) for what in ['particle.position', 'particle.velocity',
                                                     'particle.species', 'particle.mass', 'particle.radius']}
        box = system.cell.side
        N = len(system.particle)
        distinct_species = system.distinct_species

        # Convert species from ['A', 'A', 'B', ...] to [0, 0, 1, ...] when distinct_species = ['A', 'B']
        species_to_typeid = {}
        typeid = 0
        for species in distinct_species:
            species_to_typeid[species] = typeid
            typeid += 1

        pos = data['particle.position']
        vel = data['particle.velocity']
        species = data['particle.species']    # This is 'A', 'A', 'B', etc when distinct_species = ['A', 'B']
        mass = data['particle.mass']
        radius = data['particle.radius']
        convert_to_typeid = np.vectorize(lambda species: species_to_typeid[species])
        typeid = convert_to_typeid(species).astype(int)  # This is 0, 0, 1, etc when distinct_species = ['A', 'B']

        snap = gsd.hoomd.Snapshot()
        snap.configuration.box = [box[0], box[1], box[2], 0, 0, 0]  # Assume all strains 0.
        snap.configuration.step = step
        snap.particles.types = distinct_species
        snap.particles.N = N

        snap.particles.position = pos   # atooms.system and gsd both save positions from -L/2 to L/2.
        snap.particles.typeid = typeid
        if 'particle.velocity' in variables:
            snap.particles.velocity = vel
        if 'particle.mass' in variables:
            snap.particles.mass = mass
        if 'particle.diameter' in variables:
            snap.particles.diameter = 2 * radius

        self._file.append(snap)
