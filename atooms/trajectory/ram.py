"""Store trajectory in memory (can be huge)."""

import copy
from atooms.system.system import System
from atooms.system.particle import Particle
from .base import TrajectoryBase


class TrajectoryRamFull(TrajectoryBase):

    def __init__(self, fname=None, mode='w'):
        TrajectoryBase.__init__(self, fname, mode)
        self._system = []
        self.mode = mode

    def write_sample(self, system, step):
        self._system.append(copy.deepcopy(system))
        self.steps.append(step)

    def read_sample(self, frame):
        return self._system[frame]

    def __setitem__(self, i, value):
        try:
            step = self.steps[i]
        except IndexError:
            if len(self.steps) > 0:
                step = self.steps[-1]+1
            else:
                step = 0
        self.write(value, step)


class TrajectoryRam(TrajectoryBase):

    """
    Store trajectory in RAM.

    Somewhat lighter in terms of memory than TrajectoryRamFull. We do
    not store velocities.
    """

    def __init__(self, fname=None, mode='w'):
        TrajectoryBase.__init__(self, fname, mode)
        self._pos = []
        self._vel = []
        self._species = []
        self._cell = []
        self._radius = []
        self._mass = []
        self._thermostat = []
        self._barostat = []
        self._reservoir = []
        self.mode = mode

    def write_sample(self, system, step):
        try:
            # Overwrite
            ind = self.steps.index(step)
        except:
            ind = None

        particle = copy.copy(system.particle)
        if ind is not None:
            # Overwrite
            self._radius[ind] = [p.radius for p in particle]
            self._mass[ind] = [p.mass for p in particle]
            self._species[ind] = [p.species for p in particle]
            self._pos[ind] = [p.position for p in particle]
            self._vel[ind] = [p.velocity for p in particle]
            self._cell[ind] = copy.copy(system.cell)
            # Store optional system objects
            if system.thermostat is not None:
                self._thermostat[ind] = copy.copy(system.thermostat)
            if system.barostat is not None:
                self._barostat[ind] = copy.copy(system.barostat)
            if system.reservoir is not None:
                self._reservoir[ind] = copy.copy(system.reservoir)
        else:
            # Append a new frame
            self._radius.append([p.radius for p in particle])
            self._mass.append([p.mass for p in particle])
            self._species.append([p.species for p in particle])
            self._pos.append([p.position for p in particle])
            self._vel.append([p.velocity for p in particle])
            self._cell.append(copy.copy(system.cell))
            # Add optional system objects
            if system.thermostat is not None:
                self._thermostat.append(copy.copy(system.thermostat))
            if system.barostat is not None:
                self._barostat.append(copy.copy(system.barostat))
            if system.reservoir is not None:
                self._reservoir.append(copy.copy(system.reservoir))

    def read_sample(self, frame):
        particles = []
        for i in range(len(self._pos[frame])):
            particles.append(Particle(position=self._pos[frame][i],
                                      velocity=self._vel[frame][i],
                                      species=self._species[frame][i],
                                      mass=self._mass[frame][i],
                                      radius=self._radius[frame][i]))
        cell = self._cell[frame]
        s = System(particles, cell)

        # Add optional system objects
        if len(self._thermostat) > 0:
            s.thermostat = self._thermostat[frame]
        if len(self._barostat) > 0:
            s.barostat = self._barostat[frame]
        if len(self._reservoir) > 0:
            s.reservoir = self._reservoir[frame]

        return s

    def __setitem__(self, i, value):
        try:
            step = self.steps[i]
        except IndexError:
            if len(self.steps) > 0:
                step = self.steps[-1]+1
            else:
                step = 0
        self.write(value, step)
