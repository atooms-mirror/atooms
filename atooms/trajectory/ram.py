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
        self._system.append(copy.copy(system))
        self.steps.append(step)

    def read_sample(self, frame):
        return self._system[frame]


class TrajectoryRam(TrajectoryBase):

    """Lighter in terms of memory. We only store positions."""

    def __init__(self, fname=None, mode='w'):
        TrajectoryBase.__init__(self, fname, mode)
        self._pos = []
        self._species = []
        self._cell = []
        self.mode = mode

    def write_sample(self, system, step):
        try:
            # Ovewrite
            ind = self.steps.index(step)
        except:
            ind = None

        if ind is not None:
            # Overwrite
            self._species[ind] = [p.species for p in copy.copy(system.particle)]
            self._pos[ind] = [p.position for p in copy.copy(system.particle)]
            self._cell[ind] = copy.copy(system.cell)
        else:
            # Append a new frame
            self._species.append([p.species for p in copy.copy(system.particle)])
            self._pos.append([p.position for p in copy.copy(system.particle)])
            self._cell.append(copy.copy(system.cell))

    def read_sample(self, frame):
        particles = []
        for i in range(len(self._pos[frame])):
            # TODO: copy mass and radius as well!
            particles.append(Particle(position=self._pos[frame][i],
                                      species=self._species[frame][i]))
        cell = self._cell[frame]
        return System(particles, cell)

    def __setitem__(self, i, y):
        try:
            step = self.steps[i]
        except IndexError:
            if len(self.steps) > 0:
                step = self.steps[-1]+1
            else:
                step = 0
        self.write(y, step)
