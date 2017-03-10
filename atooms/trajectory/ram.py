import copy
from atooms.system.system import System
from atooms.system.particle import Particle
from atooms.system.cell import Cell
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
        self._ids = []
        self._name = []
        self._cell = []
        self.mode = mode

    def write_sample(self, system, step):
        self._name.append([p.name for p in copy.copy(system.particle)])
        self._ids.append([p.id for p in copy.copy(system.particle)])
        self._pos.append([p.position for p in copy.copy(system.particle)])
        self._cell.append(copy.copy(system.cell))

    def read_sample(self, frame):
        particles = []
        for i in range(len(self._pos[frame])):
            particles.append(Particle(position=self._pos[frame][i],
                                      id=self._ids[frame][i],
                                      name=self._name[frame][i]))
        cell = self._cell[frame]
        return System(particles, cell)
