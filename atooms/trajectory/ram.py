"""Store trajectory in memory (it can be huge)."""

import copy
from atooms.system.system import System
from atooms.system.particle import Particle
from .base import TrajectoryBase


class TrajectoryRam(TrajectoryBase):

    """
    Store trajectory in RAM

    The read_sample method of this class conforms with the normal
    Trajectory behavior, i.e. a copy of the system is returned when
    requesting the same frame multiple times.
    """

    def __init__(self, filename=None, mode='w'):
        TrajectoryBase.__init__(self, filename, mode)
        self._system = []
        self._overwrite = True

    def write_sample(self, system, step):
        if step in self.steps:
            ind = self.steps.index(step)
            self._system[ind].update(system)
        else:
            self._system.append(copy.deepcopy(system))
            self.steps.append(step)

    def read_sample(self, frame):
        return copy.deepcopy(self._system[frame])

    def __setitem__(self, i, value):
        try:
            step = self.steps[i]
        except IndexError:
            if len(self.steps) > 0:
                step = self.steps[-1]+1
            else:
                step = 0
        self.write(value, step)


class TrajectoryRamView(TrajectoryRam):

    """
    This class deviates from the normal Trajectory behavior in that it
    returns views on the System when calling read_sample(), and not
    copies. Thus modifications to the read system object will be
    propagated to the trajectory.
    """

    def read_sample(self, frame):
        return self._system[frame]

# This is maintanined for backward compatibility
TrajectoryRamFull = TrajectoryRamView
