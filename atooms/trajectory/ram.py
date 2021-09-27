"""Store trajectory in memory (it can be huge)."""

import copy
from .base import TrajectoryBase


class TrajectoryRam(TrajectoryBase):

    """
    Store trajectory in RAM

    The read_system method of this class conforms with the normal
    Trajectory behavior, i.e. a copy of the system is returned when
    requesting the same frame multiple times.
    """

    def __init__(self, filename=None, mode='w'):
        super(TrajectoryRam, self).__init__(filename, mode)
        self._system = []
        self._overwrite = True

    def write_system(self, system, step):
        if step in self.steps:
            ind = self.steps.index(step)
            self._system[ind].update(system)
        else:
            self._system.append(copy.deepcopy(system))

    def read_system(self, frame):
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
    returns views on the System when calling read_system(), and not
    copies. Thus modifications to the read system object will be
    propagated to the trajectory.
    """

    def read_system(self, frame):
        return self._system[frame]


# This is maintanined for backward compatibility
TrajectoryRamFull = TrajectoryRamView
