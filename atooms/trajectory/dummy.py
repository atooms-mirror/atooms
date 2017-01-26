# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

from atooms.trajectory.base import TrajectoryBase

class TrajectoryDummy(TrajectoryBase):

    """Dummy Trajectory to illustrate how plugins work."""

    suffix = 'dummy'

    def __init__(self, filename, mode='r'):
        TrajectoryBase.__init__(self, filename, mode)

    def read_sample(self, sample):
        return None

    def write_sample(self, system, step):
        pass

    def close(self):
        pass
