# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""PDB format, write-only."""

import os
from base import TrajectoryBase

class TrajectoryPDB(TrajectoryBase):

    """Trajectory file with PDB layout"""

    suffix = 'pdb'

    def __init__(self, filename, mode='r'):
        super(TrajectoryPDB, self).__init__(filename)
        if mode == 'w':
            self.trajectory = open(self.filename, 'w')

    def write_sample(self, system, step):
        self._system = system
        cfg = ''
        cfg += 'MODEL%9i\n' % step
        cfg += ('CRYST1' + (3*'%9.3f') + '     90     90     90 P 1           1\n') % tuple(system.cell.side)
        for i, p in enumerate(system.particle):
            if p.tag is None:
                p.tag = 1.0
            fmt = 'HETATM%5d          %4s    ' + 3*'%8.3f' + 2*'%6.2f' + '          %4s\n'
            cfg += fmt % ((i, p.name) + tuple(p.position) + (1.0, float(p.tag), p.name))
        self.trajectory.write(cfg)

    def close(self):
        self.trajectory.close()
