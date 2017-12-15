# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""PDB format, write-only."""

from .base import TrajectoryBase


class TrajectoryPDB(TrajectoryBase):

    """Trajectory file with PDB layout"""

    suffix = 'pdb'

    def __init__(self, filename, mode='w'):
        super(TrajectoryPDB, self).__init__(filename, mode)
        self.trajectory = open(self.filename, self.mode)

    def write_sample(self, system, step):
        cfg = ''
        cfg += 'MODEL%9i\n' % step
        cfg += ('CRYST1' + (3*'%9.3f') + '     90     90     90 P 1           1\n') % tuple(system.cell.side)
        for i, p in enumerate(system.particle):
            # If particle has a field property we dump it in the pdb file
            if hasattr(p, 'field'):
                x = float(p.field)
            else:
                x = 1.0
            fmt = 'HETATM%5d          %4s    ' + 3*'%8.3f' + 2*'%6.2f' + '          %4s\n'
            cfg += fmt % ((i, p.species) + tuple(p.position) + (1.0, x, p.species))
        self.trajectory.write(cfg)

    def close(self):
        self.trajectory.close()
