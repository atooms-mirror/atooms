# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import numpy
import logging

from .base import TrajectoryBase
from .xyz import TrajectorySimpleXYZ

log = logging.getLogger(__name__)

class TrajectoryField(TrajectorySimpleXYZ):

    """
    Trajectory of scalar or vector fields based on simple xyz layout.
    """

    suffix = 'xyz'

    def __init__(self, filename, mode='r'):
        TrajectoryBase.__init__(self, filename, mode)
        self.trajectory = open(self.filename, self.mode)
        self.field_name = ''
        self.fmt = ['field']
        if self.mode == 'r':
            # Internal index of lines to seek and tell.
            # We may delay setup, moving to read_init() assuming
            # self.steps becomes a property
            self._setup_index()
            self._setup_steps()

    def read_init(self):
        pass

    def read_sample(self, sample):
        meta = self._read_metadata(sample)
        npart = meta['npart']
        self.trajectory.seek(self._index_sample[sample])
        field = []
        for _ in range(npart):
            data = self.trajectory.readline().strip()
            if len(data.split()) > 1:
                data = data.split()
            field.append(data)
        try:
            return numpy.array(field, dtype=float)
        except ValueError:
            return numpy.array(field, dtype=str)

    def _comment_header(self, step, field):
        return "step:%d columns:%s" % (step, ','.join(self.fmt))

    def write_sample(self, field, step):
        self.trajectory.write("%s\n" % len(field))
        self.trajectory.write(self._comment_header(step, field) + '\n')
        try:
            ndim = len(field[0])
            if type(field[0]) is str:
                ndim = 1
        except TypeError:
            ndim = 1

        if ndim == 1:
            for f in field:            
                self.trajectory.write('%s\n' % f)
        else:
            fmt = ndim*" %g" + "\n"
            for f in field:            
                self.trajectory.write(fmt % tuple(f))
            
    def close(self):
        self.trajectory.close()
