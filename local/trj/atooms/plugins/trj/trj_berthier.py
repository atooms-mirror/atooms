#!/usr/bin/env python
# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Berthier trajectory format. Adapted from Ninarello (->Kob) format."""

import sys
import os
import numpy
import glob
import re

from atooms.trajectory import TrajectoryBase
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System

class TrajectoryBerthier(TrajectoryBase):

    # Warning: pass right basename to get steps
    def __init__(self, dirname, box, basename='co*', mode='r', x=1.0):
        dirname = dirname.strip('/')
        TrajectoryBase.__init__(self, dirname, mode)
        if not os.path.isdir(dirname):
            raise IOError("We expected this to be a dir (%s)" % dirname)
        self.__f_all = glob.glob(dirname + '/%s*' % basename)
        self.__box = box
        self.__x = x
        # TODO: sorting is copied over from lammps, refactor
        file_steps = []
        for f in self.__f_all:
            s = re.search(r'%s_(\d*)' % basename, f)
            if s:
                try:
                    step = int(s.group(1))
                except:
                    step = len(file_steps)+1
                file_steps.append((f, step))
    
        file_steps.sort(key = lambda a : a[1])
        self.__f_all = [a[0] for a in file_steps]
        self.steps = [a[1] for a in file_steps]

    def read_sample(self, sample):
        # Read particles
        particle = []
        with open(self.__f_all[sample]) as fh:            
            for l in fh:
                data = l.split()
                pos = numpy.array(data[:3], dtype=float)
                particle.append(Particle(name='A', id=1, position=pos))
                particle.diameter = float(data[3])

        # Create cell
        cell = Cell(numpy.array(3*[self.__box]))

        # Center and fold positions into central cell        
        for i, p in enumerate(particle):
            p.position -= cell.side / 2
            try:
                p.fold(cell)
            except:
                print i, sample, self.__f_all[sample]
                raise

        return System(particle, cell)

if __name__ == '__main__':
    from atooms.trajectory import  TrajectoryXYZ, TrajectoryHDF5, convert
    for f in sys.argv[2:]:
        with TrajectoryBerthier(f, box=float(sys.argv[1])) as t:
            # Box side is in st_L* file sometimes
            fout = convert(t, TrajectoryHDF5, prefix='conv-')
            #fout = convert(t, TrajectoryXYZ)
