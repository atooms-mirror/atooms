#!/usr/bin/env python
# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Ninarello trajectory format. Adapted from Kob format."""

import sys
import os
import numpy
import glob
import re

from base import TrajectoryBase
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System

class TrajectoryNinarello(TrajectoryBase):

    def __init__(self, dirname, basename='Cnf', mode='r', x=1.0, box=None):

        TrajectoryBase.__init__(self, dirname, mode)
        if not os.path.isdir(dirname):
            raise IOError("We expected this to be a dir (%s)" % dirname)
        self.__f_all = glob.glob(dirname + '/%s-*' % basename)
        self.__box = box
        self.__x = x
        # TODO: sorting is copied over from lammps, refactor
        file_steps = []
        for f in self.__f_all:
            s = re.search(r'%s-(\d*)' % basename, f)
            if s:
                step = int(s.group(1))
                file_steps.append((f, step))
        file_steps.sort(key = lambda a : a[1])
        self.__f_all = [a[0] for a in file_steps]
        self.steps = [a[1] for a in file_steps]

    def read_sample(self, sample):
        particle = []
        with open(self.__f_all[sample]) as fh:            
            first = True
            for l in fh:
                if first:
                    self.__box = float(l)
                    first = False
                else:
                    pos = numpy.array(l.split()[:3], dtype=float)
                    particle.append(Particle(name='A', id=1, position=pos))
        if self.__x < 1.0:
            for i, p in enumerate(particle):
                if i >= int(len(particle)*self.__x):
                    p.name = 'B'
                    p.id = 2
        if self.__box:
            cell = Cell(numpy.array(3*[self.__box]))
            for i, p in enumerate(particle):
                try:
                    p.fold(cell)
                except:
                    print i, sample, self.__f_all[sample]
                    raise
        else:
            cell = None
        return System(particle, cell)

if __name__ == '__main__':
    #'../data/home/medusa/kob/lj_pin/N300/T0.40/p000/run01_all/configs'
    from atooms.trajectory import TrajectoryBase, TrajectoryHDF5, convert
    #side = float(sys.argv[1])
    for f in sys.argv[1:]:
        #with TrajectoryNinarello(f, x=0.8) as t:
        with TrajectoryNinarello(f) as t:
            fout = convert(t, TrajectoryHDF5)            

        # import re
        # p=re.match('.*[tT](\d\.\d*)', f)
        # if p:
        #     T=p.group(1)
        # else:
        #     T=0.0
        # print 'converted %s (T=%s)' % (f, T)
        # with open(fout + '.pp.thermo.ave', 'w') as fh:
        #     fh.write('Temp=%s' % T)
