# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

import os
import tarfile
import numpy

from .base import TrajectoryBase
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System

from xml.etree import ElementTree


def map_label_id(names):
    lab = {}
    for i, key in enumerate(sorted(set(names))):
        lab[key] = i
    return lab


class TrajectoryHOOMD(TrajectoryBase):

    """HOOMD format"""

    suffix = 'tgz'

    def __init__(self, fname, mode='r'):

        super(TrajectoryHOOMD, self).__init__(fname, mode)
        base = os.path.basename(fname.split('.tgz')[0])

        if mode == 'r':
            self.__tmp_path = '/tmp/' + base
            try:
                os.makedirs(self.__tmp_path)
            except IOError:
                pass

            tar = tarfile.open(fname)
            file_list = sorted([f.name for f in tar.getmembers()])
            tar.extractall(path=self.__tmp_path)
            tar.close()
            self.__f_frames = [os.path.join(self.__tmp_path, f) for f in file_list]
            # cfg, box, pos, typ, vel = self.__read_one(self.__f_frames[0])

        elif mode == 'w':
            pass

        elif mode == 'w:gz':
            self._file = tarfile.open(fname, "w:gz")

    def read_steps(self):
        steps = []
        for f in sorted(self.__f_frames):
            tree = ElementTree.parse(f)
            root = tree.getroot()
            cfg = root.find('configuration')
            steps.append(int(cfg.attrib['time_step']))

        # First sort them, then subtract out the first step
        steps = sorted(steps)
        return [s - steps[0] for s in steps]

    def __read_one(self, fname):
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        cfg = root.find('configuration')
        box = cfg.find('box')
        pos = cfg.find('position')
        try:
            vel = cfg.find('velocity')
        except:
            vel = None
        typ = cfg.find('type')
        # TODO: use tipify here
        pos_list = [list([float(x) for x in r.split()]) for r in pos.text.strip().split('\n')]
        box_list = [float(box.attrib[r]) for r in ['lx', 'ly', 'lz']]
        typ_list = typ.text.strip().split('\n')
        if vel is not None:
            vel_list = [list([float(x) for x in v.split()]) for v in vel.text.strip().split('\n')]
        else:
            vel_list = None
        return cfg, box_list, pos_list, typ_list, vel_list

    def read_system(self, frame):
        cfg, box, pos, typ, vel = self.__read_one(self.__f_frames[frame])
        if vel is None:
            particle = [Particle(species=t, position=numpy.array(p)) for p, t in zip(pos, typ)]
        else:
            particle = [Particle(species=t, position=numpy.array(p), velocity=numpy.array(v))
                        for p, t, v in zip(pos, typ, vel)]
        cell = Cell(numpy.array(box))
        return System(particle, cell)

    def write_system(self, system, step):
        ndim = len(system.particle[0].position)
        n = len(system.particle)

        if self.mode == 'w':
            fname = self.filename.split('.tgz')[0] + '_step%010d.xml' % step
        else:
            base = os.path.basename(self.filename.split('.tgz')[0])
            fname = base + '_step%010d.xml' % step

        fh = open(fname, 'w')
        fh.write("""\
<?xml version="1.0" encoding="UTF-8"?>
<hoomd_xml version="1.4">
<configuration time_step="%d" dimensions="%d" natoms="%d" >
<box lx="%g" ly="%g" lz="%g"/>
""" % ((step, ndim, n) + tuple(system.cell.side)))

        fh.write("<position num=\"%d\">\n" % n)
        for p in system.particle:
            fh.write((ndim*" %g" + "\n") % tuple(p.position))
        fh.write("</position>\n")

        fh.write("<velocity num=\"%d\">\n" % n)
        for p in system.particle:
            fh.write((ndim*" %g" + "\n") % tuple(p.velocity))
        fh.write("</velocity>\n")

        fh.write("<type num=\"%d\">\n" % n)
        for p in system.particle:
            fh.write(("%s\n") % p.species)
        fh.write("</type>\n")

        fh.write("""\
</configuration>
</hoomd_xml>
""")
        fh.close()

        if self.mode == 'w:gz':
            self._file.add(fname)
            os.remove(fname)

    def close(self):
        # clean up
        if self.mode == 'r':
            import shutil
            shutil.rmtree(self.__tmp_path)
        elif self.mode == 'w:gz':
            self._file.close()
