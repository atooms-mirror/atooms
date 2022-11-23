# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

import numpy
import re
from .base import TrajectoryBase
from .utils import gopen
from atooms.core.utils import tipify
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System


class TrajectorySimpleXYZ(TrajectoryBase):

    """
    Simple implementation of the xyz layout (https://en.wikipedia.org/wiki/XYZ_file_format)

    It uses a memory light-weight indexed access.
    """

    suffix = 'xyz'

    def __init__(self, filename, mode='r'):
        super(TrajectorySimpleXYZ, self).__init__(filename, mode)
        self._cell = None
        # Trajectory file handle
        self._file = gopen(self.filename, self.mode)
        if self.mode == 'r':
            # Internal index of lines to seek and tell.
            # We may delay setup, moving to read_init() assuming
            # self.steps becomes a property
            self._setup_index()

    def read_steps(self):
        steps = []
        for frame in range(len(self._index_frame)):
            meta = self._read_comment(frame)
            try:
                steps.append(meta['step'])
            except KeyError:
                # If no step info is found, we add steps sequentially
                steps.append(frame+1)
        return steps

    def _setup_index(self):
        """Sample indexing via tell / seek"""
        self._index_frame = []
        self._index_header = []
        self._index_cell = None
        self._file.seek(0)
        while True:
            line = self._file.tell()
            data = self._file.readline().strip()

            # We break if file is over or we found an empty line
            if not data:
                break

            # The first line contains the number of particles
            # If something goes wrong, this could be the last line
            # with the cell side (Lx,Ly,Lz) and we parse it some other way
            try:
                npart = int(data)
                self._index_header.append(line)
            except ValueError:
                self._index_cell = line
                break

            # Skip npart+1 lines
            _ = self._file.readline()
            self._index_frame.append(self._file.tell())
            for _ in range(npart):
                self._file.readline()

    def _read_comment(self, frame):
        """Internal xyz method to get header metadata from comment line of given `frame`.

        We assume metadata format is a space-separated sequence of
        comma separated entries such as:

        columns:id,x,y,z step:10
        columns=id,x,y,z step=10
        """
        # Go to line and skip Npart info
        self._file.seek(self._index_header[frame])
        npart = int(self._file.readline())
        data = self._file.readline()

        # Fill metadata dictionary
        meta = {}
        meta['npart'] = npart
        for e in data.split():
            s = re.search(r'(\S+)\W*[=:]\W*(\S+)', e)
            if s is not None:
                tag, data = s.group(1), s.group(2)
                # Remove dangling commas
                data = data.strip(',')
                # If there are commas, this is a list, else a scalar.
                # We convert the string to appropriate types
                if ',' in data:
                    meta[tag] = [tipify(_) for _ in data.split(',')]
                else:
                    meta[tag] = tipify(data)
        return meta

    def read_init(self):
        # Grab cell from the end of file if it is there
        try:
            side = self._read_comment(0)['cell']
            self._cell = Cell(side)
        except KeyError:
            self._cell = self._parse_cell()

    def _parse_cell(self):
        """Internal emergency method to grab the cell."""
        cell = None
        if self._index_cell:
            self._file.seek(self._index_cell)
            side = numpy.fromstring(self._file.readline(), sep=' ')
            cell = Cell(side)
        return cell

    def read_system(self, frame):
        meta = self._read_comment(frame)
        self._file.seek(self._index_frame[frame])

        # Read particles
        particle = []
        for _ in range(meta['npart']):
            data = self._file.readline().strip().split()
            species = data[0]
            r = numpy.array(data[1:4], dtype=float)
            particle.append(Particle(species=species, position=r))

        # Read cell
        try:
            side = meta['cell']
            self._cell = Cell(side)
        except KeyError:
            pass

        return System(particle, self._cell)

    def _comment_header(self, step, system):
        fmt = "step:%d columns:id,x,y,z" % step
        if system.cell is not None:
            fmt += " cell:" + ','.join(['%s' % x for x in system.cell.side])
        return fmt

    def write_system(self, system, step):
        self._file.write("%s\n" % len(system.particle))
        self._file.write(self._comment_header(step, system) + '\n')
        ndim = len(system.particle[0].position)
        fmt = "%s" + ndim*" %14.6f" + "\n"
        for p in system.particle:
            self._file.write(fmt % ((p.species,) + tuple(p.position)))

    def close(self):
        self._file.close()
