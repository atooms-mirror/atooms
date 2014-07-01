# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import gzip
import numpy

from atooms.trajectory  import TrajectoryBase
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System

class TrajectoryXYZ(TrajectoryBase):

    """Trajectory with XYZ layout using memory leightweight indexed access."""

    suffix = 'xyz'

    def __init__(self, filename, mode='r'):
        TrajectoryBase.__init__(self, filename, mode)
        self._timestep = 1.0
        self._cell = None
        self._map_id = [] # list to map numerical ids (indexes) to chemical species (entries)

        if self.mode == 'w':
            self.trajectory = open(self.filename, 'w')

        elif self.mode == 'a':
            self.trajectory = open(self.filename, 'a')

        elif self.mode == 'r':
            ext = os.path.splitext(self.filename)[1]
            if ext == '.gz':
                self.trajectory = gzip.open(self.filename)
            else:
                self.trajectory = open(self.filename, 'r')
            self._setup_index(self.trajectory)
            # Define sample list and fix case when no step is available
            # Should it be in _setup_index()
            if not self.steps:
                self.steps = range(1,len(self._index)+1)
        else:
            raise ValueError('Specify mode (r/w) for file %s (invalid: %s)' % (self.filename, self.mode))

    def _setup_index(self, fh):
        """Sample indexing via tell / seek"""
        self._npart = []
        self._index = []
        self._index_cell = None

        # TODO: index lines at which first particle is found, not comment
        fh.seek(0)
        l = 0
        while True:
            data = fh.readline()

            # There nothing here, file is over
            if not data:
                break

            # Guess this is the number of particles
            try:
                npart = int(data)
            except ValueError:
                # This must be the side of the cell or some empty line
                # TODO: move search for cell in a separate function, more general solution 
                if data:
                    self._index_cell = l
                break

            # Store number of particles and skip npart+1 lines
            self._index.append(l)
            self._npart.append(npart)
            for i in range(npart+1):
                d = fh.readline()
            l = fh.tell()

        # Get number of columns
        self.trajectory.seek(0)
        self.trajectory.readline()
        self.trajectory.readline()
        self._ncols = len(self.trajectory.readline().split())
        self.trajectory.seek(0)

        # Get cell
        self._cell = self._parse_cell()

        # Get list of steps and samples
        self.steps = []
        for i in self._index:
            fh.seek(i)
            fh.readline() # skip Npart
            step = self._parse_step(fh.readline())
            if step is None:
                break
            else:
                self.steps.append(step)

    def rewind(self):
        self.trajectory.seek(0)

    def _parse_step(self, data):
        """Internal xyz method to grab the step index. Can be overwritten in subclasses."""
        try:
            return int(data.split()[-1])
        except:
            return None

    def _parse_cell(self):
        """Internal xyz method to grab the cell. Can be overwritten in subclasses."""
        # This is schizophrenic, move everything related to cell search here
        cell = None
        if self._index_cell:
            self.trajectory.seek(self._index_cell)
            side = numpy.fromstring(self.trajectory.readline(), sep=' ')
            cell = Cell(side)
        return cell

    @property
    def grandcanonical(self):
        if self._grandcanonical is None:
            self._grandcanonical = False
            for i in range(len(self._npart)-1):
                if self._npart[i] != self._npart[i+1]:
                    self._grandcanonical = True
                    break

        # TODO: check what the problem is with GC and species sorting
        # if self._grandcanonical:
        #     raise NotImplementedError('cannot do spcies sorting with GC')
        return self._grandcanonical

    def read_sample(self, sample):
        self.trajectory.seek(self._index[sample])
        self.trajectory.readline() # skip npart
        self.trajectory.readline() # skip comment header

        # TODO: can we fusion these two loops ?
        data = []
        for j in range(self._npart[sample]):
            data.append(self.trajectory.readline().split())

        p = []
        for d in data:
            # Get particle name and id and update local database if needed
            name = d[0]
            if not name in self._map_id:
                self._map_id.append(name)

            # Get positions, velocities, tag and append to particle list
            r = numpy.array(d[1:4], dtype=float)
            v = numpy.zeros(3)
            tag = None
            if len(d)==7:
                v = numpy.array(d[4:7], dtype=float)
            elif len(d)==5:
                tag = d[4]
            elif len(d)>=12: # hackish for RUMD
                v = numpy.array(d[7:10], dtype=float)
            p.append(Particle(name=name, id=None, position=r, velocity=v, tag=tag))

        # Assign ids to particles according to the updated database
        for pi in p:
            pi.id = self._map_id.index(pi.name)

        # TODO: I think the name of the system is irrelevant, move it to last parameter
        return System('unknown', p, self._cell)

    def _comment_header(self, step, system):
        return "Step: %d\n" % step

    def write_sample(self, system, step, ignore=[]):
        self._cell = system.cell
        self.trajectory.write("%8i\n" % len(system.particle))
        self.trajectory.write(self._comment_header(step, system))
        ndim = len(system.particle[0].position)
        if abs(system.particle[0].velocity[0]) < 1e-15 and \
           abs(system.particle[-1].velocity[-1]) < 1e-15:
            for p in system.particle:
                # Check for tag, somewhat hard hack to write voronoi polyehdron
                # TODO: could be improved 
                tag = ''
                if not p.tag is None:
                    tag = 4*" " + p.tag
                self.trajectory.write(("%s"+ndim*" %14.6f"+tag+"\n") % ((p.name,) + tuple(p.position)))
        else:
            for p in system.particle:
                self.trajectory.write(("%s"+2*ndim*" %14.6f"+"\n") % ((p.name,) + tuple(p.position) + tuple(p.velocity)))

    def close(self):
        # As a last thing, write the side of the cell
        if self.mode == 'w':
            if self._cell:
                ndim = len(self._cell.side)
                self.trajectory.write((ndim*" %14.6f"+"\n") % tuple(self._cell.side))
        self.trajectory.close()

# TODO: refactor xyz class to parse a predefined sequence of float/int/str with appropriate tags / metadata as RUMD
class TrajectoryNeighbors(TrajectoryXYZ):

    """Neighbors info"""

    def __init__(self, filename, offset=-1):
        super(TrajectoryNeighbors, self).__init__(filename)
        # TODO: determine minimum value of index automatically
        self.offset = offset # neighbors produced by voronoi are indexed from 1

    def read_sample(self, sample):
        self.trajectory.seek(self._index[sample])
        self.trajectory.readline() # skip npart
        self.trajectory.readline() # skip comment header
        p = []
        for j in range(self._npart[sample]):
            data = self.trajectory.readline().split()
            neigh = numpy.array(data, dtype=int)
            p.append(neigh+self.offset)
        return p

class TrajectoryPDB(TrajectoryBase):

    """Trajectory file with PDB layout"""

    suffix = 'pdb'

    def __init__(self, filename, mode='r'):
        super(TrajectoryPDB, self).__init__(filename)
        if mode == 'w':
            self.trajectory = open(self.filename, 'w')

    def write_sample(self, system, step, ignore=[]):
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


class TrajectoryXYZIkeda2(TrajectoryXYZ):

    """Trajectory with indexed XYZ layout from Atsushi Ikeda. Assume one component"""

    def _setup_index(self, fh):

        self._index = []
        self._index_cell = None

        fh.seek(0)
        data = self.trajectory.readline().split()
        npart = int(data[0])
        ncfg = int(data[1])
        L = numpy.array([data[2]] * 3, dtype=float)
        self._npart = [npart] * ncfg
        self._cell = Cell(L)

        for i in range(ncfg):
            # This is a float in Atsushi's format
            step = float(self.trajectory.readline())
            self.steps.append(step)
            self._index.append(fh.tell())
            for i in range(npart):
                fh.readline()

        # Atsushi stores the actual time (istep*dt).
        # To get it right we define an artificial timestep as the time difference
        # between the two first samples and then redefine steps = samples
        self._timestep = self.steps[1] - self.steps[0] # assume linear grid
        self.steps = range(len(self.steps))

    def read_sample(self, sample):

        self.trajectory.seek(self._index[sample])
        p = []
        for j in range(self._npart[sample]):
            d = self.trajectory.readline().split()
            r = numpy.array(d[1:4], dtype=float)
            p.append(Particle(name='A', id=1, position=r)) #- self._cell.side/2.0)))
        return System('unknown', p, self._cell)





