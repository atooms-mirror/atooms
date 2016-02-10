# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import gzip
import numpy
import re

from atooms.trajectory  import TrajectoryBase
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System

class TrajectoryXYZ(TrajectoryBase):

    """Trajectory with XYZ layout using memory leightweight indexed access."""

    suffix = 'xyz'

    # TODO: move all file init to write_init. Rename trajectory fh or something like this.

    def __init__(self, filename, mode='r'):
        TrajectoryBase.__init__(self, filename, mode)
        # This is the default column format.
        # TODO: Using vx, vy, vz in the header will grab the velocities
        self.fmt = ['id', 'x', 'y', 'z']
        self._timestep = 1.0
        self._cell = None
        self._map_id = [] # list to map numerical ids (indexes) to chemical species (entries)
        self._min_id = 1 # minimum integer for ids, can be modified by subclasses
        self._step = 0

        if self.mode == 'w':
            self.trajectory = open(self.filename, 'w')

        elif self.mode == 'a':
            self.trajectory = open(self.filename, 'a')

        elif self.mode == 'r':
            # TODO: allow reading xyz file from stdin. Seek wont work though, need to pass through tmp file
            ext = os.path.splitext(self.filename)[1]
            if ext == '.gz':
                self.trajectory = gzip.open(self.filename)
            else:
                self.trajectory = open(self.filename, 'r')

            self._setup_index(self.trajectory)
            # Define sample list and fix case when no step is available
            # Should it be in _setup_index()
            if len(self.steps) == 0:
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

        # Get list of steps, samples and cell
        self.steps = []
        for i in self._index:
            fh.seek(i)
            fh.readline() # skip Npart
            meta = self._parse_header(fh.readline())
            self.steps.append(meta['step'])
        
        # Grab it from the end of file in case it is there
        if meta['cell'] is not None:
            self._cell = Cell(meta['cell'])
        else:
            self._cell = self._parse_cell()

    def rewind(self):
        self.trajectory.seek(0)
        self._step = 0

    def _parse_header(self, data):
        """Internal xyz method to get header metadata."""
        meta = {'step':None, 'cell':None}
        grabber = {'step':lambda x: int(x), 
                   'cell':lambda x: map(float, x.split(','))}
        for m in meta:
            p = re.compile(r'%s\s*[=:]\s*(\S*)\s' % m, re.IGNORECASE)
            s = p.search(data)
            if s is not None:     
                meta[m] = grabber[m](s.group(1))

        # Fix step
        if meta['step'] is None:
            try:
                n = int(data.split()[-1])
            except:
                self._step += 1
                n = self._step
            meta['step'] = n
        return meta

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

        self._sampledata = []
        for j in range(self._npart[sample]):
            line = self.trajectory.readline().split()
            # Unpack into a dictionary according to the specific
            # format (self.fmt). This dict can then be used by a
            # subclass to modify or override system properties.  
            self._sampledata.append({})
            for l, f in zip(line, self.fmt):
                self._sampledata[-1][f] = l
            # Put what remains in line in the tag
            if len(line) > len(self.fmt):
                self._sampledata[-1]['tag'] = line[len(self.fmt):]
            # The special greedy '*' in the last field matches all
            # that remains so we redefine this here
            lastfmt = self.fmt[-1]
            if len(line) > len(self.fmt) and '*' in lastfmt:
                self._sampledata[-1][lastfmt] = ' '.join(line[len(self.fmt)-1:])

        p = []
        for data in self._sampledata:
            # Get particle name and id and update local database if needed
            name = data['id']
            if not name in self._map_id:
                self._map_id.append(name)

            # Get positions and velocities
            r = numpy.array([data['x'], data['y'], data['z']], dtype=float)
            try:
                v = numpy.array([data['vx'], data['vy'], data['vz']], dtype=float)
            except KeyError:
                v = numpy.zeros(3)

            # Try to grab left over as a tag.
            try:
                tag = data['tag']
            except KeyError:
                tag = None

            p.append(Particle(name=name, id=None, position=r, velocity=v, tag=tag))

        # Assign ids to particles according to the updated database
        # The minimum id is set by self._min_id
        for pi in p:
            pi.id = self._map_id.index(pi.name) + self._min_id

        return System(p, self._cell)

    def _comment_header(self, step, system):
        if system.cell is not None:
            fmt = "Step:%d Cell:%s Columns:%s\n"
            return fmt % (step, ','.join(map(lambda x: '%s' % x, system.cell.side)),
                          ','.join(self.fmt))
        else:
            fmt = "Step:%d Columns:%s\n"
            return fmt % (step, ','.join(self.fmt))

    def write_sample(self, system, step):
        self._cell = system.cell
        self.trajectory.write("%8i\n" % len(system.particle))
        self.trajectory.write(self._comment_header(step, system))
        ndim = len(system.particle[0].position)
        if (abs(system.particle[0].velocity[0]) < 1e-15 and \
           abs(system.particle[-1].velocity[-1]) < 1e-15) or \
            (not 'vx' in self.fmt):
            for p in system.particle:
                # Check for tag, somewhat hard hack to write voronoi polyehdron
                # TODO: could be improved 
                tag = ''
                if not p.tag is None and not isinstance(p.tag, list):
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

# TODO: refactor xyz class to parse a predefined sequence of float/int/str with appropriate tags / metadata as RUMD -> drop cell
class TrajectoryNeighbors(TrajectoryXYZ):

    """Neighbors info"""

    def __init__(self, filename, offset=1):
        super(TrajectoryNeighbors, self).__init__(filename)
        # TODO: determine minimum value of index automatically
        self._offset = offset # neighbors produced by voronoi are indexed from 1
        self._netwon3 = False
        self._netwon3_message = False

    def _parse_header(self, data):
        """Internal xyz method to get header metadata."""        
        meta = {'time':None, 'cell':None}
        grabber = {'time':lambda x: int(x)}
        for m in meta:
            p = re.compile(r'%s\s*[=:]\s*(\S*)\s' % m, re.IGNORECASE)
            s = p.search(data)
            if s is not None:     
                meta[m] = grabber[m](s.group(1))

        # Fix step
        meta['step'] = meta['time']
        if meta['step'] is None:
            try:
                n = int(data.split()[-1])
            except:
                self._step += 1
                n = self._step
            meta['step'] = n
        return meta

    def read_sample(self, sample):
        self.trajectory.seek(self._index[sample])
        self.trajectory.readline() # skip npart
        self.trajectory.readline() # skip comment header
        s = System()
        s.neighbors = []
        for j in range(self._npart[sample]):
            data = self.trajectory.readline().split()
            neigh = numpy.array(data, dtype=int)
            s.neighbors.append(neigh-self._offset)

        # Ensure III law Newton.
        # If this is ok on first sample we skip it for the next ones
        # if not self._netwon3:
        #     self._netwon3 = True
        #     for i, ilist in enumerate(p):
        #         for j in ilist:
        #             if not i in p[j]:
        #                 p[j].append(i)
        #                 self._netwon3 = False
        #     if not self._netwon3 and not self._netwon3_message:
        #         print 'Warning: enforcing 3rd law of Newton...'
        return s

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
        return System(p, self._cell)
