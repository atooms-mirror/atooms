# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import gzip
import numpy
import re

from base import TrajectoryBase
from atooms.utils import tipify 
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System

class TrajectoryXYZBase(TrajectoryBase):

    suffix = 'xyz'

    def __init__(self, filename, mode='r'):
        TrajectoryBase.__init__(self, filename, mode)
        self.fmt = []
        self.cbk_write = []
        self.cbk_read = []
        self.fh = open(self.filename, self.mode)

    def write_sample(self, system, step):
        # Check that all arrays in data have the same length
        # TODO: hey whats the status of this cbk thing?
        data = [[f(p) for p in system.particle] for f in self.cbk_write]
        nlines = len(data[0])
        ncols = len(data)
        lengths_ok = map(lambda x: len(x) == nlines, data)
        if not all(lengths_ok):
            raise ValueError('All arrays must have the same length')

        # Write in xyz format
        self.fh.write('%d\n' % nlines)
        # Comment line: concatenate metadata
        line = 'step: %d; ' % step 
        line += 'columns:' + ','.join(self.fmt)
        self.fh.write(line + '\n')
        # Write data. This is inefficient but
        # we cannot use numpy.savetxt because there is no append mode.
        fmt = '%s ' * len(data)
        fmt = fmt[:-1] + '\n'
        for i in range(nlines):
            self.fh.write(fmt % tuple([data[j][i] for j in range(ncols)]))


class TrajectoryXYZ(TrajectoryBase):

    """Trajectory with XYZ layout using memory leightweight indexed access."""

    suffix = 'xyz'

    # TODO: move all file init to write_init. Rename trajectory fh or something like this.

    def __init__(self, filename, mode='r', tags={}, fmt=['id','x','y','z']):
        TrajectoryBase.__init__(self, filename, mode)
        # This is the default column format.
        # TODO: Using vx, vy, vz in the header will grab the velocities
        self.tags = tags
        self.fmt = fmt
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
            # Internal index of lines to seek and tell.
            # We may delay setup, moving to read_init() assuming
            # self.steps becomes a property
            self._setup_index()
            self._setup_steps()
        else:
            raise ValueError('Specify mode (r/w) for file %s (invalid: %s)' % (self.filename, self.mode))

    def rewind(self):
        self.trajectory.seek(0)
        self._step = 0

    def _setup_index(self):
        """Sample indexing via tell / seek"""
        self._npart = []
        self._index = []
        self._index_cell = None
        fh = self.trajectory

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

    def _setup_steps(self):
        """Define sample list."""
        # Get list of steps, samples and cell
        self.steps = []
        for sample in range(len(self._index)):
            meta = self._read_metadata(sample)
            try:
                self.steps.append(meta['step'])
            except KeyError:
                # If no step info is found, we add steps sequentially
                self.steps.append(sample+1)

    def _read_metadata(self, sample):
        """Internal xyz method to get header metadata from comment line of given sample."""
        # Get comment line
        self.trajectory.seek(self._index[sample])
        self.trajectory.readline() # skip Npart
        line = self.trajectory.readline()
        # We assume metadata fmt is a space-separated sequence of
        # comma separated entries such as
        #   columns:id,x,y,z step:10
        #   columns=id,x,y,z step=10
        meta = {}
        for e in line.split():
            s = re.search(r'(\S*)[=:](\S*)', e)
            if s is not None:
                tag, data = s.group(1), s.group(2)
                # Remove dangling commas
                data = data.strip(',')
                # If there are commas, this is a list, else scalar.
                # We convert the string to appropriate types
                if ',' in data:
                    meta[tag] = map(tipify, data.split(','))
                else:
                    meta[tag] = tipify(data)

        # Apply an alias dict to tags, e.g. to add step if Cfg was found instead
        for alias, tag in self.tags.items():
            try:
                meta[tag] = meta[alias]
            except KeyError:
                pass

        return meta

    def read_init(self):
        # Grab cell from the end of file if it is there
        try:
            side = self._read_metadata(0)['cell']
            self._cell = Cell(side)
        except KeyError:
            self._cell = self._parse_cell()

    # def _parse_header(self, data):
    #     """Internal xyz method to get header metadata."""
    #     meta = {'step':None, 'cell':None, 'columns':None}
    #     grabber = {'step':lambda x: int(x),
    #                'cell':lambda x: map(float, x.split(',')),
    #                'columns':lambda x: x.split(',')}
    #     for m in meta:
    #         p = re.compile(r'%s\s*[=:]\s*(\S*)\s' % m, re.IGNORECASE)
    #         s = p.search(data)
    #         if s is not None:     
    #             meta[m] = grabber[m](s.group(1))

    #     # Fix step
    #     if meta['step'] is None:
    #         try:
    #             n = int(data.split()[-1])
    #         except:
    #             self._step += 1
    #             n = self._step
    #         meta['step'] = n
    #     return meta

    def _parse_cell(self):
        """Internal xyz method to grab the cell. Can be overwritten in subclasses."""
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

        # Fix the mass
        meta = self._read_metadata(sample)
        try:
            mass = map(float, meta['mass'])
            mass_db = {}           
            for k, e in zip(self._map_id, mass):
                mass_db[k] = e
            for pi in p:
                pi.mass = mass_db[pi.name]
        except KeyError:
            pass

        return System(p, self._cell)

    def _comment_header(self, step, system):
        if system.cell is not None:
            fmt = "step:%d cell:%s columns:%s\n"
            return fmt % (step, ','.join(map(lambda x: '%s' % x, system.cell.side)),
                          ','.join(self.fmt))
        else:
            fmt = "step:%d columns:%s\n"
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


class TrajectoryNeighbors(TrajectoryXYZ):

    """Neighbors info"""

    def __init__(self, filename, offset=1):
        super(TrajectoryNeighbors, self).__init__(filename)
        self.tags['step': 'time']
        # TODO: determine minimum value of index automatically
        self._offset = offset # neighbors produced by voronoi are indexed from 1
        self._netwon3 = False
        self._netwon3_message = False

    # TODO: possible regression here if no 'time' tag is found
    # def _parse_header(self, data):
    #     """Internal xyz method to get header metadata."""        
    #     meta = {'time':None, 'cell':None}
    #     grabber = {'time':lambda x: int(x)}
    #     for m in meta:
    #         p = re.compile(r'%s\s*[=:]\s*(\S*)\s' % m, re.IGNORECASE)
    #         s = p.search(data)
    #         if s is not None:     
    #             meta[m] = grabber[m](s.group(1))

    #     # Fix step
    #     meta['step'] = meta['time']
    #     if meta['step'] is None:
    #         try:
    #             n = int(data.split()[-1])
    #         except:
    #             self._step += 1
    #             n = self._step
    #         meta['step'] = n
    #     return meta

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

if __name__ == '__main__':
    import atooms.trajectory as trj
    t=trj.Trajectory('/home/coslo/projects/polydisperse_swap/data/misaki/const_volume/EQ/N8000/phi0.640/conv-configs/config.h5')
    t1=trj.TrajectoryXYZBase('/tmp/1.xyz', 'w')
    t1.fmt=['radius']
    t1.cbk_write=[lambda x: x.radius]
    t1.write(t[0], 0)
