# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

import numpy
import re
import logging

from .base import TrajectoryBase
from .utils import gopen
from atooms.core.utils import tipify
from atooms.system.particle import Particle, distinct_species
from atooms.system.cell import Cell
from atooms.system import System

log = logging.getLogger(__name__)


class TrajectorySimpleXYZ(TrajectoryBase):

    """
    Trajectory with simple xyz layout.

    It uses a memory light-weight indexed access.
    """

    suffix = 'xyz'

    def __init__(self, filename, mode='r'):
        TrajectoryBase.__init__(self, filename, mode)
        self._cell = None
        self.trajectory = open(self.filename, self.mode)
        if self.mode == 'r':
            # Internal index of lines to seek and tell.
            # We may delay setup, moving to read_init() assuming
            # self.steps becomes a property
            self._setup_index()
            self._setup_steps()

    def _setup_index(self):
        """Sample indexing via tell / seek"""
        self._index_frame = []
        self._index_header = []
        self._index_cell = None
        self.trajectory.seek(0)
        while True:
            line = self.trajectory.tell()
            data = self.trajectory.readline().strip()

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
            _ = self.trajectory.readline()
            self._index_frame.append(self.trajectory.tell())
            for i in range(npart):
                _ = self.trajectory.readline()

    def _read_metadata(self, frame):
        """Internal xyz method to get header metadata from comment line of given `frame`.

        We assume metadata format is a space-separated sequence of
        comma separated entries such as:

        columns:id,x,y,z step:10
        columns=id,x,y,z step=10
        """
        # Go to line and skip Npart info
        self.trajectory.seek(self._index_header[frame])
        npart = int(self.trajectory.readline())
        data = self.trajectory.readline()

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

    def _setup_steps(self):
        """Find steps list."""
        self.steps = []
        for frame in range(len(self._index_frame)):
            meta = self._read_metadata(frame)
            try:
                self.steps.append(meta['step'])
            except KeyError:
                # If no step info is found, we add steps sequentially
                self.steps.append(frame+1)

    def read_init(self):
        # Grab cell from the end of file if it is there
        try:
            side = self._read_metadata(0)['cell']
            self._cell = Cell(side)
        except KeyError:
            self._cell = self._parse_cell()

    def _parse_cell(self):
        """Internal emergency method to grab the cell."""
        cell = None
        if self._index_cell:
            self.trajectory.seek(self._index_cell)
            side = numpy.fromstring(self.trajectory.readline(), sep=' ')
            cell = Cell(side)
        return cell

    def read_sample(self, frame):
        meta = self._read_metadata(frame)
        self.trajectory.seek(self._index_frame[frame])

        # Read particles
        particle = []
        for _ in range(meta['npart']):
            data = self.trajectory.readline().strip().split()
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

    def write_sample(self, system, step):
        self.trajectory.write("%s\n" % len(system.particle))
        self.trajectory.write(self._comment_header(step, system) + '\n')
        ndim = len(system.particle[0].position)
        fmt = "%s" + ndim*" %14.6f" + "\n"
        for p in system.particle:
            self.trajectory.write(fmt % ((p.species,) + tuple(p.position)))

    def close(self):
        self.trajectory.close()


# Format callbacks

def update_radius(particle, data, meta):
    particle.radius = float(data[0])
    return data[1:]

def update_tag(particle, data, meta):
    # TODO: what is this???
    particle.tag = data[0:]
    return data[1:]

def update_species(particle, data, meta):
    particle.species = data[0]
    return data[1:]

def update_position(particle, data, meta):
    ndim = meta['ndim']
    # It is crucial to assing position, not to use the slice!
    # Otherwise we get a reference, not a copy.
    particle.position = numpy.array(data[0:ndim], dtype=float)
    return data[ndim:]

def update_velocity(particle, data, meta):
    ndim = meta['ndim']
    particle.velocity = numpy.array(data[0:ndim], dtype=float)
    return data[ndim:]

def _optimize_fields(fields):
    if 'x' in fields:
        fields[fields.index('x')] = 'pos'
    if 'vx' in fields:
        fields[fields.index('vx')] = 'vel'
    for tag in ['y', 'z', 'vy', 'vz']:
        if tag in fields:
            fields.remove(tag)
    return fields


class TrajectoryXYZ(TrajectoryBase):

    """
    Trajectory with XYZ layout using memory leightweight indexed
    access.
    """

    suffix = 'xyz'
    callback_read = {'species': update_species,
                     'type': update_species,  # alias
                     'name': update_species,  # alias
                     'id': update_species,  # alias
                     'tag': update_tag,
                     'radius': update_radius,
                     'pos': update_position,
                     'vel': update_velocity}

    def __init__(self, filename, mode='r', alias=None, fields=None):
        TrajectoryBase.__init__(self, filename, mode)
        if alias is None:
            alias = {}
        # TODO: actualize fields on reading if found and not given on input
        # TODO: clarify fields / _fields handling
        if fields is None:
            fields = ['id', 'pos']
        self.fields = fields
        self._fields = None
        self._fields_float = True
        self._done_format_setup = False
        self.alias = alias
        self.shortcuts = {'pos': 'position',
                          'x': 'position[0]',
                          'y': 'position[1]',
                          'z': 'position[2]',
                          'vel': 'velocity',
                          'vx': 'velocity[0]',
                          'vy': 'velocity[1]',
                          'vz': 'velocity[2]',
                          'id': 'species',
                          'type': 'species'}
        self._cell = None
        self.trajectory = gopen(self.filename, self.mode)
        if self.mode == 'r':
            # Internal index of lines to seek and tell.
            # We may delay setup, moving to read_init() assuming
            # self.steps becomes a property
            self._setup_index()
            # Warning: setting up steps require aliases to be defined in
            # init and not later.
            self._setup_steps()

    def _setup_format(self):
        if not self._done_format_setup:
            self._done_format_setup = True
            # %g allows to format both float and int but it's 2x slower.
            # This switch is for performance
            if self._fields_float:
                _fmt = '%.' + str(self.precision) + 'f'
            else:
                _fmt = '%g'

            def array_fmt(arr):
                """Remove commas and [] from numpy array repr."""
                # Passing a scalar will trigger an error (gotcha: even
                # when casting numpy array to list, the elements remain of
                # numpy type and this function gets called! (4% slowdown)
                try:
                    return ' '.join([_fmt % x for x in arr])
                except:
                    return _fmt % arr
                    # except:
                    #     return numpy.array2string(arr, precision=self.precision, separator=',')[1:-1]
            numpy.set_string_function(array_fmt, repr=False)

    def _setup_index(self):
        """Sample indexing via tell / seek"""
        self._index_frame = []
        self._index_header = []
        self._index_cell = None
        self.trajectory.seek(0)
        while True:
            line = self.trajectory.tell()
            data = self.trajectory.readline().strip()

            # We break if file is over or we found an empty line
            if not data:
                break

            # The first line contains the number of particles
            # If something went wrong, this could be the last line
            # with the cell side (Lx,Ly,Lz) and we parse it some other way
            try:
                npart = int(data)
                self._index_header.append(line)
            except ValueError:
                self._index_cell = line
                break

            # Skip npart+1 lines
            _ = self.trajectory.readline()
            self._index_frame.append(self.trajectory.tell())
            for i in range(npart):
                _ = self.trajectory.readline()

    def _setup_steps(self):
        """Find steps list."""
        self.steps = []
        for frame in range(len(self._index_frame)):
            meta = self._read_metadata(frame)
            try:
                self.steps.append(meta['step'])
            except KeyError:
                # If no step info is found, we add steps sequentially
                self.steps.append(frame+1)

    def _expand_shortcuts(self):
        _fields = []
        for field in self.fields:
            try:
                _fields.append(self.shortcuts[field])
            except KeyError:
                _fields.append(field)
        return _fields

    def _read_metadata(self, frame):
        """
        Internal xyz method to get header metadata from comment line of
        given `frame`.

        We assume metadata fmt is a space-separated sequence of comma
        separated entries such as:

        columns:id,x,y,z step:10
        columns=id,x,y,z step=10
        """
        # Go to line and skip Npart info
        self.trajectory.seek(self._index_header[frame])
        npart = int(self.trajectory.readline())
        data = self.trajectory.readline()

        # Remove spaces around : or =
        data = re.sub(r'\W*[=:]\W*', ':', data)

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

        # Apply an alias dict to tags, e.g. to add step if Cfg was found instead
        for alias, tag in self.alias.items():
            try:
                meta[tag] = meta[alias]
            except KeyError:
                pass

        # Fix dimensions based on side of cell.
        # Fallback to ndim metadata or 3.
        try:
            if 'ndim' not in meta:
                meta['ndim'] = len(meta['cell'])
        except KeyError:
            meta['ndim'] = 3  # default

        return meta

    def read_init(self):
        # Grab cell from the end of file if it is there
        try:
            side = self._read_metadata(0)['cell']
            self._cell = Cell(side)
        except KeyError:
            self._cell = self._parse_cell()

    def read_sample(self, frame):
        meta = self._read_metadata(frame)
        if 'columns' in meta:
            # Use columns if they are found in the header
            fields = meta['columns']
            # Fix single column
            if not isinstance(fields, list):
                fields = [fields]
        else:
            # Stick to the default
            fields = self.fields
        fields = _optimize_fields(fields)

        # Read frame now
        self.trajectory.seek(self._index_frame[frame])
        particle = []
        for i in range(meta['npart']):
            p = Particle()
            # Note: we cannot optimize by shifting an index instead of
            # cropping lists all the time
            data = self.trajectory.readline().split()
            for key in fields:
                # If the key is associated to a explicit callback, go
                # for it. Otherwise we throw the field in an particle
                # attribute named key.
                if key in self.callback_read:
                    data = self.callback_read[key](p, data, meta)
                else:
                    p.__dict__[key] = tipify(data[0])
                    data = data[1:]
            particle.append(p)

        # Fix the masses.
        # We assume masses read from the header are sorted by species name.
        # The mass metadata must be adjusted to the given frame.
        if 'mass' in meta:
            species = distinct_species(particle)
            if len(species) == 1:
                for p in particle:
                    p.mass = float(meta['mass'])
            else:
                # We must have as many mass entries as species
                if len(species) != len(meta['mass']):
                    raise ValueError('mass metadata issue %s, %s' %
                                     (species, meta['mass']))
                db = {}
                for key, value in zip(species, meta['mass']):
                    db[key] = value
                for p in particle:
                    p.mass = float(db[p.species])

        # Check if we have a cell
        try:
            cell = Cell(meta['cell'])
        except KeyError:
            cell = self._cell

        return System(particle, cell)

    def read_timestep(self):
        meta = self._read_metadata(0)
        if 'dt' in meta:
            return meta['dt']
        elif 'timestep' in meta:
            return meta['timestep']
        else:
            return 1.0

    def _comment_header(self, step, system):
        # Comment line: concatenate metadata
        line = 'step:%d ' % step
        line += 'columns:' + ','.join(self.fields)
        if system.cell is not None:
            line += " cell:" + ','.join(['%s' % x for x in system.cell.side])
        try:
            line += " dt:%g" % self.timestep
        except:
            pass
        return line

    def write_sample(self, system, step):
        self._setup_format()
        if self._fields is None:
            self._fields = self._expand_shortcuts()
        self.trajectory.write('%d\n' % len(system.particle))
        self.trajectory.write(self._comment_header(step, system) + '\n')
        fmt = ' '.join(['{0.' + field + '}' for field in self._fields]) + '\n'
        for i, p in enumerate(system.particle):
            p._index = i
            p._step = step
            self.trajectory.write(fmt.format(p))

    def _parse_cell(self):
        """Internal xyz method to grab the cell. Can be overwritten in subclasses."""
        cell = None
        if self._index_cell:
            self.trajectory.seek(self._index_cell)
            side = numpy.fromstring(self.trajectory.readline(), sep=' ')
            cell = Cell(side)
        return cell

    def close(self):
        self.trajectory.close()


class TrajectoryNeighbors(TrajectoryXYZ):

    """Neighbors trajectory."""

    def __init__(self, filename, mode='r', offset=1):
        super(TrajectoryNeighbors, self).__init__(filename, mode=mode, alias={'time': 'step'})
        # TODO: determine minimum value of index automatically
        # TODO: possible regression here if no 'time' tag is found
        self._offset = offset  # neighbors produced by voronoi are indexed from 1
        self._netwon3 = False
        self._netwon3_message = False

    def read_sample(self, frame):
        meta = self._read_metadata(frame)
        self.trajectory.seek(self._index_frame[frame])
        s = System()
        s.neighbors = []
        for _ in range(meta['npart']):
            data = self.trajectory.readline().split()
            neigh = numpy.array(data, dtype=int)
            s.neighbors.append(neigh-self._offset)

        # Ensure III law Newton.
        # If this is ok on first frame we skip it for the next ones
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
