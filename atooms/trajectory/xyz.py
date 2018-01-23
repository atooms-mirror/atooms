# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

import numpy
import re
import logging
from .simple import TrajectorySimpleXYZ
from .base import TrajectoryBase
from .utils import gopen
from atooms.core.utils import tipify
from atooms.system.particle import Particle, distinct_species
from atooms.system.cell import Cell
from atooms.system import System

log = logging.getLogger(__name__)


# Format callbacks

def _update_radius(particle, data, meta):
    particle.radius = float(data[0])
    return data[1:]

def _update_tag(particle, data, meta):
    # Kept for backward compatibility. We should rather rely on
    # tipified properties.
    particle.tag = data[0:]
    return data[1:]

def _update_species(particle, data, meta):
    particle.species = data[0]
    return data[1:]

def _update_position(particle, data, meta):
    ndim = meta['ndim']
    # It is crucial to assing position, not to use the slice!
    # Otherwise we get a reference, not a copy.
    particle.position = numpy.array(data[0:ndim], dtype=float)
    return data[ndim:]

def _update_velocity(particle, data, meta):
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

def _expand_fields(fields, shortcuts):
    """Expand `shortcuts` present in `fields`"""
    _fields = []
    for field in fields:
        try:
            _fields.append(shortcuts[field])
        except KeyError:
            _fields.append(field)
    return _fields


class TrajectoryXYZ(TrajectoryBase):

    """
    Trajectory with XYZ layout using memory leightweight indexed
    access.
    """

    suffix = 'xyz'
    callback_read = {'species': _update_species,
                     'type': _update_species,  # alias
                     'name': _update_species,  # alias
                     'id': _update_species,  # alias
                     'tag': _update_tag,
                     'radius': _update_radius,
                     'pos': _update_position,
                     'vel': _update_velocity,
                     'position': _update_position,
                     'velocity': _update_velocity}

    def __init__(self, filename, mode='r', alias=None, fields=None):
        TrajectoryBase.__init__(self, filename, mode)
        self.alias = {} if alias is None else alias
        self.fields = ['id', 'pos'] if fields is None else fields

        # Internal variables
        self._cell = None
        self._fields_float = True
        self._done_format_setup = False
        self._shortcuts = {'pos': 'position',
                           'x': 'position[0]',
                           'y': 'position[1]',
                           'z': 'position[2]',
                           'vel': 'velocity',
                           'vx': 'velocity[0]',
                           'vy': 'velocity[1]',
                           'vz': 'velocity[2]',
                           'id': 'species',
                           'type': 'species'}

        # Trajectory file handle
        self.trajectory = gopen(self.filename, self.mode)

        # Internal index of lines via seek and tell.
        if self.mode == 'r':
            self._setup_index()

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
            data = self.trajectory.readline()

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

    def read_steps(self):
        """Find steps list."""
        steps = []
        for frame in range(len(self._index_frame)):
            meta = self._read_metadata(frame)
            try:
                steps.append(meta['step'])
            except KeyError:
                # If no step info is found, we add steps sequentially
                steps.append(frame+1)
        return steps

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
        # Redefine fields
        if 'columns' in meta:
            # Use columns as fields if they are found in the header
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
                # This block below takes ~50% of the time
                # --
                # If the key is associated to a explicit callback, go
                # for it. Otherwise we throw the field in an particle
                # attribute named key.
                if key in self.callback_read:
                    data = self.callback_read[key](p, data, meta)
                else:
                    p.__dict__[key] = tipify(data[0])
                    data = data[1:]
                # --
            particle.append(p)

        # Fix the masses.
        # We assume masses read from the header are sorted by species name.
        # The mass metadata must be adjusted to the given frame.
        if 'mass' in meta:
            if isinstance(meta['mass'], list) or isinstance(meta['mass'], tuple):
                species = distinct_species(particle)
                # We must have as many mass entries as species
                if len(species) != len(meta['mass']):
                    raise ValueError('mass metadata issue %s, %s' %
                                     (species, meta['mass']))
                db = {}
                for key, value in zip(species, meta['mass']):
                    db[key] = value
                for p in particle:
                    p.mass = float(db[p.species])
            else:
                for p in particle:
                    p.mass = float(meta['mass'])

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
        # Concatenate metadata in comment line
        line = 'step:{} '.format(step)
        line += 'columns:{} '.format(','.join(self.fields))
        if self.timestep is not None:
            line += 'dt:{:g} '.format(self.timestep)
        if system.cell is not None:
            line += 'cell:{} '.format(','.join([str(x) for x in system.cell.side]))
        return line

    def write_sample(self, system, step):
        # Make sure fields are expanded
        self._setup_format()
        self.trajectory.write('%d\n' % len(system.particle))
        self.trajectory.write(self._comment_header(step, system) + '\n')
        # Expand shortcut fields now (the comment header keeps the shortcuts)
        fields = _expand_fields(self.fields, self._shortcuts)
        fmt = ' '.join(['{0.' + field + '}' for field in fields]) + '\n'
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
        # Restore numpy formatting defaults
        numpy.set_string_function(None, False)
        numpy.set_string_function(None, True)


def _update_neighbors_consume(particle, data, meta):
    # Consume all entries in data
    particle.neighbors = numpy.array(data, dtype=int)
    return []

def _update_neighbors(particle, data, meta):
    # Extract comma separated entries in the first element of data
    particle.neighbors = numpy.array(data[0].split(','), dtype=int)
    return data[1:]

def _add_neighbors_to_system(system, offset):
    for p in system.particle:
        p.neighbors -= offset
    system.neighbors = [p.neighbors for p in system.particle]
    return system


class TrajectoryNeighbors(TrajectoryXYZ):

    """Neighbors trajectory."""

    def __init__(self, filename, mode='r', fields=None, offset=1):
        super(TrajectoryNeighbors, self).__init__(filename, mode=mode,
                                                  alias={'time':
                                                         'step'})
        self._offset = offset
        self.fields = ['neighbors*'] if fields is None else fields
        self.callback_read = {'neighbors': _update_neighbors,
                              'neighbors*': _update_neighbors_consume}
        self.add_callback(_add_neighbors_to_system, self._offset)
