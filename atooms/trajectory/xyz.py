# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

import numpy
import re
import logging
from copy import copy
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

def _update_diameter(particle, data, meta):
    particle.radius = float(data[0]) / 2
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

def _update_neighbors_consume(particle, data, meta):
    # Consume all entries in data
    particle.neighbors = numpy.array(data, dtype=int)
    return []

def _update_neighbors(particle, data, meta):
    # Extract comma separated entries in the first element of data
    if len(data) > 0 and ',' in data[0]:
        if data[0] == ',':
            particle.neighbors = numpy.array([], dtype=int)
        else:
            particle.neighbors = numpy.array(data[0].split(','), dtype=int)
        return data[1:]
    else:
        particle.neighbors = numpy.array(data, dtype=int)
        return []

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
                     'diameter': _update_diameter,
                     'pos': _update_position,
                     'vel': _update_velocity,
                     'position': _update_position,
                     'velocity': _update_velocity,
                     'neighbors': _update_neighbors,
                     'neighbors*': _update_neighbors_consume}

    def __init__(self, filename, mode='r', alias=None, fields=None):
        super(TrajectoryXYZ, self).__init__(filename, mode)
        self._fields_default = ['id', 'pos']
        self.fields = copy(self._fields_default) if fields is None else fields
        self.alias = {} if alias is None else alias

        # Internal variables
        self._cell = None
        self._fields_float = True
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
        # Note: there is little gain in using the read_len() method to just
        # go through the file instead of setting up the index. So we stick to it.
        if self.mode == 'r':
            self._setup_index()
            assert len(self._index_frame) > 0, 'empty file {}'.format(self.trajectory)
            assert len(self._index_header) > 0, 'empty file {}'.format(self.trajectory)
            # Read metadata
            self.metadata = self._read_comment(0)
            # Redefine fields
            self._setup_fields()

    def _setup_fields(self):
        """
        Redefine fields.

        If set by the user, self.fields takes precedence on columns metadata.
        - fields is default and columns is present: read everything
        - fields is default and columns not present: read from fields
        - fields is set by user and columns is present: pad
        - fields is set by user but columns not present: read from fields
        """
        if self.fields == self._fields_default:
            # Default fields: metadata has priority
            if 'columns' in self.metadata:
                self.fields = self.metadata['columns']
        else:
            # Fields set by user: pad missing fields if columns metadata are present
            if 'columns' in self.metadata:
                fields = []
                for key in self.metadata['columns']:
                    if key in self.fields + _optimize_fields(self.fields):
                        fields.append(key)
                    else:
                        fields.append(None)
                self.fields = fields

        # Replace some column keys with more efficient ones
        self.fields = _optimize_fields(self.fields)

        # Remove skippable fields
        for key in self.fields[-1::-1]:
            if key is None:
                self.fields.pop()
            else:
                break

        # Store a copy of fields in a cache variable
        # We can avoid setup if the fields have not changed
        self.__cache_fields = copy(self.fields)

    def _setup_format(self):
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
                # Note: numpy.array2string is MUCH slower
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
            # TODO: drop compatibility with xyz files having cell as last line.
            try:
                npart = int(data)
                self._index_header.append(line)
            except ValueError:
                self._index_cell = line
                break

            # Skip npart+1 lines, making sure we have read precisely
            # that number of lines
            _ = self.trajectory.readline()
            line = self.trajectory.tell()
            for i in range(npart):
                _ = self.trajectory.readline()
            # Store first line /after/ we have read the frame
            # making sure the last we read was not emtpy
            # Note that readline() returns an empty string on EOF
            if len(_) > 0:
                self._index_frame.append(line)
            else:
                raise IOError('malformed xyz file [%s]', self.filename)

    def read_steps(self):
        """Find steps list."""
        steps = []
        for frame in range(len(self._index_frame)):
            meta = self._read_comment(frame)
            try:
                steps.append(meta['step'])
            except KeyError:
                # If no step info is found, we add steps sequentially
                steps.append(frame+1)
        return steps

    def _read_comment(self, frame):
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
        meta = {}
        meta['npart'] = npart

        # Read metadata line
        if not ('=' in data or ':' in data):
            # The comment line contains a list of unspecified fields
            # We add them to the dict using a sequential integer key
            for i, value in enumerate(data.split()):
                meta[i] = tipify(value)

        else:
            # The comment line contains self descriptive fields
            # TODO: accept extended xyz format
            # Remove spaces around : or = and replace = by :
            data = re.sub(r'\s*[=:]\s*', ':', data)

            # Fill metadata dictionary
            for e in data.split():
                s = re.search(r'(\S+):(\S+)', e)
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
        except TypeError:
            meta['ndim'] = 1  # it is a scalar
        except KeyError:
            meta['ndim'] = 3  # default

        # Make sure columns is a list
        if 'columns' in meta:
            if not isinstance(meta['columns'], list):
                meta['columns'] = [meta['columns']]

        return meta

    def read_sample(self, frame):
        # Setup fields again, in case they have changed
        if self.__cache_fields != self.fields:
            self._setup_fields()
        # Read metadata of this frame
        meta = self._read_comment(frame)
        # Define read callbacks list before reading lines
        callbacks_read = []
        def _skip(p, data, meta):
            return data[1:]
        for key in self.fields:
            # If the key is associated to a explicit callback, go
            # for it. Otherwise we throw the field in an particle
            # attribute named key. If the key is None it means we skip
            # this column.
            if key is None:
                callbacks_read.append(_skip)
            elif key in self.callback_read:
                callbacks_read.append(self.callback_read[key])
            else:
                # Trick. We instantiate dynamically a fallback function
                # to avoid adding `key` to the other callbacks' interface
                namespace = {}
                exec("""
from atooms.core.utils import tipify
def fallback(p, data, meta):
    p.__dict__['%s'] = tipify(data[0])
    return data[1:]
""" % key, namespace)
                callbacks_read.append(namespace['fallback'])

        # Read frame now
        self.trajectory.seek(self._index_frame[frame])
        particle = []
        for i in range(meta['npart']):
            p = Particle()
            # Note: we cannot optimize by shifting an index instead of
            # cropping lists all the time
            data = self.trajectory.readline().split()
            for cbk in callbacks_read:
                data = cbk(p, data, meta)
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

        # Add cell info
        if 'cell' in meta:
            cell = Cell(meta['cell'])
        else:
            cell = None

        return System(particle, cell)

    def read_timestep(self):
        if 'dt' in self.metadata:
            return self.metadata['dt']
        elif 'timestep' in self.metadata:
            return self.metadata['timestep']
        else:
            return 1.0

    def _comment(self, step, system):
        # Concatenate metadata in comment line
        line = 'step:{} '.format(step)
        line += 'columns:{} '.format(','.join(self.fields))
        if self.timestep is not None:
            line += 'dt:{:g} '.format(self.timestep)
        if system.cell is not None:
            line += 'cell:{} '.format(','.join([str(x) for x in system.cell.side]))
        for entry in self.metadata:
            line += '{}:{} '.format(entry, self.metadata[entry])
        return line

    def write_sample(self, system, step):
        # Make sure fields are expanded
        self._setup_format()
        self.trajectory.write('%d\n' % len(system.particle))
        self.trajectory.write(self._comment(step, system) + '\n')
        # Expand shortcut fields now (the comment header keeps the shortcuts)
        fields = _expand_fields(self.fields, self._shortcuts)
        fmt = ' '.join(['{0.' + field + '}' for field in fields]) + '\n'
        for i, p in enumerate(system.particle):
            p._index = i
            p._step = step
            self.trajectory.write(fmt.format(p))

    def close(self):
        self.trajectory.close()
        # Restore numpy formatting defaults
        numpy.set_string_function(None, False)
        numpy.set_string_function(None, True)


def _add_neighbors_to_system(system, offset):
    for p in system.particle:
        p.neighbors -= offset
    system.neighbors = [p.neighbors for p in system.particle]
    return system


class TrajectoryNeighbors(TrajectoryXYZ):

    """
    Neighbors trajectory.

    By default, for reading we assume an xyz file with space-separated
    integers indicating the particles indices. For writing, we use
    comma-separated entries (no space).
    """

    _fields_default = ['neighbors*']

    def __init__(self, filename, mode='r', fields=None, offset=1):
        super(TrajectoryNeighbors, self).__init__(filename, mode=mode,
                                                  alias={'time': 'step'})
        self._offset = offset
        if mode == 'w':
            self.fields = ['neighbors'] if fields is None else fields
        self.callback_read['neighbors'] = _update_neighbors
        self.callback_read['neighbors*'] = _update_neighbors_consume
        self.add_callback(_add_neighbors_to_system, self._offset)
        self._fields_float = False
