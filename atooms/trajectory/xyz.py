# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

import warnings
import numpy
import re
import logging
from copy import copy
from .simple import TrajectorySimpleXYZ
from .base import TrajectoryBase, canonicalize
from .utils import gopen
from atooms.core.utils import tipify
from atooms.system.particle import Particle, distinct_species
from atooms.system.cell import Cell
from atooms.system import System

log = logging.getLogger(__name__)

# Lock to avoid race conditions on numpy array formatting 
_numpy_fmt_lock = False

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

def _update_position_unfolded(particle, data, meta):
    from atooms.system.particle import _periodic_vector_unfolded
    ndim = meta['ndim']
    r = numpy.array(data[0:ndim], dtype=float)
    particle.position_unfolded = r
    particle.position = _periodic_vector_unfolded(r, meta['cell'])
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

def _optimize_variables(variables):
    # TODO: this should be more clever, say we had a 2d simulation
    #       embedded in 3d we want only to read x and y but this will
    #       force to read a full vector and it may even cause problems
    new_variables = []
    for variable in variables:
        if variable == 'x':
            new_variables.append('particle.position')
        elif variable == 'vx':
            new_variables.append('particle.velocity')
        elif variable in ['y', 'z', 'vy', 'vz']:
            pass
        else:
            new_variables.append(variable)
    return new_variables


class TrajectoryXYZ(TrajectoryBase):

    """
    Trajectory with XYZ layout using memory leightweight indexed
    access.
    """

    suffix = 'xyz'
    callback_read = {'species': _update_species,
                     'particle.species': _update_species,
                     'type': _update_species,
                     'name': _update_species,
                     'id': _update_species,
                     'tag': _update_tag,
                     'radius': _update_radius,
                     'diameter': _update_diameter,
                     'pos': _update_position,
                     'vel': _update_velocity,
                     'position': _update_position,
                     'position_unfolded': _update_position_unfolded,
                     'particle.position': _update_position,
                     'particle.position_unfolded': _update_position_unfolded,
                     'velocity': _update_velocity,
                     'particle.velocity': _update_velocity,
                     'neighbors': _update_neighbors,
                     'neighbors*': _update_neighbors_consume}

    def __init__(self, filename, mode='r', alias=None, fields=None):        
        super(TrajectoryXYZ, self).__init__(filename, mode)
        self.variables = [
            'particle.species',
            'particle.position'
        ]
        self.constants = {
            'particle.mass',
            'cell.side'
        }
        if fields is not None:
            self.variables = fields
            warnings.warn('fields is deprecated, set trajectory.variables instead', FutureWarning)

        self.alias = {} if alias is None else alias
        # Internal variables
        self._cell = None
        # Trajectory file handle
        self.trajectory = gopen(self.filename, self.mode)

        if self.mode == 'r':
            # Internal index of lines via seek and tell.
            self._setup_index()
            # Read metadata
            self.metadata = self._read_comment(0)
            # Update schema
            self._setup_schema()
            self._active_read_callbacks = None

    def _setup_format(self):
        _fmt_any = '{}'
        _fmt_float = '{:.' + str(self.precision) + 'f}'
            
        def array_fmt(arr):
            """Remove commas and [] from numpy array repr."""
            # Passing a scalar will trigger an error (gotcha: even
            # when casting numpy array to list, the elements remain of
            # numpy type and this function gets called! (4% slowdown)
            if arr.dtype == 'float':
                _fmt = _fmt_float
            else:
                _fmt = _fmt_any
            if len(arr.shape) == 1:
                return ' '.join([_fmt.format(x) for x in arr])
            elif len(arr.shape) == 0:
                return _fmt.format(arr)
            else:
                raise ValueError('cannot write multi-dimensional array')
                
        # Note: numpy.array2string is MUCH slower
        numpy.set_string_function(array_fmt, repr=False)

    def _setup_index(self):
        """Sample indexing via tell / seek"""
        self._index_frame = []
        self._index_header = []
        self.trajectory.seek(0)
        while True:
            line = self.trajectory.tell()
            data = self.trajectory.readline()

            # We break if file is over or we found an empty line
            if not data:
                break

            # The first line contains the number of particles
            # If something goes wrong, this could be the last line
            # with the cell side (Lx,Ly,Lz). We do not support this anymore
            try:
                npart = int(data)
                self._index_header.append(line)
            except ValueError:
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

        # Sanity tests
        assert len(self._index_frame) > 0, 'empty file {}'.format(self.trajectory)
        assert len(self._index_header) > 0, 'empty file {}'.format(self.trajectory)

    def _setup_schema(self):
        """
        Update variables and constants from trajectory metadata
        """
        # Default fields: metadata has priority
        if 'columns' in self.metadata:
            self.variables = self.metadata['columns']
        # Replace some column keys with more efficient ones
        # Some subclasses may not want to do that
        self.variables = _optimize_variables(self.variables)
        # Store a copy of fields in a cache variable
        # This way we can avoid setup if the fields have not changed
        self._original_variables = copy(self.variables)

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
            # Remove spaces around : or = and replace = by :
            data = re.sub(r'\s*[=:]\s*', ':', data)
            # Tolerate spaces between entries of vectors
            data = re.sub(r'\s*,\s*', ',', data)

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

    def _setup_read_callbacks(self):
        """
        Return a list of callbacks to read the variables given in self.variables
        """
        # We only call this method once, when the we read the first sample
        # If self.variables change after setup, we ignore the change
        # Similarly, we assume that columns do not change
        if self._active_read_callbacks is not None:
            return self._active_read_callbacks

        # Check that we have not reset the variables
        if self.variables is None:
            self.variables = self._original_variables
        
        # Define active callbacks before reading the sample.
        # We handle the case in which the user has removed some
        # variables, and self.variables is a subset of
        # self._original_variables
        callbacks = []
        test_particle = Particle()
        for key in self._original_variables:
            if key not in self.variables:
                # If the key is not the in original list of variables
                # the user has removed it from the list and we skip it
                nsafe = 20
                cbk = self.callback_read[key]
                consumed_entries = nsafe - len(cbk(test_particle, range(nsafe), self.metadata))
                # Trick. We instantiate dynamically a skip function
                # with the right number of consumed entries
                namespace = {}
                exec("""
def skip(p, data, meta):
    return data[{}:]
""".format(consumed_entries), namespace)
                callbacks.append(namespace['skip'])
            elif key in self.callback_read:
                # The key is associated to a explicit callback
                callbacks.append(self.callback_read[key])
            else:
                # We throw the field in an particle.
                # Trick. We instantiate dynamically a fallback function
                # to avoid adding `key` to the other callbacks' interface
                namespace = {}
                exec("""
from atooms.core.utils import tipify
def fallback(p, data, meta):
    p.__dict__['%s'] = tipify(data[0])
    return data[1:]
""" % key, namespace)
                callbacks.append(namespace['fallback'])

        self._active_read_callbacks = callbacks
                
    def read_sample(self, frame):
        # Define actual list of callbacks
        self._setup_read_callbacks()
        callbacks = self._active_read_callbacks
        
        # Read metadata of this frame
        db = self._read_comment(frame)

        # Read frame now
        self.trajectory.seek(self._index_frame[frame])
        particle = []
        for i in range(db['npart']):
            p = Particle()
            # Note: we cannot optimize by shifting an index instead of
            # cropping lists all the time
            data = self.trajectory.readline().split()
            for cbk in callbacks:
                data = cbk(p, data, db)
            particle.append(p)

        # Fix the masses.
        # We assume masses read from the header are sorted by species name.
        # The mass dbdata must be adjusted to the given frame.
        if 'mass' in db:
            if isinstance(db['mass'], list) or isinstance(db['mass'], tuple):
                species = distinct_species(particle)
                # We must have as many mass entries as species
                if len(species) != len(db['mass']):
                    raise ValueError('mass metadata issue %s, %s' %
                                     (species, db['mass']))
                db_species = {}
                for key, value in zip(species, db['mass']):
                    db_species[key] = value
                for p in particle:
                    p.mass = float(db_species[p.species])
            else:
                for p in particle:
                    p.mass = float(db['mass'])

        # Add cell info
        if 'cell' in db:
            cell = Cell(db['cell'])
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
        line += 'columns:{} '.format(','.join(self.variables))
        if self.timestep is not None:
            line += 'dt:{:g} '.format(self.timestep)
        if system.cell is not None:
            line += 'cell:{} '.format(','.join([str(x) for x in system.cell.side]))
        for entry in self.metadata:
            line += '{}:{} '.format(entry, self.metadata[entry])
        return line

    def write_sample(self, system, step):
        # Make sure fields are expanded
        global _numpy_fmt_lock
        _numpy_fmt_lock = True
        self._setup_format()
        self.trajectory.write('%d\n' % len(system.particle))
        self.trajectory.write(self._comment(step, system) + '\n')
        # Canonicalize variables now.
        # The comment header keeps the non-canonical fields
        variables = canonicalize(self.variables)
        fmt = ' '.join(['{0.' + variable.split('particle.')[-1] + '}' for variable in variables]) + '\n'
        for i, p in enumerate(system.particle):
            p._index = i
            p._step = step
            self.trajectory.write(fmt.format(p))
        _numpy_fmt_lock = False

    def close(self):
        self.trajectory.close()
        # Restore numpy formatting defaults
        if not _numpy_fmt_lock:
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
