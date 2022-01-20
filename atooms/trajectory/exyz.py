import logging
import numpy

from .xyz import TrajectoryXYZ
from .utils import gopen
from atooms.core.utils import tipify
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System

log = logging.getLogger(__name__)


class TrajectoryEXYZ(TrajectoryXYZ):

    """
    Extended XYZ layout (https://github.com/libAtoms/extxyz)
    """

    suffix = 'xyz'
    callback_read = {}

    def __init__(self, filename, mode='r'):
        super(TrajectoryEXYZ, self).__init__(filename, mode)

        # Trajectory file handle
        self.precision = 12
        self._file = gopen(self.filename, self.mode)
        self.alias = {'pos': 'position',
                      'vel': 'velocity'}
        self.variables = ['species', 'pos']
        # Internal index of lines via seek and tell.
        if self.mode == 'r':
            self._setup_index()
            assert len(self._index_frame) > 0, 'empty file {}'.format(self._file)
            assert len(self._index_header) > 0, 'empty file {}'.format(self._file)
            # Read metadata
            self.metadata = self._read_comment(0)
            # Note: currently, changing variables in read mode will
            # not affect which attributes are parsed
            self.variables = [key for key, _, _ in self.metadata['Properties']]

    def _read_comment(self, frame):
        """
        Internal xyz method to get header metadata from comment line of
        given `frame`.
        """
        # Go to line and skip Npart info
        self._file.seek(self._index_header[frame])
        npart = int(self._file.readline())
        data = self._file.readline()
        meta = {}

        # We first gather all keys
        keys = []
        for entry in data.split("=")[:-1]:
            keys.append(entry.split()[-1])

        # Now we extract the values for each key using a dynamic regexp
        import re
        regexp = ''
        for key in keys:
            regexp += '{}=(.*)'.format(key)
        match = re.match(regexp, data)
        for i, key in enumerate(keys):
            meta[key] = match.groups()[i].strip()

        # The Properties key is special so we deal with it first
        entries = meta['Properties'].split(':')
        assert len(entries) % 3 == 0, len(entries)
        properties = []
        for i in range(0, len(entries), 3):
            key, fmt, ndims = entries[i: i+3]
            properties.append((key, fmt, ndims))

        meta['Properties'] = properties

        # Go through keys, listify and tipify them
        for key in meta:
            if key != 'Properties':
                if meta[key].startswith('"'):
                    meta[key] = meta[key].strip('"')
                    meta[key] = [tipify(_) for _ in meta[key].split()]
                else:
                    meta[key] = tipify(meta[key])

        return meta

    def read_steps(self):
        """Find steps list."""
        steps = []
        for frame in range(len(self._index_frame)):
            meta = self._read_comment(frame)
            try:
                steps.append(meta['Step'])
            except KeyError:
                # If no step info is found, we add steps sequentially
                steps.append(frame+1)
        return steps

    def read_system(self, frame):
        # Read metadata of this frame
        meta = self._read_comment(frame)

        # Get number of particles
        self._file.seek(self._index_header[frame])
        npart = int(self._file.readline())

        # Read frame now
        self._file.seek(self._index_frame[frame])
        particle = []
        for _ in range(npart):
            p = Particle()
            data = self._file.readline().split()
            i = 0
            for key, fmt, ndims in meta['Properties']:
                ndims = int(ndims)
                if key in self.alias:
                    key = self.alias[key]
                if ndims == 1:
                    if fmt == 'R':
                        setattr(p, key, float(data[i]))
                    elif fmt == 'I':
                        setattr(p, key, int(data[i]))
                    elif fmt == 'S':
                        setattr(p, key, data[i])
                    else:
                        raise ValueError('unknown format key')
                else:
                    if fmt == 'R':
                        setattr(p, key, numpy.array(data[i: i+ndims], dtype=numpy.float64))
                    elif fmt == 'I':
                        setattr(p, key, numpy.array(data[i: i+ndims], dtype=numpy.int64))
                    elif fmt == 'S':
                        setattr(p, key, numpy.array(data[i: i+ndims]))
                    else:
                        raise ValueError('unknown format key')
                i += ndims
            particle.append(p)

        side = meta["Lattice"]
        # TODO: remove hard coded
        cell = Cell([side[0], side[4], side[8]])
        return System(particle, cell)

    def read_timestep(self):
        for key in ['dt', 'Dt', 'timestep', 'Timestep']:
            if key in self.metadata:
                return self.metadata[key]
        return 1.0

    def _comment(self, step, system):
        # TODO: improve
        L = []
        for i in range(len(system.cell.side)):
            for j in range(len(system.cell.side)):
                if i == j:
                    L.append(system.cell.side[i])
                else:
                    L.append(0.0)

        # Reformat properties
        # TODO: improve
        _properties = ''
        for i in range(len(self._properties)):
            _properties += ':'.join([str(_) for _ in self._properties[i]]) + ':'
        _properties = _properties.strip(':')

        line = 'Properties={} '.format(_properties)
        line += 'Lattice="{}" '.format(' '.join([str(_) for _ in L]))
        line += 'Timestep={} '.format(self.timestep)
        line += 'Step={} '.format(step)
        return line.strip()

    def write_system(self, system, step):
        from atooms.core.utils import is_array

        def detect_format(arg):
            if isinstance(arg, float):
                return 'R'
            elif isinstance(arg, int):
                return 'I'
            else:
                return 'S'

        properties = []
        for variable in self.variables:
            attr = variable.split('.')[1]
            # Guess attribute features (ndim, fmt)
            val = getattr(system.particle[0], attr)
            if is_array(val):
                ndim, fmt = len(val), detect_format(val[0])
            else:
                ndim, fmt = 1, detect_format(val)
            properties.append([attr, fmt, ndim])
        # TODO: this should be aliased to stick with shortcuts
        self._properties = properties

        # Write header
        self._file.write('{}\n'.format(len(system.particle)))
        self._file.write(self._comment(step, system) + '\n')

        # Formatters
        alias_fmt = {'S': '{} ',
                     'R': '{{:.{precision}f}} '.format(precision=self.precision),
                     'I': '{} '}

        # Write particles
        for p in system.particle:
            line = ''
            for i in range(len(properties)):
                key, fmt, ndims = properties[i]
                val = getattr(p, key)
                if ndims == 1:
                    line += alias_fmt[fmt].format(val)
                else:
                    # Note: numpy.array2string is MUCH slower
                    line += ' '.join([alias_fmt[fmt].format(x) for x in val])
            self._file.write(line.strip() + '\n')
