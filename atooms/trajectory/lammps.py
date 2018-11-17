# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""LAMMPS trajectory format."""

import numpy

from .base import TrajectoryBase
from .folder import TrajectoryFolder
from atooms.system.particle import Particle, distinct_species
from atooms.system.cell import Cell
from atooms.system import System


# Formatting callbacks

def _parse_type(data, particle, cell):
    particle.species = data

def _parse_x(data, particle, cell):
    particle.position[0] = float(data)

def _parse_y(data, particle, cell):
    particle.position[1] = float(data)

def _parse_z(data, particle, cell):
    particle.position[2] = float(data)

def _parse_xs(data, particle, cell):
    particle.position[0] = (float(data) - 0.5) * cell.side[0]

def _parse_ys(data, particle, cell):
    particle.position[1] = (float(data) - 0.5) * cell.side[1]

def _parse_zs(data, particle, cell):
    particle.position[2] = (float(data) - 0.5) * cell.side[2]

def _parse_vx(data, particle, cell):
    particle.velocity[0] = float(data)

def _parse_vy(data, particle, cell):
    particle.velocity[1] = float(data)

def _parse_vz(data, particle, cell):
    particle.velocity[2] = float(data)


class TrajectoryLAMMPS(TrajectoryBase):

    """
    Trajectory layout for LAMMPS.

    In write mode, an additional .inp file is used as startup file.
    """

    suffix = 'atom'
    _cbk = {'x': _parse_x, 'y': _parse_y, 'z': _parse_z,
            'xu': _parse_x, 'yu': _parse_y, 'zu': _parse_z,
            'xs': _parse_xs, 'ys': _parse_ys, 'zs': _parse_zs,
            'xsu': _parse_xs, 'ysu': _parse_ys, 'zsu': _parse_zs,
            'vx': _parse_vx, 'vy': _parse_vy, 'vz': _parse_vz,
            'type': _parse_type}

    def __init__(self, filename, mode='r'):
        TrajectoryBase.__init__(self, filename, mode)
        self._fh = open(self.filename, self.mode)
        if mode == 'r':
            self._setup_index()

    def _setup_index(self):
        """Sample indexing via tell / seek"""
        from collections import defaultdict
        self._fh.seek(0)
        self._index_db = defaultdict(list)
        while True:
            line = self._fh.tell()
            data = self._fh.readline()
            # We break if file is over or we found an empty line
            if not data:
                break
            if data.startswith('ITEM:'):
                for block in ['TIMESTEP', 'NUMBER OF ATOMS',
                              'BOX BOUNDS', 'ATOMS']:
                    if data[6:].startswith(block):
                        # entry contains whatever is found after block
                        entry = data[7+len(block):]
                        self._index_db[block].append((line, entry))
                        break
        self._fh.seek(0)

    def read_steps(self):
        steps = []
        for idx, _ in self._index_db['TIMESTEP']:
            self._fh.seek(idx)
            self._fh.readline()
            step = int(self._fh.readline())
            steps.append(step)
        self._fh.seek(0)
        return steps

    def read_sample(self, frame):
        # Read number of particles
        idx, _ = self._index_db['NUMBER OF ATOMS'][frame]
        self._fh.seek(idx)
        self._fh.readline()
        data = self._fh.readline()
        npart = int(data)

        # Read box
        idx, data = self._index_db['BOX BOUNDS'][frame]
        self._fh.seek(idx)
        self._fh.readline()
        ndim = len(data.split())  # line is ITEM: BOX BONDS pp pp pp
        L, center = [], []
        for i in range(ndim):
            data = [float(x) for x in self._fh.readline().split()]
            L.append(data[1] - data[0])
            center.append((data[1] + data[0]) / 2)
        cell = Cell(numpy.array(L), center=numpy.array(center))

        # Read atoms data
        idx, data = self._index_db['ATOMS'][frame]
        fields = data.split()  # fields on a line
        _ = self._fh.readline()
        particles = [Particle() for i in range(npart)]
        for i in range(npart):
            data = self._fh.readline().split()
            # Accept unsorted particles by parsing their id
            if 'id' in fields:
                idx = int(data[0]) - 1
            else:
                idx = i
            # Read fields
            for j, field in enumerate(fields):
                if field in self._cbk:
                    self._cbk[field](data[j], particles[idx], cell)
                else:
                    # We should store these fields in particle anyway
                    pass

        return System(particle=particles, cell=cell)

    def write_init(self, system):
        f = open(self.filename + '.inp', 'w')
        np = len(system.particle)
        L = system.cell.side
        sp = distinct_species(system.particle)

        # LAMMPS header
        h = '\n'
        h += "%i atoms\n" % np
        h += "%i atom types\n" % len(sp)
        h += "%g %g  xlo xhi\n" % (-L[0]/2, L[0]/2)
        h += "%g %g  ylo yhi\n" % (-L[1]/2, L[1]/2)
        h += "%g %g  zlo zhi\n" % (-L[2]/2, L[2]/2)
        f.write(h + '\n')

        # LAMMPS body
        # Masses of species
        m = "\nMasses\n\n"
        for isp in range(len(sp)):
            # Iterate over particles. Find instances of species and get masses
            for p in system.particle:
                if p.species == sp[isp]:
                    m += '%s %g\n' % (isp+1, p.mass)
                    break

        # Atom coordinates
        r = "\nAtoms\n\n"
        v = "\nVelocities\n\n"
        for i, p in enumerate(system.particle):
            r += '%s %s %g %g %g\n' % tuple([i+1, sp.index(p.species)+1] + list(p.position))
            v += '%s    %g %g %g\n' % tuple([i+1] + list(p.velocity))

        f.write(m)
        f.write(r)
        f.write(v)
        f.close()

    def write_sample(self, system, step):
        pass

    def close(self):
        self._fh.close()


class TrajectoryFolderLAMMPS(TrajectoryFolder):

    """
    Trajectory layout for LAMMPS.

    In write mode, an additional .inp file is used as startup file.
    """

    suffix = '.tgz'

    def __init__(self, filename, mode='r', file_pattern='*', step_pattern=r'[a-zA-Z\.]*(\d*)'):
        TrajectoryFolder.__init__(self, filename, mode=mode,
                                  file_pattern=file_pattern,
                                  step_pattern=step_pattern)
        # We force reading steps from lammps file
        self._steps = None

    def read_steps(self):
        steps = []
        for filename in self.files:
            with TrajectoryLAMMPS(filename, 'r') as th:
                steps.append(th.steps[0])
        return steps

    def read_sample(self, frame):
        with TrajectoryLAMMPS(self.files[frame], 'r') as th:
            return th[0]

    def write_sample(self, system, step):
        # We cannot write
        raise NotImplementedError('cannot write lammps folder trajectory')

# Note: to get the tabulated potential from a dump of potential.x do
# > { echo -e "\nPOTENTIAL\nN 10000\n"; grep -v '#' /tmp/kalj.ff.potential.1-1 | \
#    awk '{printf "%i %g %12e %12e\n", NR, $1, $2, -$3}' ; }
