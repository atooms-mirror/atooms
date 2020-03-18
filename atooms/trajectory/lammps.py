# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""LAMMPS trajectory format."""

import numpy

from .base import TrajectoryBase
from .folder import TrajectoryFolder
from atooms.system import System, Particle, Cell
from atooms.system.particle import distinct_species
from atooms.interaction import Interaction

import sys

# Redefine range for python
if sys.version_info[0] == 2:
    range = xrange

# Formatting callbacks

def _parse_type(data, idx, system):
    system.particle[idx].species = data

def _parse_x(data, idx, system):
    system.particle[idx].position[0] = float(data)

def _parse_y(data, idx, system):
    system.particle[idx].position[1] = float(data)

def _parse_z(data, idx, system):
    system.particle[idx].position[2] = float(data)

def _parse_xs(data, idx, system):
    system.particle[idx].position[0] = (float(data) - 0.5) * system.cell.side[0]

def _parse_ys(data, idx, system):
    system.particle[idx].position[1] = (float(data) - 0.5) * system.cell.side[1]

def _parse_zs(data, idx, system):
    system.particle[idx].position[2] = (float(data) - 0.5) * system.cell.side[2]

def _parse_vx(data, idx, system):
    system.particle[idx].velocity[0] = float(data)

def _parse_vy(data, idx, system):
    system.particle[idx].velocity[1] = float(data)

def _parse_vz(data, idx, system):
    system.particle[idx].velocity[2] = float(data)

def _parse_fx(data, idx, system):
    system.interaction.forces[idx, 0] = float(data)

def _parse_fy(data, idx, system):
    system.interaction.forces[idx, 1] = float(data)

def _parse_fz(data, idx, system):
    system.interaction.forces[idx, 2] = float(data)


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
            'fx': _parse_fx, 'fy': _parse_fy, 'fz': _parse_fz,
            'type': _parse_type}

    def __init__(self, filename, mode='r', single_frame=False,
                 first_particle=-1, last_particle=-1):
        TrajectoryBase.__init__(self, filename, mode)
        self.precision = 14  # default to double precision
        self.single_frame = single_frame
        self.first_particle = first_particle
        self.last_particle = last_particle
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
                # Avoid reading after ATOOMS block has been found.
                # We assume it is the last block in the file.
                # The single_frame variable is a hint that there are
                # no more frames in the file
                if block == 'ATOMS' and self.single_frame:
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

        # Build the system
        system = System()
        system.particle = []
        for i in range(npart):
            if self.first_particle > 0 and i < self.first_particle:
                continue
            if self.last_particle > 0 and i >= self.last_particle:
                break
            system.particle.append(Particle())

        # Add cell
        idx, data = self._index_db['BOX BOUNDS'][frame]
        self._fh.seek(idx)
        self._fh.readline()
        ndim = len(data.split())  # line is ITEM: BOX BONDS pp pp pp
        L, center = [], []
        for i in range(ndim):
            data = [float(x) for x in self._fh.readline().split()]
            L.append(data[1] - data[0])
            center.append((data[1] + data[0]) / 2)
        system.cell = Cell(numpy.array(L), center=numpy.array(center))

        # Read atoms data
        idx, data = self._index_db['ATOMS'][frame]
        fields = data.split()  # fields on a line
        _ = self._fh.readline()

        # Add interaction if forces are present
        # In atooms, forces belong to the interaction, not to particles
        if 'fx' in fields or 'fy' in fields or 'fz' in fields:
            # TODO: this won't work with first and last particles
            system.interaction = Interaction([])  # empty list of potentials
            system.interaction.forces = numpy.ndarray((npart, ndim))
        else:
            interaction = None

        for i in range(npart):
            # Limit reading the ATOMS section if requested
            if self.first_particle > 0 and i < self.first_particle:
                continue
            if self.last_particle > 0 and i >= self.last_particle:
                break
            data = self._fh.readline().split()
            # Accept unsorted particles by parsing their id
            if 'id' in fields:
                idx = int(data[0]) - 1
            else:
                idx = i
            # Populate particle's attributes by reading fields
            for j, field in enumerate(fields):
                if field in self._cbk:
                    self._cbk[field](data[j], idx, system)
                else:
                    # We should store these fields in particle anyway
                    pass

        return system

    def write_init(self, system):
        f = open(self.filename + '.inp', 'w')
        np = len(system.particle)
        L = system.cell.side
        species_db = distinct_species(system.particle)

        # LAMMPS header
        h = '\n'
        h += "{:d} atoms\n".format(np)
        h += "{:d} atom types\n".format(len(species_db))
        h += "{:.{prec}f} {:.{prec}f} xlo xhi\n".format(-L[0]/2, L[0]/2, prec=self.precision)
        h += "{:.{prec}f} {:.{prec}f} ylo yhi\n".format(-L[1]/2, L[1]/2, prec=self.precision)
        h += "{:.{prec}f} {:.{prec}f} zlo zhi\n".format(-L[2]/2, L[2]/2, prec=self.precision)

        # LAMMPS body
        # Masses of species
        m = "\nMasses\n\n"
        for isp in range(len(species_db)):
            # Iterate over particles. Find instances of species and get masses
            for p in system.particle:
                if p.species == species_db[isp]:
                    m += '{:d} {:{prec}f}\n'.format(isp + 1, p.mass, prec=self.precision)
                    break

        # Atom coordinates
        r = "\nAtoms\n\n"
        v = "\nVelocities\n\n"
        for i, p in enumerate(system.particle):
            isp = species_db.index(p.species) + 1
            r += '{:d} {:d} {:{prec}} {:{prec}} {:{prec}}\n'.format(i+1, isp, *p.position, prec=self.precision)
            v += '{:d}      {:{prec}} {:{prec}} {:{prec}}\n'.format(i+1, *p.velocity, prec=self.precision)

        f.write(h)
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

    suffix = 'tgz'

    def __init__(self, filename, mode='r', file_pattern='*',
                 step_pattern=r'[a-zA-Z\.]*(\d*)', first_particle=-1, last_particle=-1):
        TrajectoryFolder.__init__(self, filename, mode=mode,
                                  file_pattern=file_pattern,
                                  step_pattern=step_pattern)
        self.first_particle = first_particle
        self.last_particle = last_particle
        # Small trick to force reading steps from lammps file
        self._steps = None
        # Sort frames according to step read in lammps file
        sorted_steps = sorted(self.steps)
        files_with_steps = [(x, y) for x, y in zip(self.files, self.steps)]
        files_with_steps.sort(key=lambda x: sorted_steps.index(x[1]))
        files = [_[0] for _ in files_with_steps]
        self.files = files
        self._steps = sorted_steps

    def read_steps(self):
        steps = []
        for filename in self.files:
            with TrajectoryLAMMPS(filename, 'r', single_frame=True,
                                  first_particle=self.first_particle,
                                  last_particle=self.last_particle) as th:
                steps.append(th.steps[0])
        return steps

    def read_sample(self, frame):
        with TrajectoryLAMMPS(self.files[frame], 'r', single_frame=True,
                              first_particle=self.first_particle,
                              last_particle=self.last_particle) as th:
            return th[0]

    def write_sample(self, system, step):
        # We cannot write
        raise NotImplementedError('cannot write lammps folder trajectory')

# Note: to get the tabulated potential from a dump of potential.x do
# > { echo -e "\nPOTENTIAL\nN 10000\n"; grep -v '#' /tmp/kalj.ff.potential.1-1 | \
#    awk '{printf "%i %g %12e %12e\n", NR, $1, $2, -$3}' ; }
