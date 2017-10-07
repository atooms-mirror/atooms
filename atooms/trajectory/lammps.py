# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""LAMMPS trajectory format."""

import numpy

from .base import TrajectoryBase
from .folder import TrajectoryFolder
from atooms.system.particle import Particle, distinct_species
from atooms.system.cell import Cell
from atooms.system import System


def _read_item(t, item):
    """Parse generic `item`."""
    data = t.readline()
    if item not in data:
        raise ValueError('expecting "%s" got "%s" on %s' % (item, data, t.name))
    return data

def _lammps_parse_step(finp):
    """Read step from `finp`."""
    data = _read_item(finp, 'ITEM: TIMESTEP')
    data = finp.readline()
    return int(float(data))

def _lammps_parse_system(finp):
    """Read system from `finp`."""
    # Read number of particles
    data = _read_item(finp, 'ITEM: NUMBER OF ATOMS')
    data = finp.readline()
    npart = int(data)

    # Read box
    data = _read_item(finp, "ITEM: BOX BOUNDS")
    ndim = len(data.split()[3:])  # line is ITEM: BOX BONDS pp pp pp
    if len(data.split()) == 3:
        ndim = 3
    L, offset = [], []
    for i in range(ndim):
        data = [float(x) for x in finp.readline().split()]
        L.append(data[1] - data[0])
        offset.append(data[0])
    c = Cell(numpy.array(L))

    # Read positions and velocities
    data = _read_item(finp, "ITEM: ATOMS").split()[2:]
    # Determine how many variables are there
    n = len(data)
    ind = data.index('type')
    # Index of x coordinate
    try:
        ix = data.index('x')
    except ValueError:
        ix = data.index('xu')
    # Index of vx coordinate, if there at all
    try:
        ivx = data.index('vx')
    except ValueError:
        ivx = None

    # Rather well optimized now. Ignore velocities, do not fold particles.
    from itertools import islice
    data = ''.join(list(islice(finp, npart)))
    d = numpy.fromstring(data, sep=' ').reshape((npart, n))
    if ivx is not None:
        p = [Particle(species=d[i, ind], mass=1.0,
                      position=d[i, ix: ix+3], velocity=d[i, ivx: ivx+3])
             for i in xrange(npart)]

    else:
        p = [Particle(species=d[i, ind], mass=1.0,
                      position=d[i, ix: ix+3]) for i in xrange(npart)]

    for pi in p:
        pi.fold(c)

    p.sort(key=lambda a: a.species)
    return System(particle=p, cell=c)


class TrajectoryLAMMPS(TrajectoryBase):

    """
    Trajectory layout for LAMMPS.

    In write mode, an additional .inp file is used as startup file.
    """

    def __init__(self, filename, mode='r'):
        TrajectoryBase.__init__(self, filename, mode)
        self.steps = [0]

    def read_sample(self, frame):
        # TODO: respect input frame
        with open(self.filename, 'r') as fh:
            _ = _lammps_parse_step(fh)
            s = _lammps_parse_system(fh)
        return s

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


class TrajectoryFolderLAMMPS(TrajectoryFolder):

    """
    Trajectory layout for LAMMPS.

    In write mode, an additional .inp file is used as startup file.
    """

    suffix = '.tgz'

    def _get_step(self, fileinp):
        with open(fileinp) as fh:
            step = _lammps_parse_step(fh)
        return step

    def read_sample(self, frame):
        with open(self.files[frame], 'r') as fh:
            _ = _lammps_parse_step(fh)
            s = _lammps_parse_system(fh)
        return s

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
        # We cannot write
        return

# Note: to get the tabulated potential from a dump of potential.x do
# > { echo -e "\nPOTENTIAL\nN 10000\n"; grep -v '#' /tmp/kalj.ff.potential.1-1 | \
#    awk '{printf "%i %g %12e %12e\n", NR, $1, $2, -$3}' ; }
