# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""LAMMPS trajectory."""

import os
import re
import tarfile
import numpy

from atooms.trajectory import TrajectoryBase
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System

def _read_item(t, item):
    data = t.readline()
    if not item in data:
        raise ValueError('expecting "%s" got "%s" on %s' % (item, data, t.name))
    return data

def _lammps_parse_step(finp):
    # Read step
    data = _read_item(finp, 'ITEM: TIMESTEP')
    data = finp.readline()
    return int(float(data))

def _lammps_parse_system(finp):
    # Read Number of atoms
    data = _read_item(finp, 'ITEM: NUMBER OF ATOMS')
    data = finp.readline()
    npart = int(data)

    # Read box
    data = _read_item(finp, "ITEM: BOX BOUNDS")
    ndim = len(data.split()[3:]) # line is ITEM: BOX BONDS pp pp pp
    if len(data.split()) ==3 :
        ndim = 3
    L = []; offset = []
    for i in range(ndim):
        data = map(float, finp.readline().split())
        L.append(data[1] - data[0])
        offset.append(data[0])
    c = Cell(numpy.array(L))

    # Read positions and velocities
    imap = {1:'A', 2:'B', 3:'C', 4:'D'}
    data = _read_item(finp, "ITEM: ATOMS").split()
    # Determine how many variables are there
    n = len(data)-2
    ind = data.index('type')-2
    try:
        ix = data.index('x')-2
    except ValueError:
        ix = data.index('xu')-2
    
    # Rather well optimized now. Ignore velocities, do not fold particles.
    from itertools import islice
    data = ''.join(list(islice(finp, npart)))
    d = numpy.fromstring(data, sep=' ').reshape((npart, n))
    p = [Particle(int(d[i,ind]), imap[int(d[i,ind])], 1.0, d[i,ix:ix+3]) for i in xrange(npart)] #.fold(c))
    for pi in p:
        pi.fold(c)

    p.sort(key = lambda a : a.id)
    return System(particle=p, cell=c)

#@profilehooks.profile
def _lammps_parse_system_update(finp, system):
    # Read Number of atoms

    data = _read_item(finp, 'ITEM: NUMBER OF ATOMS')
    data = finp.readline()
    npart = int(data)

    if len(system.particle) == 0:
        system.particle = [Particle() for i in range(npart)]

    # Read box
    data = _read_item(finp, "ITEM: BOX BOUNDS")
    ndim = len(data.split()[3:]) # line is ITEM: BOX BONDS pp pp pp
    L = []; offset = []
    for i in range(ndim):
        data = map(float, finp.readline().split())
        L.append(data[1] - data[0])
        offset.append(data[0])
    c = Cell(numpy.array(L))

    system.cell = c

    # Read positions and velocities
    imap = {1:'A', 2:'B', 3:'C', 4:'D'}
    data = _read_item(finp, "ITEM: ATOMS")
    #p = []
    for i in range(npart):
        data = finp.readline().split()
        ids = int(data[0])
        r = numpy.array(map(float, data[1:4]))
        if (len(data)==7):
            v = numpy.array(map(float, data[4:7]))
        else:
            v = numpy.zeros(3)
        system.particle[i].id = ids
        system.particle[i].name = imap[ids]
        system.particle[i].position = r
        system.particle[i].velocity = v
        system.particle[i].fold(c)
        #map[ids] append(Particle(ids, imap[ids], 1.0, r, v).fold(c))
    system.particle.sort(key = lambda a : a.id)
    return  system

class TrajectoryLAMMPS(TrajectoryBase):

    """Trajectory layout for LAMMPS with additional .data file"""

    def __init__(self, fname, mode='r'):

        TrajectoryBase.__init__(self, fname, mode)

        # Assume cfgs are packed in some compressed tar file
        # They will be extracted inplace and deleted upon trajectory close()
        # TODO: this may be common to many Trajectories: isolate it
        base, ext = os.path.splitext(self.filename) # keep?
        dirname = os.path.dirname(self.filename) # keep?
        try:
            tar = tarfile.open(self.filename)
            tar.extractall(path=dirname)
            self._file_list = [os.path.join(dirname, f.name) for f in tar.getmembers()]
            self.tarred = True
        except:
            self.tarred = False
            self._file_list = [self.filename]

        file_steps = []
        for f in self._file_list:
            if os.path.isdir(f):
                continue
            fh = open(f, 'r')
            file_steps.append((f, _lammps_parse_step(fh)))
            fh.close()
        file_steps.sort(key = lambda a : a[1])
        self._file_samples = [a[0] for a in file_steps]
        self.steps = [a[1] for a in file_steps]
        #self._s = System()
        if self.tarred:
            tar.close()

    def read_sample(self, sample, unfolded=True):
        fh = open(self._file_samples[sample], 'r')
        s = _lammps_parse_step(fh)
        #self._s = _lammps_parse_system_update(fh, self._s)
        s = _lammps_parse_system(fh)
        fh.close()
        #return self._s
        return s

    # to get the tabulated potential from a dump of potential.x do 
    # > { echo -e "\nPOTENTIAL\nN 10000\n"; grep -v '#' /tmp/kalj.ff.potential.1-1 | awk '{printf "%i %g %12e %12e\n", NR, $1, $2, -$3}' ; }

    def write_init(self, system):
        # TODO: check this lammps write
        f = open(self.filename + '.data', 'w')
        np = len(system.particle)
        L = system.cell.side
        nsp = system.number_of_species

        # LAMMPS header
        h = '\n'
        h += "%i atoms\n" % np
        h += "%i atom types\n" % nsp
        h += "%g %g  xlo xhi\n" % (-L[0]/2, L[0]/2)
        h += "%g %g  ylo yhi\n" % (-L[1]/2, L[1]/2)
        h += "%g %g  zlo zhi\n" % (-L[2]/2, L[2]/2)
        f.write(h + '\n')

        # LAMMPS body

        # Masses of species
        m = "\nMasses\n\n"
        for isp in range(nsp):
            # Iterate over particles. Find instances of species and get masses
            for p in system.particle:
                if p.id == isp+1:
                    m += '%i %g\n' % (isp+1, p.mass)
                    break

        # Atom coordinates
        r = "\nAtoms\n\n"
        v = "\nVelocities\n\n"
        for i, p in enumerate(system.particle):
            r += '%i %i %g %g %g\n' % tuple([i+1, p.id] + list(p.position))
            v += '%i    %g %g %g\n' % tuple([i+1]       + list(p.velocity))

        f.write(m)
        f.write(r)
        f.write(v)
        f.close()
    
    def close(self):
        # Redefine it because we do not necessarily write cell here
        # Moreover we must delete files if it was open for writing
        # We remove possible directories at the end
        if not self.tarred:
            return
        d = []
        for f in self._file_list:
            if os.path.isdir(f):
                d.append(f)
            else:
                os.remove(f)
        for di in d:
            os.removedirs(di)
        return
