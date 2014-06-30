# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Trajectory formats for 3rd party simulation packages"""

import os
import re
import tarfile
import numpy

from atooms.trajectory import TrajectoryBase, SuperTrajectory
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

class TrajectoryLammps(TrajectoryBase):

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

        self._timestep = 1.0

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

import os
import glob
import tarfile
import random
from xml.etree import ElementTree

def map_label_id(names):
    lab = {}
    for i, t in enumerate(sorted(set(names))):
        lab[t] = i
    return lab

class TrajectoryHOOMDXML(TrajectoryBase):

    suffix = 'tgz'

    def __init__(self, fname, mode='r'):

        TrajectoryBase.__init__(self, fname, mode)
        base = os.path.basename(fname.split('.tgz')[0])
        
        if mode == 'r':
            self.__tmp_path = '/tmp/' + base
            try:
                os.makedirs(self.__tmp_path)
            except:
                pass

            tar = tarfile.open(fname)
            file_list = sorted([f.name for f in tar.getmembers()])
            tar.extractall(path=self.__tmp_path)
            tar.close()

            self.__f_samples = [os.path.join(self.__tmp_path, f) for f in file_list]

            cfg, box, pos, typ, vel = self.__read_one(self.__f_samples[0])

            self._timestep = 1.0
            self.steps = []
            for f in sorted(self.__f_samples):
                tree = ElementTree.parse(f)
                root = tree.getroot()
                cfg = root.find('configuration')
                self.steps.append(int(cfg.attrib['time_step']))

            # First sort them, then subtract out the first step
            self.steps.sort()
            self.steps = [s - self.steps[0] for s in self.steps]

        elif mode == 'w':
            pass

        elif mode == 'w:gz':
            self._tar = tarfile.open(fname, "w:gz")

    # TODO: rewind automatically when unfolding and restarting from a previous sample
    def rewind(self):
        self.trajectory.seek(0)

    def close(self):
        # clean up
        if self.mode == 'r':
            import shutil
            shutil.rmtree(self.__tmp_path)
        elif self.mode == 'w:gz':
            self._tar.close()

    def __read_one(self, fname):
        tree = ElementTree.parse(fname)
        root = tree.getroot()
        cfg = root.find('configuration')
        box = cfg.find('box')
        pos = cfg.find('position')
        try:
            vel = cfg.find('velocity')
        except:
            vel = None
        typ = cfg.find('type')
        pos_list = [map(float, r.split()) for r in pos.text.strip().split('\n')]
        box_list = [float(box.attrib[r]) for r in ['lx', 'ly', 'lz']]
        typ_list = typ.text.strip().split('\n')
        if not vel is None:
            vel_list = [map(float, v.split()) for v in vel.text.strip().split('\n')]
        else:
            vel_list = None
        return cfg, box_list, pos_list, typ_list, vel_list

    def read_sample(self, sample):
        cfg, box, pos, typ, vel = self.__read_one(self.__f_samples[sample])
        lab = map_label_id(typ)

        if vel is None:
            particle = [Particle(name=t, id=lab[t], position=numpy.array(p)) for p, t in zip(pos, typ)]
        else:
            particle = [Particle(name=t, id=lab[t], position=numpy.array(p), velocity=numpy.array(v)) 
                        for p, t, v in zip(pos, typ, vel)]
        cell = Cell(numpy.array(box))
        return System('unknown', particle, cell)

    def write_sample(self, system, step, ignore=[]):
        ndim = len(system.particle[0].position)
        n = len(system.particle)

        if self.mode == 'w':
            fname = self.filename.split('.tgz')[0] + '_step%010d.xml' % step
        else:
            base = os.path.basename(self.filename.split('.tgz')[0])
            fname = base + '_step%010d.xml' % step

        fh = open(fname, 'w')
        fh.write("""\
<?xml version="1.0" encoding="UTF-8"?>
<hoomd_xml version="1.4">
<configuration time_step="%d" dimensions="%d" natoms="%d" >
<box lx="%g" ly="%g" lz="%g"/>
""" % ((step, ndim, n) + tuple(system.cell.side)))

        fh.write("<position num=\"%d\">\n" % n)
        for p in system.particle:
            fh.write((ndim*" %g" +"\n") % tuple(p.position))
        fh.write("</position>\n")

        fh.write("<velocity num=\"%d\">\n" % n)
        for p in system.particle:
            fh.write((ndim*" %g" +"\n") % tuple(p.velocity))
        fh.write("</velocity>\n")

        fh.write("<type num=\"%d\">\n" % n)
        for p in system.particle:
            fh.write(("%s\n") % p.name)
        fh.write("</type>\n")

        fh.write("""\
</configuration>
</hoomd_xml>
""")
        fh.close()

        if self.mode == 'w:gz':
            self._tar.add(fname)
            os.remove(fname)


from xyz import TrajectoryXYZ

class TrajectoryXYZRUMD(TrajectoryXYZ):
    # TODO: allow reading unfolded configuration by parsing the box image integers

    # def __header_dict(self, line):
    #     params = {}
    #     for key, value in [d.split('=') for d in line.split()]:
    #         params[key] = value
    #     # Array entry have comma separated elements, split them into lists

    def __init__(self, filename, mode='r', basename='block'):
        """basename: prefix of RUMD configurations."""
        # Use an internal counter for ioformat=2
        self._step = 0
        self.basename = basename

        super(TrajectoryXYZRUMD, self,).__init__(filename, mode)
        
        # Redefine samples and steps to make sure these are the absolute steps and samples
        # This is important when trajectories are written in blocks.
        # To extract the block index we look at the filename indexing.
        # If the name is different the block index is set to zero and steps have no offset
        s = re.search(r'%s(\d*)' % basename, filename)
        if s:
            iblock = int(s.group(1))
        else:
            return

        # Redefine available steps to account for block offset
        dt = self.steps[-1]
        self.steps = [i+dt*iblock for i in self.steps]

    def _parse_step(self, data):
        s = re.search(r'timeStepIndex=(\d*)', data)
        if s is None:
            self._step += 1
            n = self._step
        else:
            n = s.group(1)
        return int(n)

    def _parse_cell(self):
        self.trajectory.seek(0)
        self.trajectory.readline()
        data = self.trajectory.readline()
        s = re.search(r'boxLengths=(\S*)', data)
        # Fix compatibility issue with more recent file format
        if s is None:
            s = re.search(r'sim_box=(\S*)', data)
        L = s.group(1)
        # TODO: improve parsing of timestep dt in xyz indexed files, we put it into _parse_cell() for the moment. We could have a parse metadata that returns a dict.
        s = re.search(r'dt=(\S*)', data)
        if s is None:
            self._timestep = 1.0
        else:
            self._timestep = float(s.group(1))
        return Cell(numpy.array(L.split(',')[1:], dtype=float))

    def _comment_header(self, step, system):

        def first_of_species(system, isp):
            for i, p in enumerate(system.particle):
                if p.id == isp:
                    return i
            raise ValueError('no species %d found in system' % isp)
            
        nsp = set([p.id for p in system.particle])
        mass = [system.particle[first_of_species(system, i)].mass for i in nsp]
        hdr = 'ioformat=1 dt=%g timeStepIndex=%d boxLengths=' + '%.8g,%.8g,%.8g' + ' numTypes=%d mass=' + '%.8g,'*(len(nsp)) + ' columns=type,x,y,z,vx,vy,vz\n'
        return hdr % tuple([self.timestep, step] + list(system.cell.side) + [len(nsp)] + mass)

    def write_sample(self, system, step, ignore=[]):
        # We need to redfine the id, because it expects numerical ids from 0 to nsp-1
        # We get the smallest species id, which we will then subtract.
        # TODO: cache id_min for efficiency
        id_min = min([p.id for p in system.particle])
        self.trajectory.write("%d\n" % len(system.particle))
        self.trajectory.write(self._comment_header(step, system))
        ndim = len(system.particle[0].position)
        for p in system.particle:
            self.trajectory.write(("%s"+ndim*" %14.6f" + ndim*" %g " + "\n") % ((p.id-id_min,) + tuple(p.position) + tuple(p.velocity)))

    def close(self):
        # We do not write the cell here in this format
        self.trajectory.close()

import glob

class SuperTrajectoryXYZRUMD(SuperTrajectory):
    
    def __new__(self, d, basename='block'):
        """ Takes a directory as input and get all block*gz files in there """
        if not os.path.isdir(d):
            raise IOError("We expected this to be a dir (%s)" % d)
        f_all = glob.glob(d + '/%s*gz' % basename)
        return SuperTrajectory([TrajectoryXYZRUMD(f, basename=basename) for f in f_all])
