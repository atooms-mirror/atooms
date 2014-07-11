# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Trajectory formats for 3rd party simulation packages"""

import os
import re
import numpy

from atooms.trajectory.xyz import TrajectoryXYZ
from atooms.trajectory import SuperTrajectory
from atooms.system.cell import Cell


class TrajectoryRUMD(TrajectoryXYZ):
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

        super(TrajectoryRUMD, self,).__init__(filename, mode)
        
        # The minimum id for RUMD is 0
        self._min_id = 0
        
        if basename == 'block':
            # Redefine samples and steps to make sure these are the absolute steps and samples
            # This is important when trajectories are written in blocks.
            # To extract the block index we look at the filename indexing.
            # If the name is different the block index is set to zero and steps have no offset
            s = re.search(r'%s_(\d*)' % basename, filename)
            if s:
                iblock = int(s.group(1))
                # Redefine available steps to account for block offset
                dt = self.steps[-1]
                self.steps = [i+dt*iblock for i in self.steps]
        else:
            # In case of non native RUMD filename, we assume the step
            # is written after the basename.
            s = re.search(r'%s_(\d*)' % basename, filename)
            if s:
                self.steps = [int(s.group(1))]

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
        side = s.group(1)
        # TODO: improve parsing of timestep dt in xyz indexed files, we put it into _parse_cell() for the moment. We could have a parse metadata that returns a dict.
        s = re.search(r'dt=(\S*)', data)
        if s is None:
            self._timestep = 1.0
        else:
            self._timestep = float(s.group(1))
        return Cell(numpy.array(side.split(',')[1:], dtype=float))

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

class SuperTrajectoryRUMD(SuperTrajectory):
    
    def __new__(self, inp, basename='block'):
        """ Takes a directory as input and get all block*gz files in there """
        if not os.path.isdir(inp):
            raise IOError("We expected this to be a dir (%s)" % inp)
        f_all = glob.glob(inp + '/%s*gz' % basename)
        return SuperTrajectory([TrajectoryRUMD(f, basename=basename) for f in f_all])
