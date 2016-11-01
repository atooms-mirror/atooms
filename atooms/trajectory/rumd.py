# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""RUMD trajectory format."""

import os
import re
import glob

from atooms.trajectory.xyz import TrajectoryXYZ
from atooms.trajectory import SuperTrajectory


class TrajectoryRUMD(TrajectoryXYZ):
    # TODO: allow reading unfolded configuration by parsing the box image integers

    def __init__(self, filename, mode='r'):
        # Use an internal counter for ioformat=2
        self._timestep = 1.0
        alias = {'timeStepIndex': 'step', 'boxLengths':'cell', 'sim_box':'cell'}
        super(TrajectoryRUMD, self).__init__(filename, mode, alias=alias)
        # The minimum id for RUMD is 0
        self._min_id = 0

    def _setup_steps(self):
        super(TrajectoryRUMD, self)._setup_steps()

        # RUMD specific stuff
        basename = os.path.basename(self.filename)
        s = re.search(r'([a-zA-Z0-9]+)_(\d+)', basename)
        if s:
            base, step = s.group(1), s.group(2)
            # For native rumd trajectories we add the block offset
            if base == 'block':
                # Redefine samples and steps to make sure these are the absolute steps and samples
                # This is important when trajectories are written in blocks.
                # To extract the block index we look at the filename indexing.
                # If the name is different the block index is set to zero and steps have no offset
                iblock = int(step)
                dt = self.steps[-1]
                self.steps = [i+dt*iblock for i in self.steps]
            else:
                # In case of non native RUMD filename, we assume the step
                # is written after the basename.
                if len(self.steps) == 1:
                    self.steps = [int(step)]

    def _read_metadata(self, sample):
        meta = super(TrajectoryRUMD, self)._read_metadata(sample)
        # RUMD specific stuff that can't be handled as aliases.
        if 'integrator' in meta:
            self._timestep = meta['integrator'][1]
            meta['dt'] = meta['integrator'][1]
        if 'sim_box' in meta:
            # After sim_box there is a keyword for the box type which we ignore
            meta['cell'] = meta['sim_box'][1:]
        return meta

    def read_timestep(self):
        meta = self._read_metadata(0)
        if 'dt' in meta:
            return meta['dt']
        else:
            return 1.0

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

    def write_sample(self, system, step):
        # We need to redfine the id, because it expects numerical ids from 0 to nsp-1
        # We get the smallest species id, which we will then subtract.
        id_min = min([p.id for p in system.particle])
        self.trajectory.write("%d\n" % len(system.particle))
        self.trajectory.write(self._comment_header(step, system))
        ndim = len(system.particle[0].position)
        for p in system.particle:
            self.trajectory.write(("%s"+ndim*" %14.6f" + ndim*" %g " + "\n") % ((p.id-id_min,) + tuple(p.position) + tuple(p.velocity)))

    def close(self):
        # We do not write the cell here in this format
        self.trajectory.close()


class SuperTrajectoryRUMD(SuperTrajectory):

    def __new__(self, inp, variable=False, periodic=True, basename='block'):
        """ Takes a directory as input and get all block*gz files in there """
        if not os.path.isdir(inp):
            raise IOError("We expected this to be a dir (%s)" % inp)
        f_all = glob.glob(inp + '/%s*gz' % basename)
        return SuperTrajectory(f_all, TrajectoryRUMD, variable=variable, periodic=periodic)
