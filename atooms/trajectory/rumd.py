import os
import re
import glob

from atooms.system.particle import distinct_species
from atooms.trajectory.xyz import TrajectoryXYZ
from atooms.trajectory import SuperTrajectory


class TrajectoryRUMD(TrajectoryXYZ):

    """RUMD trajectory format (https://rumd.org)"""

    # TODO: allow reading unfolded configuration by parsing the box image integers

    def __init__(self, filename, mode='r', fields=None):
        # We do not currently support writing with custom fields
        # The format is: species, position, velocity
        if mode == 'w' and fields is not None:
            raise ValueError('TrajectoryRUMD cannot write custom fields')
        super(TrajectoryRUMD, self).__init__(filename, mode,
                                             alias={'timeStepIndex': 'step',
                                                    'boxLengths': 'cell',
                                                    'sim_box': 'cell'},
                                             fields=fields)
        # The minimum id for RUMD is 0
        self._min_id = 0

    def read_steps(self):
        steps = super(TrajectoryRUMD, self).read_steps()

        # RUMD specific stuff
        basename_ext = os.path.basename(self.filename)
        basename = basename_ext.split('.')[0]

        # For native rumd trajectories we add the block offset
        s = re.search(r'(trajectory|block)(\d+)', basename)
        if s:
            # Redefine samples and steps to make sure these are the absolute step indexes
            # This is important when trajectories are written in blocks.
            _, block = s.group(1), s.group(2)
            iblock = int(block)
            dt = steps[-1]
            steps = [i + dt*iblock for i in steps]

        # If we follow a folder based logic, files are named according
        # to the step and contain a single step like
        #     folder/0000001.xyz.gz
        # We grab the step from
        # the file name.
        s = re.search(r'^(\d+)$', basename)
        if s and len(steps) == 1:
            steps = [int(basename)]
        return steps

    def _read_comment(self, frame):
        meta = super(TrajectoryRUMD, self)._read_comment(frame)

        # RUMD specific stuff that can't be handled as aliases.
        if 'integrator' in meta:
            meta['dt'] = meta['integrator'][1]
        if 'sim_box' in meta:
            # After sim_box there is a keyword for the box type which we ignore
            meta['cell'] = meta['sim_box'][1:]
            meta['ndim'] = len(meta['cell'])

        return meta

    def _comment(self, step, system):

        def first_of_species(system, species):
            for i, p in enumerate(system.particle):
                if p.species == species:
                    return i
            raise ValueError('no species %d found in system' % species)

        sp = distinct_species(system.particle)
        mass = [system.particle[first_of_species(system, isp)].mass for isp in sp]
        hdr = 'ioformat=1 dt=%g timeStepIndex=%d boxLengths=' + \
              '%.12f,%.12f,%.12f' + \
              ' numTypes=%d mass=' + '%g,'*(len(sp))
        hdr = hdr.strip(',')
        hdr += ' columns=type,x,y,z,vx,vy,vz\n'
        return hdr % tuple([self.timestep, step] + list(system.cell.side) + [len(sp)] + mass)

    def write_system(self, system, step):
        self._setup_format()
        sp = distinct_species(system.particle)
        self._file.write("%d\n" % len(system.particle))
        self._file.write(self._comment(step, system))
        ndim = len(system.particle[0].position)
        for p in system.particle:
            # We get the integer index corresponding to species Ex.:
            # if species are 'A', 'B' we get 0 and 1. Note that in
            # general getting the sample back via read_system() will
            # not preserve the species.
            isp = sp.index(p.species)
            self._file.write("{0} {1.position} {1.velocity}\n".format(isp, p))


class SuperTrajectoryRUMD(SuperTrajectory):

    """SuperTrajectory for RUMD format"""

    # TODO: why new here? init should be enough and will make it possible dynamic extension
    def __new__(cls, inp, mode='r', basename='trajectory*.gz'):
        """ Takes a directory as input and get all block*gz files in there """
        if not os.path.isdir(inp):
            raise IOError("We expected this to be a dir (%s)" % inp)
        f_all = glob.glob(inp + '/' + basename)
        if len(f_all) == 0:
            # Let's try with 00000.xyz.gz lie files
            f_all = glob.glob(inp + '/[0-9]*gz')
            f_all.sort()
        elif len(f_all) > 0:
            # Avoid last block because rumd does not write the last cfg!
            f_all.sort()
            f_all = f_all[:-1]
        return SuperTrajectory(f_all, TrajectoryRUMD)
