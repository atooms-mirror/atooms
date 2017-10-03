# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Minimal simulation backend for LAMMPS (http://lammps.sandia.gov).
"""

import os
import subprocess
from atooms import trajectory
from atooms import system
from atooms.trajectory import TrajectoryLAMMPS

try:
    _ = subprocess.check_output('lammps < /dev/null', shell=True, stderr=subprocess.STDOUT)
    _version = _.split('\n')[0]
except subprocess.CalledProcessError:
    raise ImportError('lammps not installed')


class System(system.System):

    def __init__(self, filename, commands):
        """
        The input trajectory file `filename` can be any trajectory format
        recognized by atooms. Lammps `commands` are passed as a string
        and should not contain dump and run commands.
        """
        self._filename = filename
        self._commands = commands
        if os.path.exists(filename):
            # We accept any trajectory format, but if the format is
            # not recognized we force lammps native (atom) format
            try:
                with trajectory.Trajectory(filename) as t:
                    s = t[0]
            except:
                with trajectory.TrajectoryLAMMPS(filename) as t:
                    s = t[0]

            super(System, self).__init__(s.particle, s.cell,
                                         s.interaction, thermostat=s.thermostat)
        else:
            super(System, self).__init__()

    def potential_energy(self, normed=False):
        """Full calculation of potential energy."""
        return 0.0


class LAMMPS(object):

    """LAMMPS simulation backend."""

    def __init__(self, fileinp, commands):
        """
        The input trajectory file `fileinp` can be any trajectory format
        recognized by atooms. Lammps `commands` are passed as a string
        and should not contain dump and run commands.
        """
        self.fileinp = fileinp
        self.commands = commands
        self.verbose = False
        if os.path.exists(commands):
            with open(commands) as fh:
                self.commands = fh.read()
        self.system = System(fileinp, self.commands)
        self.trajectory = TrajectoryLAMMPS

    def __str__(self):
        return _version

    @property
    def rmsd(self):
        return 0.0

    def read_checkpoint(self):
        pass

    def write_checkpoint(self):
        pass

    def run(self, steps):
        # TODO: remove hard coded paths
        file_tmp = '/tmp/out.atom'
        # Update lammps startup file using self.system
        # This will write the .inp startup file
        file_inp = file_tmp + '.inp'
        with TrajectoryLAMMPS(file_tmp, 'w') as th:
            th.write(self.system, 0)

        # Do things in lammps order: units, read, commands, run. A
        # better approach would be to parse commands and place
        # read_data after units then pack commands again. Even better
        # using PyLammps...
        cmd = """\
units		lj
atom_style	atomic
read_data %s
%s
run %s
write_dump all custom %s id type x y z vx vy vz modify sort id
""" % (file_inp, self.commands, steps, file_tmp)

        # see https://stackoverflow.com/questions/163542/python-how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument
        p = subprocess.Popen(['lammps'], shell=True,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE)
        stdout = p.communicate(input=cmd)[0]
        code = p.returncode
        if code != 0:
            raise RuntimeError(stdout)
        if self.verbose:
            print stdout.decode()

        # Update internal reference to self.system
        self.system = System(file_tmp, self.commands)
