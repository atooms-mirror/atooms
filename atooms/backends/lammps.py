# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Minimal simulation backend for LAMMPS (http://lammps.sandia.gov).
"""

import os
import subprocess
import tempfile
from atooms import trajectory
from atooms import system
from atooms import interaction
from atooms.trajectory import TrajectoryLAMMPS
from atooms.core.utils import rmd

try:
    _ = subprocess.check_output('lammps < /dev/null', shell=True, stderr=subprocess.STDOUT)
    _version = _.decode().split('\n')[0][8:-1]
except subprocess.CalledProcessError:
    raise ImportError('lammps not installed')


def _run_lammps_command(cmd):
    # see https://stackoverflow.com/questions/163542/python-how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument
    p = subprocess.Popen(['lammps'], shell=True,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE)
    stdout = p.communicate(input=cmd.encode('utf8'))[0]
    code = p.returncode
    if code != 0:
        raise RuntimeError(stdout)
    return stdout.decode()


class Interaction(interaction.Interaction):

    """
    Interaction wrapper for LAMMPS.

    For the time being, it assumes `self.potential` is a string
    containing appropriate lammps commands that define the
    interaction.
    """

    # TODO: assign interaction to system based on pair_style entries in cmd

    def compute(self, observable, particle, cell):
        # We use self.potential as lammps commands
        dirout = tempfile.mkdtemp()
        file_tmp = os.path.join(dirout, 'lammps.atom')
        file_inp = os.path.join(dirout, 'lammps.atom.inp')
        # Update lammps startup file using self.system
        # This will write the .inp startup file
        with TrajectoryLAMMPS(file_tmp, 'w') as th:
            th.write(system.System(particle, cell), 0)

        # Do things in lammps order: units, read, commands, run. A
        # better approach would be to parse commands and place
        # read_data after units then pack commands again. Even better
        # using PyLammps...
        cmd = """\
units		lj
atom_style	atomic
read_data %s
%s
run 0
""" % (file_inp, self.potential)

        stdout = _run_lammps_command(cmd)
        found = False
        for line in stdout.split('\n'):
            if 'Step' in line:
                found = True
            elif found:
                self.energy = float(line.split()[2])
                break

        # Clean up
        rmd(dirout)


class System(system.System):

    """System wrapper for LAMMPS."""

    def __init__(self, filename, commands):
        """
        The input file `filename` must be in LAMMPS format or match a
        trajectory format recognized by atooms. LAMMPS `commands` must
        be a string and should not contain dump or run commands.
        """
        self._filename = filename
        self._commands = commands
        if os.path.exists(filename):
            # We accept any trajectory format, but if the format is
            # not recognized we force lammps native (atom) format
            try:
                with trajectory.Trajectory(filename) as t:
                    s = t[0]
            except ValueError:
                with trajectory.TrajectoryLAMMPS(filename) as t:
                    s = t[0]

            super(System, self).__init__(s.particle, s.cell,
                                         s.interaction, thermostat=s.thermostat)
        else:
            super(System, self).__init__()

        # Assign all commands as interaction potential, they should be stripped
        self.interaction = Interaction(commands)


class LAMMPS(object):

    """LAMMPS simulation backend."""

    version = _version

    def __init__(self, fileinp, commands):
        """
        The input file `filename` must be in LAMMPS format or match a
        trajectory format recognized by atooms. LAMMPS `commands` must
        be a string and should not contain dump or run commands.
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
        return 'LAMMPS'

    @property
    def rmsd(self):
        return 0.0

    def read_checkpoint(self):
        pass

    def write_checkpoint(self):
        pass

    def run(self, steps):
        dirout = tempfile.mkdtemp()
        file_tmp = os.path.join(dirout, 'lammps.atom')
        file_inp = os.path.join(dirout, 'lammps.atom.inp')
        # Update lammps startup file using self.system
        # This will write the .inp startup file
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

        stdout = _run_lammps_command(cmd)
        if self.verbose:
            print(stdout)

        # Update internal reference to self.system
        self.system = System(file_tmp, self.commands)

        # Clean up
        rmd(dirout)
