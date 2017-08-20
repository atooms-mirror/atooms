# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Lammps simulation backend."""

import os
import copy
import subprocess
from atooms import trajectory
from atooms import simulation
from atooms import system
from atooms.system import thermostat
from atooms.trajectory import TrajectoryLAMMPS

try:
    out = subprocess.check_output('lammps < /dev/null', shell=True, stderr=subprocess.STDOUT)
    _version = out.split('\n')[0]
except subprocess.CalledProcessError:
    raise ImportError('lammps not installed')


class System(system.System):

    def __init__(self, filename, commands):
        """We accept xyz format on input."""
        self.filename = filename
        self.commands = commands
        if os.path.exists(filename):
            try:
                with trajectory.TrajectoryXYZ(filename) as t:
                    s = t[0]
            except:
                with trajectory.TrajectoryLAMMPS(filename) as t:
                    s = t[0]
    
            super(System, self).__init__(s.particle, s.cell, s.interaction, thermostat=s.thermostat)
        else:
            super(System, self).__init__()
               
    def potential_energy(self):
        """Full calculation of potential energy."""
        return 0.0

class LammpsBackend(object):

    def __init__(self, fileinp, commands):
        self.fileinp = fileinp
        self.commands = commands
        self.verbose = False
        if os.path.exists(commands):
            with open(commands) as fh:
                self.commands = fh.read()
        self.system = System(fileinp, self.commands)
        self.trajectory = TrajectoryLAMMPS
        self.steps = 0

    def __str__(self):
        return _version

    @property
    def rmsd(self):
        return 0.0

    def write_checkpoint(self):
        pass

    def run_pre(self, restart):
        return

    def run_until(self, steps):
        n = steps - self.steps
        # TODO: remove hard coded paths
        file_tmp = '/tmp/out.atom'
        # Update input file with current system
        file_inp = file_tmp + '.inp'
        # This will write the .inp startup file
        with TrajectoryLAMMPS(file_tmp, 'w') as th:
            th.write(self.system, 0)

        # Do things in lammps order: units, read, commands, run
        # A better approach would be to parse commands and place read_data after units
        # then pack commands again. Even better using PyLammps...
        cmd = """\
units		lj
atom_style	atomic
read_data %s
# read_restart
""" % file_inp
        cmd += self.commands
        cmd += """
run %s
write_dump all custom %s id type x y z vx vy vz modify sort id
""" % (n, file_tmp)

        # see https://stackoverflow.com/questions/163542/python-how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument
        p = subprocess.Popen(['lammps'], shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out = p.communicate(input=cmd)[0]
        if self.verbose:
            print out.decode()

        # Update internal reference system
        self.system = System(file_tmp, self.commands)
        self.steps = steps
