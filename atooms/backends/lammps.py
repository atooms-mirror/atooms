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
from atooms.trajectory.decorators import change_species
from atooms.core.utils import rmd


# Local parallel environment
mpi_tasks = 1

# Check if lammps is installed
try:
    _ = subprocess.check_output('mpirun -n 1 lammps < /dev/null', shell=True,
                                stderr=subprocess.STDOUT, executable='/bin/bash')
    _version = _.decode().split('\n')[0][8:-1]
except subprocess.CalledProcessError:
    raise ImportError('lammps not installed')


def _run_lammps_command(cmd):
    # see https://stackoverflow.com/questions/163542/python-how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument
    p = subprocess.Popen(['mpirun -n {} lammps'.format(mpi_tasks)],
                         shell=True,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         executable='/bin/bash')
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
                self.energy = float(line.split()[2]) * len(particle)
                break

        # Clean up
        rmd(dirout)


# We use the base system class
System = system.System


class LAMMPS(object):

    """LAMMPS simulation backend."""

    version = _version

    def __init__(self, inp, commands):
        """
        We initialize the backend from `inp`, which can be a `System`, a
        `Trajectory` or path to a trajectory. LAMMPS `commands` must
        be a string or a file and should not contain dump or run
        commands.
        """
        self.verbose = False

        # Initialize commands
        self.commands = commands
        if os.path.exists(commands):
            with open(commands) as fh:
                self.commands = fh.read()

        # Define the initial system
        if isinstance(inp, system.System):
            # If we pass a system there is nothing to do
            self.system = inp

        elif isinstance(inp, trajectory.base.TrajectoryBase):
            # It is trajectory, we get the last frame
            self.system = inp[-1]

        elif os.path.exists(inp):
            # We accept any trajectory format, but if the format is
            # not recognized we force lammps native (atom) format
            try:
                with trajectory.Trajectory(inp) as t:
                    # We enforce fortran species layout
                    t.add_callback(change_species, 'F')
                    s = t[-1]
            except ValueError:
                with trajectory.TrajectoryLAMMPS(inp) as t:
                    s = t[-1]
            self.system = s

        else:
            raise ValueError('could not initialize system from {}'.format(inp))

        # Default trajectory format
        self.trajectory = TrajectoryLAMMPS

        # Assign commands as potentials, they should be stripped
        self.system.interaction = Interaction(commands)

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

        # Set fixes from the system if we find thermostat / barostat
        if self.system.thermostat is not None and self.system.barostat is not None:
            # NPT ensemble
            fix = 'fix 1 all npt temp {0.temperature} {0.temperature} {0.relaxation_time} iso {1.pressure} {1.pressure} {1.relaxation_time}'.format(self.system.thermostat, self.system.barostat)
        if self.system.thermostat is not None:
            # NVT ensemble
            fix = 'fix 1 all nvt temp {0.temperature} {0.temperature} {0.relaxation_time}'.format(self.system.thermostat)
        elif not 'fix' in self.commands:
            # NVE ensemble
            fix = 'fix 1 all nve'
        else:
            # The integrator is already contained in the commands
            fix = ''

        # Do things in lammps order: units, read, commands, run. A
        # better approach would be to parse commands and place
        # read_data after units then pack commands again. Even better
        # using PyLammps...
        cmd = """\
units		lj
atom_style	atomic
read_data {}
{}
{}
run {}
write_dump all custom {} id type x y z vx vy vz modify sort id
""".format(file_inp, self.commands, fix, steps, file_tmp)

        stdout = _run_lammps_command(cmd)
        if self.verbose:
            print(stdout)

        # Update internal reference to self.system
        # Note that the thermostat and barostat are not touched
        new_system = TrajectoryLAMMPS(file_tmp)[-1]
        for i in range(len(self.system.particle)):
            self.system.particle[i] = new_system.particle[i]

        # Clean up
        rmd(dirout)
