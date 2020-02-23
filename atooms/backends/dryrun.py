# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Dry run backend.

It just exposes a minimal backend interface.
"""

import copy


class DryRun(object):

    "A simulation backend that performs no simulation at all."

    version = '0.1.0'

    def __init__(self, system=None):
        self.system = system
        if self.system is None:
            self.system = System()
        self.trajectory_class = Trajectory
        self.steps = 0

    def __str__(self):
        return 'dryrun'

    def write_checkpoint(self, output_path):
        pass

    def read_checkpoint(self, output_path):
        pass

    @property
    def rmsd(self):
        return 0.0

    def run(self, steps):
        pass


class System(object):

    """
    The `System` class monkey patches the relevant methods of
    `atooms.system.System` required for a valid simulation backend.
    """

    def __init__(self, particle=None, cell=None, thermostat=None,
                 barostat=None, reservoir=None):
        self.particle = particle if particle is not None else []
        self.cell = None
        self.thermostat = Thermostat()
        self.barostat = None
        self.reservoir = None

    def potential_energy(self, per_particle=False, normed=False, cache=False):
        return 0.

    def kinetic_energy(self, per_particle=False, normed=False):
        return 0.

    def total_energy(self, per_particle=False, normed=False, cache=False):
        return 0.

    @property
    def density(self):
        return 0.

    @property
    def temperature(self):
        return 0.

    def set_temperature(self, T):
        pass

    def scale_velocities(self, factor):
        for p in self.particle:
            p.velocity *= -1.0

    def update(self, other, full=False, exclude=None, only=None):
        for key in other.__dict__:
            if exclude is not None or only is not None:
                if (exclude is not None and key not in exclude) or \
                   (only is not None and key in only):
                    self.__dict__[key] = copy.deepcopy(other.__dict__[key])
            else:
                if full or other.__dict__[key] is not None:
                    self.__dict__[key] = copy.deepcopy(other.__dict__[key])

    def report(self):
        return ''

class Thermostat(object):

    """A place-holder for the system thermostat and its temperature."""

    def __init__(self):
        self.temperature = 1.0


class Trajectory(object):

    """A minimal trajectory class interface."""

    suffix = 'dry'

    def __init__(self, filename, mode='r'):
        pass

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def write(self, system, step):
        pass

    def close(self):
        pass


class EnergyMinimization(object):

    "An optimization backend that performs no optimization at all."

    version = '0.1.0'

    def __init__(self, system):
        self.system = system
        self.trajectory = Trajectory
        self.output_path = None
        self.method = 'cg'
        self.tolerance = 1e-10
        self.max_iterations = 100000

    def run(self):
        pass
