# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Dry run backend.

It just exposes a minimal backend interface.
"""


class DryRun(object):

    "A simulation backend that performs no simulation at all."
    
    def __init__(self, system=None):
        self.system = system
        if self.system is None:
            self.system = System()
        self.trajectory = Trajectory
        self.output_path = None
        self.steps = 0

    def write_checkpoint(self):
        pass

    def read_checkpoint(self):
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

    def __init__(self):
        self.particle = []
        self.thermostat = Thermostat()

    def potential_energy(self, normed=False):
        return 0.

    def kinetic_energy(self, normed=False):
        return 0.

    def total_energy(self, normed=False):
        return 0.

    @property
    def density(self):
        return 0.

    @property
    def temperature(self):
        return 0.

    def mean_square_displacement(self, reference):
        return 0.


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
