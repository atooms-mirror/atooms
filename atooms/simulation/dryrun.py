# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Dry run backend.

It just exposes a minimal backend interface.
"""


class DryRunBackend(object):

    def __init__(self, system=None):
        self.system = system
        if self.system is None:
            self.system = System()
        self.trajectory = Trajectory
        self.output_path = None
        self.steps = 0

    def write_checkpoint(self):
        pass

    @property
    def rmsd(self):
        return 0.0

    def run_pre(self, restart):
        pass

    def run_until(self, steps):
        self.steps = steps


class System(object):

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

    def __init__(self):
        self.temperature = 1.0


class Trajectory(object):

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
