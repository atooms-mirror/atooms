# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

class DryRunBackend(object):

    def __init__(self, system):
        self.system = system
        self.trajectory = Trajectory
        self.output_path = None

    def write_checkpoint(self):
        pass

    @property
    def rmsd(self):
        return 0.0

    def run_pre(self, restart):
        pass

    def run_until(self, n):
        self.steps = n
        pass


class System(object):
    
    def potential_energy(self):
        return 0.

    def temperature(self):
        return 0.

    def mean_square_displacement(self, reference):
        return 0.

    @property
    def thermostat(self):
        pass

class Trajectory(object):

    suffix = 'dry'

    def __init__(self, filename, mode='r'):
        pass

    def write(self, system, step):
        pass

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
    
