# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich


class Thermostat(object):

    def __init__(self):
        self._temperature = 1.0

    def _get_temperature(self):
        return self._temperature

    def _set_temperature(self, value):
        self._temperature = value

    temperature = property(_get_temperature, _set_temperature, 'Temperature')


class System(object):

    def potential_energy(self):
        return 0.

    def kinetic_energy(self):
        return 0.

    def temperature(self):
        return 0.

    def mean_square_displacement(self, reference):
        return 0.

    @property
    def thermostat(self):
        return Thermostat()


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
                    

class DryRunBackend(object):

    def __init__(self, system=System()):
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


