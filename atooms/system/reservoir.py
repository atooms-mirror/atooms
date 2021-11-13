"""
Thermodynamic reservoirs to impose temperature, pressure, or
particle numbers in a simulation.

When added to a `System` object, they will determine the statistical
ensemble of a simulation:

    #!python
    system = System()
    system.thermostat = Thermostat(temperature=1.0)
    system.barostat = Barostat(pressure=1.0)
    print system.ensemble == 'NPT'

Note: in an actual simulation backend, these reservoirs will have
additional degrees of freedom, e.g. s and pi in a Nose-like
thermostat.
"""

import numpy


class Thermostat(object):

    """Thermostat to control the temperature during a simulation."""

    def __init__(self, temperature, name='', mass=1.0, coordinate=1.0,
                 momentum=1.0, relaxation_time=1.0,
                 collision_period=None):
        self.name = name
        self.temperature = temperature
        self.coordinate = numpy.asarray(coordinate)
        self.momentum = numpy.asarray(momentum)
        self.mass = numpy.asarray(mass)
        self.relaxation_time = relaxation_time
        if collision_period is not None:
            self.relaxation_time = collision_period

    @property
    def collision_period(self):
        return self.relaxation_time


class Barostat(object):

    """Barostat to control the pressure during a simulation."""

    def __init__(self, pressure, name='', mass=1.0, coordinate=1.0,
                 momentum=1.0, relaxation_time=1.0):
        self.name = name
        self.pressure = pressure
        self.coordinate = numpy.asarray(coordinate)
        self.momentum = numpy.asarray(momentum)
        self.mass = numpy.asarray(mass)
        self.relaxation_time = relaxation_time


class Reservoir(object):

    """Reservoir to control the number of particles during a simulation."""

    def __init__(self, chemical_potential, name='', mass=1.0):
        self.name = name
        self.chemical_potential = chemical_potential
        self.mass = mass
