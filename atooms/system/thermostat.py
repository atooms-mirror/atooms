# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Thermostat to control the temperature during a simulation."""

class Thermostat(object):

    def __init__(self, name, temperature, mass=1.0, collision_period=-1):
        self.name = name
        self.temperature = temperature
        self.mass = mass
        self.collision_period = collision_period
