# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import copy
import numpy
from particle import position_cm, velocity_cm, fix_cm, total_kinetic_energy

class System(object):

    """
    System class
    """

    def __init__(self,name='Unknown',particle=[],cell=None,interaction=None,matrix=None,thermostat=None):
        self.name     = name
        self.particle = particle
        self.interaction = interaction
        self.cell     = cell
        self.matrix   = matrix
        self.thermostat = thermostat

    # TODO: drop copy(), use deepcopy directly
    def copy(self, other):
        """ Copy instance variables of an object to another """
        # TODO: copy attributes pythonically
        self.name     = copy.deepcopy(other.name)
        self.particle = copy.deepcopy(other.particle)
        self.interaction = copy.deepcopy(other.interaction)
        self.cell     = copy.deepcopy(other.cell)
        self.matrix   = copy.deepcopy(other.matrix)
        self.thermostat = copy.deepcopy(other.thermostat)

    @property
    def number_of_dimensions(self):
        return len(self.particle[0].position)
        #len(self.cell.side)

    @property
    def number_of_species(self):
        return len(set(p.id for p in self.particle))

    def add_porous_matrix(self,matrix):
        self.matrix = copy.deepcopy(matrix)

    def temperature(self):
        # TODO: add translation invariance check 
        ndf = (len(self.particle)-1) * len(self.particle[0].position)
        return 2.0 / ndf * total_kinetic_energy(self.particle) 

    def kinetic_energy(self):
        return total_kinetic_energy(self.particle) 

    @property
    def velocity_cm(self):
        # Wrapper to particle method
        return velocity_cm(self.particle)

    @property
    def position_cm(self):
        # Wrapper to particle method
        return position_cm(self.particle)

    def fix_cm(self):
        # Wrapper to particle method
        fix_cm(self.particle)
        
    def maxwellian(self, T):
        """ Reset velocities to a Maxwellian distribution with fixed CM """
        for p in self.particle:
            p.maxwellian(T)
        fix_cm(self.particle)
        # After fixing the CM the temperature is not exactly the targeted one
        # Therefore we scale the velocities so as to get to the right T
        T_old = self.temperature()
        fac = (T/T_old)**0.5
        for p in self.particle:
            p.velocity *= fac
        
    def dump(self, what, dim=slice(None), pslice=slice(None)):
        """ Thrown pos or vel into a big (N, ndim) numpy array. 
        It accepts particles slice although this should be handled via trajectory decorators"""
        if what == 'pos':
            return numpy.array([p.position[dim] for p in self.particle[pslice]])
        elif what == 'vel':
            return numpy.array([p.velocity[dim] for p in self.particle[pslice]])
