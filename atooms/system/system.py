# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import copy
import numpy
from particle import position_cm, velocity_cm, fix_cm, total_kinetic_energy

ndim = 3

class System(object):

    """System class."""

    def __init__(self, particle=[], cell=None, interaction=None, matrix=None, thermostat=None, dynamics=None):
        self.particle = particle
        self.interaction = interaction
        self.cell     = cell
        self.matrix   = matrix
        self.thermostat = thermostat
        self.dynamics = dynamics
        self._potential_energy = 0.0

    @property
    def number_of_dimensions(self):
        return len(self.particle[0].position)

    @property
    def number_of_species(self):
        return len(set(p.id for p in self.particle))

    def add_porous_matrix(self,matrix):
        self.matrix = copy.deepcopy(matrix)

    @property
    def density(self):
        return len(self.particle) / self.cell.volume

    def temperature(self, ndof=None):
        # The n. of degrees of freedom can be passed to this method
        # This way we can correct for missing translational invariance.
        # Ideally, one could determine this via some additional attribute.
        if ndof is None:
            ndof = (len(self.particle)-1) * self.number_of_dimensions
        return 2.0 / ndof * total_kinetic_energy(self.particle) 

    def kinetic_energy(self):
        return total_kinetic_energy(self.particle) 

    def potential_energy(self):
        return self._potential_energy

    def mean_square_displacement(self, reference):
        return 0.0

    @property
    def velocity_cm(self):
        return velocity_cm(self.particle)

    @property
    def position_cm(self):
        return position_cm(self.particle)

    def fix_cm(self):
        fix_cm(self.particle)

    def evolve(self):
        # Time evolution is a behavior of a system.
        # But to avoid passing the whole object, we must unpack it
        if not self.dynamics is None:
            self.dynamics.evolve(self.particle,self.cell,self.interaction,self.thermostat)
        
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
        elif what == 'sigma':
            return numpy.array([p.radius*2 for p in self.particle[pslice]])

    def scale(self, factor):
        """Affine displacement"""
        for p in self.particle:
            p.position *= factor
        self.cell.side *= factor
