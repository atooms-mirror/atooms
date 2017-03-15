# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import copy
import numpy
from .particle import position_cm, velocity_cm, fix_cm, total_kinetic_energy

class System(object):

    """System class."""

    def __init__(self, particle=None, cell=None, interaction=None, matrix=None, thermostat=None, dynamics=None):
        if particle is None:
            particle = []
        self.particle = particle
        self.interaction = interaction
        self.cell = cell
        self.matrix = matrix
        self.thermostat = thermostat
        self.dynamics = dynamics
        self._potential_energy = 0.0

    @property
    def number_of_dimensions(self):
        return len(self.particle[0].position)

    @property
    def number_of_species(self):
        return len(set(p.id for p in self.particle))

    def add_porous_matrix(self, matrix):
        self.matrix = copy.deepcopy(matrix)

    @property
    def density(self):
        return len(self.particle) / self.cell.volume

    @density.setter
    def density(self, rho):
        # Note This will fail if cell is None
        factor = (self.density / rho)**(1./3)
        for particle in self.particle:
            particle.position *= factor
        self.cell.side *= factor

    def temperature(self, ndof=None):
        # The n. of degrees of freedom can be passed to this method explicitly.
        # This way we can correct for missing translational invariance.
        # Ideally, one could determine this via some additional attribute.
        if ndof is None:
            ndof = (len(self.particle)-1) * self.number_of_dimensions
        return 2.0 / ndof * total_kinetic_energy(self.particle)

    def kinetic_energy(self):
        return total_kinetic_energy(self.particle)

    def kinetic_energy_per_particle(self):
        return total_kinetic_energy(self.particle) / len(self.particle)

    def potential_energy(self):
        return self._potential_energy

    def potential_energy_per_particle(self):
        return self._potential_energy / len(self.particle)

    def mean_square_displacement(self, reference):
        """
        Compute the mean square displacement of the system's particles
        with respect to those in the *reference* system.
        """
        displ = []
        for pi, pj in zip(self.particle, reference.particle):
            rij = numpy.array(pi.distance(pj, self.cell))
            displ.append(numpy.dot(rij, rij))
        return sum(displ) / len(self.particle)

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
            self.dynamics.evolve(self.particle, self.cell, self.interaction, self.thermostat)

    def maxwellian(self, temperature):
        """Reset velocities to a Maxwellian distribution with fixed CM."""
        T = temperature
        for p in self.particle:
            p.maxwellian(T)
        fix_cm(self.particle)
        # After fixing the CM the temperature is not exactly the targeted one
        # Therefore we scale the velocities so as to get to the right T
        T_old = self.temperature()
        fac = (T/T_old)**0.5
        for p in self.particle:
            p.velocity *= fac

    def dump(self, what, pslice=slice(None), order='C', dtype=None):
        """Dump system properties into numpy arrays.

        `what` can be either a string or a list. In the latter case a
        dict is returned. Each entry in what should be of the form
        particle.<attribute> or cell.<attribute>. The following
        aliases are allowed: pos, vel, ids, box.
        
        Particles' coordinates are thrown into (N, ndim)
        arrays if `order` is C or (ndim, N) arrays if `order` is F.

        It accepts particle slices `pslice` to filter out some
        particles. This can be also handled via trajectory decorators.
        """
        # Listify input variables
        if type(what) is str:
            what_list = [what]
        else:
            what_list = what
        if dtype is None:
            dtype_list = [None] * len(what_list)

        aliases = {'pos': 'particle.position', 
                   'vel': 'particle.velocity',
                   'ids': 'particle.id'}

        dump_db = {}
        for what, dtype in zip(what_list, dtype_list):
            # Accept some aliases
            if what in aliases:
                what_aliased = aliases[what]
            else:
                what_aliased = what
            # Extract the requested attribute
            attr = what_aliased.split('.')[-1]
            # Make array of attributes    
            if what_aliased.startswith('particle'):
                data = numpy.array([p.__getattribute__(attr) for p in self.particle], dtype=dtype)
            else:
                raise ValueError('Unknown attribute %s' % what_aliased)
            # We transpose the array if F order is requested (only meaningful for 2d arrays)
            if order == 'F':
                data = numpy.transpose(data)
            dump_db[what] = data

        # If what is a string or we only have one entry we return an
        # array, otherwise we return the whole dict
        if len(what_list) == 1:
            return dump_db.values()[0]
        else:
            return dump_db

    def scale(self, factor):
        """Rescale cell and particles' coordinates by *factor*."""
        for p in self.particle:
            p.position *= factor
        self.cell.side *= factor

    def report(self):
        txt =  "number of particles: %d\n" % len(self.particle)
        if self.cell is not None:
            txt += "cell side: %s\n" % self.cell.side
        return txt
