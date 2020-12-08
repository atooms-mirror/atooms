# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
The physical system at hand.

The system of interest in a classical atomistic simulations is
composed of interacting point particles, usually enclosed in a
simulation cell. The system may be in contact with a thermostat, a
barostat or a particle reservoir.
"""

import copy
import numpy
from .particle import cm_position, cm_velocity, fix_total_momentum


class System(object):

    """System class."""

    def __init__(self, particle=None, cell=None, interaction=None,
                 thermostat=None, barostat=None, reservoir=None):
        if particle is None:
            particle = []
        self.particle = particle
        """A list of `Particle` instances."""
        self.interaction = interaction
        self.cell = cell
        self.thermostat = thermostat
        self.barostat = barostat
        self.reservoir = reservoir
        self.matrix = None
        # Internal data dictionary for array dumps
        self._data = None

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        # Note: when making a shallow copy the _data cache
        # may need to be cleared. We might set self._data to None here.
        return result

    def update(self, other, full=False, exclude=None, only=None):
        """
        Update current system attributes in-place using the `other` System
        as source.
        
        The default behavior is to make deep copies only of all the
        `other` system attributes that are not None, i.e. which are
        set.

        To overwrite all attributes, even those set in `self` but not
        in `other`, use `full=True`.

        The lists `exclude` and `only` can be used to avoid updating
        some system attributes or to update only some system
        attributes, respectively.

        This method can be used by subclasses to deal with performace
        or memory dellocation issues.
        """
        for key in other.__dict__:
            if exclude is not None or only is not None:
                if (exclude is not None and key not in exclude) or \
                   (only is not None and key in only):
                    self.__dict__[key] = copy.deepcopy(other.__dict__[key])
            else:
                if full or other.__dict__[key] is not None:
                    self.__dict__[key] = copy.deepcopy(other.__dict__[key])

        # Internal data dictionary for array dumps
        self._data = None

    @property
    def number_of_dimensions(self):
        """
        Number of spatial dimensions, guessed from the length of
        `particle[0].position`.
        """
        if len(self.particle) > 0:
            return len(self.particle[0].position)
        else:
            return 0

    def distinct_species(self):
        """Sorted list of distinct chemical species in the system."""
        return sorted(set(p.species for p in self.particle))

    @property
    def density(self):
        """
        Density of the system.

        It `cell` is None, an estimate is given.
        """
        if len(self.particle) == 0:
            return 0.0
        
        if self.cell is None:
            # Estimate density assuming V is the product of the
            # largest distances along each dimension
            volume = 1.0
            for axis in range(len(self.particle[0].position)):
                x_0 = min([p.position[axis] for p in self.particle])
                x_1 = max([p.position[axis] for p in self.particle])
                volume *= (x_1 - x_0)
            return len(self.particle) / volume
        else:
            return len(self.particle) / self.cell.volume

    @density.setter
    def density(self, rho):
        self.set_density(rho)

    def set_density(self, rho):
        """Set the system density to `rho` by rescaling the coordinates."""
        if self.cell is None:
            return ValueError('cannot compute density without a cell')
        factor = (self.density / rho)**(1./len(self.cell.side))
        for particle in self.particle:
            particle.position *= factor
        self.cell.side *= factor

    @property
    def packing_fraction(self):
        from math import pi
        return pi / 6 * sum([(2 * p.radius)**3 for p in self.particle]) / self.cell.volume

    @property
    def temperature(self):
        """
        Kinetic temperature.
        """
        # TODO: determine translational invariance via some additional attribute.
        if len(self.particle) > 1:
            # Number of degrees of freddom is (N-1)*ndim
            ndof = (len(self.particle)-1) * self.number_of_dimensions
            return 2.0 / ndof * self.kinetic_energy()
        elif len(self.particle) == 1:
            # Pathological case
            return 2.0 / self.number_of_dimensions * self.kinetic_energy()
        else:
            # Empty particle list
            return 0.0

    @temperature.setter
    def temperature(self, T):
        self.set_temperature(T)

    def set_temperature(self, temperature):
        """Reset velocities to a Maxwellian distribution with fixed CM."""
        T = temperature
        for p in self.particle:
            p.maxwellian(T)
        fix_total_momentum(self.particle)
        # After fixing the CM the temperature is not exactly the targeted one
        # Therefore we scale the velocities so as to get to the right T
        T_old = self.temperature
        if T_old > 0.0:
            fac = (T/T_old)**0.5
            for p in self.particle:
                p.velocity *= fac

    def scale_velocities(self, factor):
        """Scale particles' velocities by `factor`."""
        for p in self.particle:
            p.velocity *= factor

    def kinetic_energy(self, per_particle=False, normed=False):
        """
        Return the total kinetic energy of the system.

        If `per_particle` or `normed` is `True`, return the kinetic
        energy per particle.
        """
        ekin = sum([p.kinetic_energy for p in self.particle])
        if per_particle or normed:
            return ekin / len(self.particle)
        else:
            return ekin

    def compute_interaction(self, what):
        """
        Compute interaction in the system
        """
        if self.interaction is not None:
            self.interaction.compute(what, self.particle, self.cell)

    def potential_energy(self, per_particle=False, normed=False, cache=False):
        """
        Return the total potential energy of the system.

        If `per_particle` or `normed` is `True`, return the potential
        energy per particle.
        """
        if self.interaction is not None:
            if not cache:
                self.compute_interaction('forces')
            if per_particle or normed:
                return self.interaction.energy / len(self.particle)
            else:
                return self.interaction.energy
        else:
            return 0.0

    def total_energy(self, per_particle=False, normed=None, cache=False):
        """
        Return the total energy of the system.

        If `per_particle` or `normed` is `True`, return the total
        energy per particle.
        """
        if normed is not None:
            per_particle = normed
        return self.potential_energy(per_particle, cache) + self.kinetic_energy(per_particle)

    def force_norm(self, per_particle=True, cache=False):
        """
        Return the norm of the force vector.

        If `per_particle` is `True`, return the force norm per particle.
        """
        if self.interaction is not None:
            if not cache:
                self.compute_interaction('forces')
            if per_particle:
                return numpy.sum(self.interaction.forces**2)**0.5 / len(self.particle)
            else:
                return numpy.sum(self.interaction.forces**2)**0.5
        else:
            return 0.0

    def force_norm_square(self, per_particle=True, cache=False):
        """
        Return the squared norm of the force vector.

        If `per_particle` is `True`, return the force squared norm per particle.
        """
        if self.interaction is not None:
            if not cache:
                self.compute_interaction('forces')
            if per_particle:
                return numpy.sum(self.interaction.forces**2) / len(self.particle)
            else:
                return numpy.sum(self.interaction.forces**2)
        else:
            return 0.0

    def virial(self, per_particle=True, cache=False):
        """
        Return the virial of the system.

        If `per_unit_volume` is `True`, return the virial per particle.
        """
        if self.interaction is not None:
            if not cache:
                self.compute_interaction('forces')
            if per_particle:
                return self.interaction.virial / len(self.particle)
            else:
                return self.interaction.virial
        else:
            return 0.0

    @property
    def pressure(self):
        """
        Return the pressure of the system.

        It assumes that `self.interaction` has already been computed.
        """
        if self.thermostat:
            T = self.thermostat.temperature
        else:
            T = self.temperature
        return (len(self.particle) * T + self.interaction.virial / self.number_of_dimensions) / self.cell.volume

    @property
    def cm_velocity(self):
        """Center-of-mass velocity."""
        return cm_velocity(self.particle)

    @property
    def cm_position(self):
        """Center-of-mass position."""
        return cm_position(self.particle)

    def fix_momentum(self):
        """Subtract out the the center-of-mass motion."""
        fix_total_momentum(self.particle)

    def fold(self):
        """Fold positions into central cell."""
        for p in self.particle:
            p.fold(self.cell)

    def dump(self, what=None, order='C', dtype=None, view=False, clear=False, flat=False):
        """
        Return a numpy array with system properties specified by `what`.

        If `what` is a string, it must be of the form
        `particle.<attribute>` or `cell.<attribute>`. The following
        aliases are allowed:

        - `pos` (`particle.position`)
        - `vel` (`particle.velocity`)
        - `spe` (`particle.species`)

        If `what` is a list of strings of the form above, a dict of
        numpy arrays is returned with `what` as keys.

        If `view` is `True`, the requested particle property is set as
        a view on the corresponding portion of the dump array. This
        allows one to modify the dump array efficiently and keep the
        particle properties in sync. Numpy arrays must be modified
        in-place to keep everything in sync.

        If `clear` is True, a new view is created on the requested
        particle property.

        If `flat` is True, the resulting array is flatted following
        `order` using numpy.flatten().

        Particles' coordinates are returned as (N, ndim) arrays if
        `order` is `C` or (ndim, N) arrays if `order` is `F`.

        Examples:
        --------
        These numpy arrays are element-wise identical

            #!python
            pos = system.dump('particle.position')
            pos = system.dump('position')
            pos = system.dump('pos')

        Return a dict with both positions and velocities

            #!python
            dump = system.dump(['pos', 'vel'])
        """
        # Unless a view is requested the default behavior is to always
        # create a new dump. This allows changes in the particles'
        # properties or particle number to be reflected in the dump.
        if not view:
            clear = True

        # Setup data dictionary
        if self._data is None or clear:
            self._data = {}

        # Return immediately if we only want to clear the dump
        if what is None and clear:
            return

        # Listify input variables
        if type(what) is str:
            what_list = [what]
            dtype_list = [dtype]
        else:
            what_list = what
            dtype_list = dtype
            if dtype is None:
                dtype_list = [None] * len(what_list)

        # Accepts some aliases
        aliases = {
            'box': 'cell.side',
            'pos': 'particle.position',
            'vel': 'particle.velocity',
            'spe': 'particle.species',
            'rad': 'particle.radius',
            'radius': 'particle.radius',
            'position': 'particle.position',
            'velocity': 'particle.velocity',
            'species': 'particle.species'}
        for i, what in enumerate(what_list):
            if what in aliases:
                what_list[i] = aliases[what]

        for what, dtype in zip(what_list, dtype_list):

            # Extract the requested attribute
            attr = what.split('.')[-1]

            # Make array of attributes
            if what.startswith('particle'):

                # Skip if it has been dumped already
                # and the number of particles has not changed
                # Note: if particles are reassigned the dump will
                # not be updated. It can be fixed by keeping a list of ids
                # although testing it would bring a overhead.
                if what in self._data:
                    if not flat:
                        if order == 'F':
                            npart = self._data[what].shape[-1]
                        else:
                            npart = self._data[what].shape[0]
                    else:
                        ndims = self.number_of_dimensions
                        npart = int(len(self._data[what]) / ndims)
                        
                    if npart == len(self.particle):
                        continue
                    
                data = numpy.array([getattr(p, attr) for p in self.particle], dtype=dtype)
                
                # We transpose the array if F order is requested
                if order == 'F':
                    data = numpy.transpose(data)
                if not flat:
                    # If view is True, we set the particle property as a
                    # view on the dump array. Pay attention of C / F order.
                    # To check if particle properties and dump are associated:
                    # numpy.may_share_memory(p[0].position, pos[:, 0])
                    if view and order == 'C':
                        for i, p in enumerate(self.particle):
                            setattr(p, attr, data[i, ...])
                    if view and order == 'F':
                        for i, p in enumerate(self.particle):
                            setattr(p, attr, data[..., i])
                else:
                    data = data.flatten(order=order)
                    if view and attr in ['position', 'velocity']:
                        ndim = self.number_of_dimensions
                        for i, p in enumerate(self.particle):
                            setattr(p, attr, data[i*ndim: (i+1)*ndim])

            elif what.startswith('cell'):
                data = numpy.array(getattr(self.cell, attr), dtype=dtype)
                if view:
                    setattr(self.cell, attr, data)
            else:
                raise ValueError('Unknown attribute %s' % what)

            self._data[what] = data

        # Always keep track of the number of particles
        self._data['npart'] = len(self.particle)
            
        # If what is a string or we only have one entry we return an
        # array, otherwise we return a dict with the requested keys
        if len(what_list) == 1:
            return self._data[what_list[0]]
        else:
            db = {}
            for key in what_list:
                db[key] = self._data[key]
            return db
                
    def report(self):
        # Summary
        txt = ''
        try:
            if self.particle:
                txt += 'system composed by {0} particles\n'.format(len(self.particle))
            if self.cell:
                txt += 'enclosed in a {0.shape} box at number density rho={1:.6f}\n'.format(self.cell, self.density)
            if self.thermostat:
                txt += 'in contact with a thermostat at T={0.temperature}\n'.format(self.thermostat)
            if self.barostat:
                txt += 'in contact with a barostat at P={0.pressure}\n'.format(self.barostat)
            if self.reservoir:
                txt += 'in contact with a reservoir at mu={0.chemical_potential}\n'.format(self.reservoir)
            if self.interaction:
                txt += '\n'
                txt += self.interaction.report()
        except:
            pass
        return txt

    def show(self, backend='matplotlib', *args, **kwargs):
        from .particle import show_ovito
        from .particle import show_matplotlib
        from .particle import show_3dmol
        if backend == 'matplotlib':
            _show = show_matplotlib
        elif backend == 'ovito':
            _show = show_ovito
        elif backend == '3dmol':
            _show = show_3dmol
        else:
            raise ValueError('unknown backend for visualization')
        return _show(self.particle, self.cell, *args, **kwargs)
