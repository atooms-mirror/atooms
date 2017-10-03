# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""Point particles in a cartesian reference frame."""

import numpy
import random
import copy
from atooms.core import ndim as _ndim


class Particle(object):

    def __init__(self, species='A', mass=1.0,
                 position=numpy.zeros(_ndim),
                 velocity=numpy.zeros(_ndim), radius=0.5):
        self.species = species
        """The chemical species of the particle."""
        self.mass = mass
        self.radius = radius
        self.position = numpy.asarray(position)
        self.velocity = numpy.asarray(velocity)

    @property
    def diameter(self):
        """Particle diameter."""
        return self.radius * 2

    def nearest_image(self, particle, cell, copy=False):
        """
        Return the nearest image of `particle` in the given `cell`.

        If `copy` is `False`, the particle is transformed into to its
        nearest image, otherwise the fucction returns a copy of the
        nearest image particle and leave the original particle as is.
        """
        rij = self.position - particle.position
        _periodic_vector(rij, cell.side)
        if copy:
            from copy import deepcopy
            image = deepcopy(self)
            image.position = particle.position + rij
            return image
        else:
            self.position = particle.position + rij
            return self

    def distance(self, particle, cell=None, folded=True):
        """
        Return the distance from another `particle`.

        If `cell` is given, return the distance from the nearest image
        of `particle`. In this case, if `folded` is True, the
        coordinates are assumed to be folded into the central cell (or
        in a cell just next to it), otherwise they will be assumed to
        lie in an arbitrary periodic cell.

        If `cell` is None, periodic boundary conditions are *not*
        applied and folded is irrelevant.
        """
        r = self.position - particle.position
        if cell is not None:
            # Apply periodic boundary conditions
            if folded:
                r = _periodic_vector(r, cell.side)
            else:
                r = _periodic_vector_unfolded(r, cell.sidebox)
        return r

    def fold(self, cell):
        """Fold self into central cell."""
        self.position = _periodic_vector_unfolded(self.position, cell.side)
        return self

    def maxwellian(self, T):
        """
        Assign the velocity to particle according to a Maxwell-Boltzmann
        distribution at temperature `T`.
        """
        vx = random.gauss(0, numpy.sqrt(T / self.mass))
        vy = random.gauss(0, numpy.sqrt(T / self.mass))
        vz = random.gauss(0, numpy.sqrt(T / self.mass))
        self.velocity = numpy.array((vx, vy, vz))

    @property
    def kinetic_energy(self):
        """Kinetic energy."""
        return 0.5 * self.mass * numpy.dot(self.velocity, self.velocity)


# Utility functions

def _periodic_vector(vec, box):
    for i in xrange(vec.shape[0]):
        if vec[i] > box[i] / 2:
            vec[i] += - box[i]
        elif vec[i] < -box[i] / 2:
            vec[i] += box[i]
    # Compact version
    # return numpy.where(abs(a) > box/2, a-numpy.copysign(box, a), a)
    return vec

def _periodic_vector_unfolded(vec, box):
    return vec - numpy.rint(vec / box) * box
    # Optimized version
    # return vec - numpy.rint(vec * invbox) * box

def fix_total_momentum(particles):
    """
    Subtract out the center of mass velocity from a list of
    `particles`.
    """
    vcm = cm_velocity(particles)
    for p in particles:
        p.velocity -= vcm
    return particles

def cm_velocity(particle):
    """Velocity of the center of mass of a list of particles."""
    vcm = numpy.zeros_like(particle[0].velocity)
    mtot = 0.0
    for p in particle:
        vcm += p.velocity * p.mass
        mtot += p.mass
    return vcm / mtot

def cm_position(particle):
    """Center-of-mass of a list of particles."""
    rcm = numpy.zeros_like(particle[0].position)
    mtot = 0.0
    for p in particle:
        rcm += p.position * p.mass
        mtot += p.mass
    return rcm / mtot

def distinct_species(particles):
    """Return sorted list of distinct `species` of `particles`."""
    return list(sorted(set([p.species for p in particles])))

def composition(particles):
    """
    Return a dictionary containing the number of particles of each
    species appearing the input `particles` list.
    """
    from collections import defaultdict
    comp = defaultdict(int)
    for p in particles:
        comp[p.species] += 1
    return comp

def rotate(particle, cell):
    """
    Return a list of particles rotated around the main symmetry axis
    of the set.
    """
    import math

    def norm2(vec):
        return numpy.dot(vec, vec)

    def rotation(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = numpy.asarray(axis)
        theta = numpy.asarray(theta)
        axis = axis / math.sqrt(numpy.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return numpy.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    p = copy.deepcopy(particle)
    dist = [sum(p[0].distance(pi, cell)**2) for pi in p]
    dmax = max(dist)
    imax = dist.index(dmax)
    z_axis = numpy.array([0, -1, 0])
    pr_axis = p[0].position - p[imax].position
    ro_axis = numpy.cross(pr_axis, z_axis)
    theta = math.acos(numpy.dot(pr_axis, z_axis) / (norm2(pr_axis) * norm2(z_axis))**0.5)
    for pi in p:
        pi.position = numpy.dot(rotation(ro_axis, theta), pi.position)


def overlaps(particle, cell):
    """Check presence of overlaps between particles."""
    x = []
    for i, pi in enumerate(particle):
        for j, pj in enumerate(particle):
            if j <= i:
                continue
            d = sum(pi.distance(pj, cell)**2)**0.5
            if d < (pi.radius + pj.radius):
                x.append((i, j))
    return len(x) > 0, x
