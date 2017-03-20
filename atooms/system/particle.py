# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Point particles in a cartesian reference frame."""

import numpy
import random
import copy
from atooms.core import ndim


class Particle(object):

    def __init__(self, id=1, name='A', mass=1.0,
                 position=numpy.zeros(ndim),
                 velocity=numpy.zeros(ndim), radius=0.5, tag=None):
        self.id = id
        """An integer chemical id of the particle."""
        self.name = name
        """The name of the chemical species of the particle."""
        self.mass = mass
        self.radius = radius
        self.position = numpy.asarray(position)
        self.velocity = numpy.asarray(velocity)
        self.tag = tag

    def nearest_image(self, particle, cell):
        """
        Transform self into the nearest image of `particle` in the
        specified cell.
        """
        rij = self.position - particle.position
        periodic_vector(rij, cell.side)
        self.position = particle.position + rij
        return self

    def nearest_image_copy(self, particle, cell):
        """Return the nearest image of `particle` in the specified `cell`."""
        from copy import deepcopy
        rij = self.position - particle.position
        periodic_vector(rij, cell.side)
        image = deepcopy(self)
        image.position = particle.position + rij
        return image

    def distance(self, particle, cell=None):
        """
        Return distance from `particle`.

        If `cell` is provided, return distance from the nearest image
        of `particle`.
        """
        r = self.position - particle.position
        if cell:
            periodic_vector(r, cell.side)
        return r

    def fold(self, cell):
        """Fold self into central cell."""
        self.position = periodic_vector_safe(self.position, cell.side)
        return self

    @property
    def diameter(self):
        """Particle diameter."""
        return self.radius * 2

    def maxwellian(self, T):
        """
        Assign the velocity to particle according to a Maxwell-Boltzmann
        distribution at temperature `T`.
        """
        vx = random.gauss(0, numpy.sqrt(T / self.mass))
        vy = random.gauss(0, numpy.sqrt(T / self.mass))
        vz = random.gauss(0, numpy.sqrt(T / self.mass))
        self.velocity = numpy.array((vx, vy, vz))


# Utility functions

def periodic_vector(vec, box):
    for i in xrange(vec.shape[0]):
        if vec[i] > box[i] / 2:
            vec[i] += - box[i]
        elif vec[i] < -box[i] / 2:
            vec[i] += box[i]
    # return numpy.where(abs(a) > box/2, a-numpy.copysign(box, a), a)
    return vec

def periodic_vector_safe(vec, box):
    return vec - numpy.rint(vec / box) * box

def periodic_vector_safe_opti(vec, box, invbox):
    return vec - numpy.rint(vec * invbox) * box

def fix_cm(particles):
    """
    Subtract out the center of mass velocity from a list of
    `particles`.
    """
    vcm = velocity_cm(particles)
    for p in particles:
        p.velocity -= vcm
    return particles

def velocity_cm(particle):
    """Velocity of the center of mass of a list of particles."""
    vcm = numpy.zeros_like(particle[0].velocity)
    mtot = 0.0
    for p in particle:
        vcm += p.velocity * p.mass
        mtot += p.mass
    return vcm / mtot

def position_cm(particle):
    """Position of the center of mass of a list of particles."""
    rcm = numpy.zeros_like(particle[0].position)
    mtot = 0.0
    for p in particle:
        rcm += p.position * p.mass
        mtot += p.mass
    return rcm / mtot

def total_kinetic_energy(particles):
    """Total kinetic energy of a list of `particles`."""
    ekin = 0.0
    for p in particles:
        ekin += p.mass * numpy.dot(p.velocity, p.velocity)
    return 0.5 * ekin

def temperature(particles, ndof=None):
    """Kinetic temperature of a list of `particles`."""
    if ndof is None:
        ndof = ndim * (len(particles) - 1)
    return 2 * total_kinetic_energy(particles) / ndof

def species(particles):
    """Return list of distinct species (`id`) of `particles`."""
    return list(set([p.id for p in particles]))

def composition(particles, nsp=None):
    """
    Return a tuple containing the number of particles of each species
    appearing the input `particles` list.
    """
    # TODO: check id normalization
    if nsp is None:
        x = max([p.id for p in particles]) * [0]
    else:
        x = nsp * [0]
    for p in particles:
        x[p.id - 1] += 1
    return tuple(x)

def rotated(particle, cell):
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
