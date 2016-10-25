# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import numpy 
import random
from atooms.core import ndim

def periodic_vector(a, box):
    #return numpy.where(abs(a) > box/2, a-numpy.copysign(box, a), a)
    for i in xrange(a.shape[0]):
        if a[i] > box[i]/2:
            a[i] += - box[i]
        elif a[i] < -box[i]/2:
            a[i] += box[i]
    return a

def periodic_vector_safe(a, box):
    return a - numpy.rint(a/box) * box

def periodic_vector_safe_opti(a, box, invbox):
    return a - numpy.rint(a * invbox) * box
    
def fix_cm(particle):
    """ Subtract out the motion of the CM """
    vcm = velocity_cm(particle)
    for p in particle:
        p.velocity -= vcm
    return p

def velocity_cm(particle):
    """ Velocity of the center of mass of a list of particles """
    vcm = numpy.zeros_like(particle[0].velocity)
    mtot = 0.0
    for p in particle:
        vcm += p.velocity * p.mass
        mtot += p.mass
    return vcm / mtot

def position_cm(particle):
    """ Position of the center of mass of a list of particles """
    rcm = numpy.zeros_like(particle[0].position)
    mtot = 0.0
    for p in particle:
        rcm += p.position * p.mass
        mtot += p.mass
    return rcm / mtot

def total_kinetic_energy(particle):
    ekin = 0.0
    for p in particle:
        ekin += p.mass * numpy.dot(p.velocity, p.velocity)
    return 0.5 * ekin

def temperature(particle, ndof=None):
    if ndof is None:
        ndof = ndim * (len(particle) - 1)
    return 2 * total_kinetic_energy(particle) / ndof

def composition(particle):
    """Return a tuple containing the number of particles 
    of each species appearing the input particle list"""
    # TODO: check id normalization
    x = max([p.id for p in particle]) * [0]
    for p in particle:
        x[p.id-1] += 1
    return tuple(x)

def rotated(particle, cell):
    """Return a rotated list of particles around the main symmetry axis of the set."""
    import numpy as np
    import math

    def norm2(x):
        return numpy.dot(x, x)

    def rotation(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        theta = np.asarray(theta)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

    p = copy.deepcopy(particle)
    dist = [sum(p[0].distance(pi, cell)**2) for pi in p]
    dmax = max(dist)
    imax = dist.index(dmax)
    z_axis = numpy.array([0,-1,0])
    pr_axis = p[0].position - p[imax].position
    ro_axis = numpy.cross(pr_axis, z_axis)
    theta = math.acos(numpy.dot(pr_axis, z_axis) / (norm2(pr_axis) * norm2(z_axis))**0.5)
    for pi in p:
        pi.position = numpy.dot(rotation(ro_axis, theta), pi.position)            

class Particle(object):

    """
    Particle class
    """

    def __init__(self,
                 id       = 1,
                 name     = 'A',
                 mass     = 1.0,
                 position = numpy.zeros(ndim),
                 velocity = numpy.zeros(ndim),                 
                 radius   = 0.5, # sigma=1.0
                 tag      = None):
        self.id       = id
        self.name     = name
        self.mass     = mass
        self.radius   = radius
        self.position = position
        self.velocity = velocity
        self.tag      = tag
    
    def move(self):
        pass

    def nearest_image(self, p, cell):
        """Transform self into the nearest image of p in the specified cell."""
        rij = self.position - p.position
        periodic_vector(rij, cell.side)
        self.position = p.position + rij
        return self

    def nearest_image_copy(self, p, cell):
        """Return the nearest image of p in the specified cell."""
        from copy import deepcopy
        rij = self.position - p.position
        periodic_vector(rij, cell.side)
        image = deepcopy(self)
        image.position = p.position + rij
        return image

    def distance(self, p, cell=None):
        """
        Return distance from particle p. 
        If cell is provided compute distance from periodic image of p
        """
        r = self.position - p.position
        if cell:
            periodic_vector(r, cell.side)
        return r

    def fold(self, cell):
        """Fold self into central cell."""
        self.position = periodic_vector_safe(self.position, cell.side)
        return self

    def maxwellian(self, T):
        vx = random.gauss(0, numpy.sqrt(T / self.mass))
        vy = random.gauss(0, numpy.sqrt(T / self.mass))
        vz = random.gauss(0, numpy.sqrt(T / self.mass))
        self.velocity = numpy.array((vx, vy, vz))

if __name__ == '__main__':
    particle = Particle()
    for key in  particle.__dict__.keys():
        print key
