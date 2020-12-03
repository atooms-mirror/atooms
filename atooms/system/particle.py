# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""Point particles in a cartesian reference frame."""

import logging
import numpy
import random
from copy import deepcopy
from atooms.core import ndim as _ndim

_log = logging.getLogger(__name__)


class Particle(object):

    def __init__(self, species='A', mass=1.0, position=None,
                 velocity=None, radius=0.5):
        # TODO: position should have been first
        self.species = species
        """The chemical species of the particle."""
        self.mass = mass
        self.radius = radius
        if position is None:
            self.position = numpy.zeros(_ndim)
        else:
            self.position = numpy.asarray(position)
        if velocity is None:
            self.velocity = numpy.zeros(_ndim)
        else:
            self.velocity = numpy.asarray(velocity)

    @property
    def diameter(self):
        """Particle diameter."""
        return self.radius * 2

    def nearest_image(self, particle, cell, copy=False, folded=False):
        """
        Return the nearest image of `particle` in the given `cell`.

        If `copy` is `False`, the particle is transformed into to its
        nearest image, otherwise the fucction returns a copy of the
        nearest image particle and leave the original particle as is.

        If `folded` is True, the coordinates are assumed to be folded
        into the central cell (or in a cell just next to it),
        otherwise they will be assumed to lie in an arbitrary periodic
        cell.
        """
        rij = self.position - particle.position
        if folded:
            rij = _periodic_vector(rij, cell.side)
        else:
            rij = _periodic_vector_unfolded(rij, cell.side)
        if copy:
            image = deepcopy(self)
            image.position = particle.position + rij
            return image
        else:
            self.position[:] = particle.position + rij
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
                r = _periodic_vector_unfolded(r, cell.side)
        return r

    def __repr__(self):
        return 'Particle(species={0.species}, mass={0.mass}, ' \
            'position={0.position}, velocity={0.velocity}, ' \
            'radius={0.radius})'.format(self)

    def fold(self, cell):
        """Fold self into central cell."""

        # Move the center to 0
        self.position -= cell.center
        self.position[:] = _periodic_vector_unfolded(self.position, cell.side)

        # Restore the center
        self.position += cell.center
        return self

    def maxwellian(self, T):
        """
        Assign the velocity to particle according to a Maxwell-Boltzmann
        distribution at temperature `T`.
        """
        vx = random.gauss(0, numpy.sqrt(T / self.mass))
        vy = random.gauss(0, numpy.sqrt(T / self.mass))
        vz = random.gauss(0, numpy.sqrt(T / self.mass))
        self.velocity[:] = numpy.array([vx, vy, vz])

    @property
    def kinetic_energy(self):
        """Kinetic energy."""
        return 0.5 * self.mass * numpy.dot(self.velocity, self.velocity)


# Utility functions

def _periodic_vector(vec, box):
    # TODO: what about particle distances precisely equal to L/2 or -L/2?
    for i in range(vec.shape[0]):
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


def _periodic_vector_delta_unfolded(vec, box):
    """
    Return periodic vector delta to enable in-place modification of
    folded particles.
    """
    # Optimized version
    # return numpy.rint(vec * invbox) * box
    return numpy.rint(vec / box) * box


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
    if len(particle) == 0:
        return 0.0
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

    p = deepcopy(particle)
    dist = [sum(p[0].distance(pi, cell)**2) for pi in p]
    dmax = max(dist)
    imax = dist.index(dmax)
    z_axis = numpy.array([0, -1, 0])
    pr_axis = p[0].position - p[imax].position
    ro_axis = numpy.cross(pr_axis, z_axis)
    theta = math.acos(numpy.dot(pr_axis, z_axis) / (norm2(pr_axis) * norm2(z_axis))**0.5)
    for pi in p:
        pi.position = numpy.dot(rotation(ro_axis, theta), pi.position)
    return p


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


def gyration_radius(particles, cell=None, weight=None, center=None,
                    method="N1"):
    """
    Gyration radius of a list of `particles` in a `cell` with an
    optional `weight`.

    The optional `center` variable is the index of the central
    particle for pbc unfolding. If `weight` is given but `center` is
    not, the latter will be the index of the first largest element of
    `weight`.

    Possible values of `methods` are `N1` (default), `N2` and `min`.
    """
    if method not in ['N1', 'N2', 'min']:
        raise ValueError('unknown method %s' % method)

    # Order N^2 calculation
    if method == 'N2':
        rg = 0.0
        for pi in particles:
            for pj in particles:
                if pi is pj:
                    continue
                dr = pi.distance(pj, cell)
                rg += numpy.dot(dr, dr)
        return (rg / (2*len(particles)**2))**0.5

    # Order N^1 calculation
    elif method == 'N1':
        # Assign weights
        if weight is None:
            weight = numpy.ones(len(particles))
            if center is None:
                center = 0
        else:
            if len(weight) != len(particles):
                raise ValueError('n. weights differs from n. particles')
            weight = numpy.asarray(weight)
            if center is None:
                center = weight.argmax()

        # Unfold cluster across pbc
        p_central = particles[center]
        if cell is None:
            cluster = particles
        else:
            cluster = []
            for p in particles:
                cluster.append(p.nearest_image(p_central, cell, copy=True))

        def weighted_cm_position(particle, weight):
            """Weighted center-of-mass of a list of particles."""
            rcm = numpy.zeros_like(particle[0].position)
            wtot = sum(weight)
            for i, p in enumerate(particle):
                rcm += p.position * weight[i]
            return rcm / wtot

        # Compute gyration radius
        rcm = weighted_cm_position(cluster, weight)
        rg = 0.0
        for i, p in enumerate(cluster):
            dr = p.position - rcm
            rg += numpy.dot(dr, dr) * weight[i]
        rg /= len(cluster)
        return rg**0.5

    # Minimize over possible centers
    elif method == 'min':
        rg = None
        for center in range(len(particles)):
            rg_new = gyration_radius(particles, cell=cell,
                                     weight=weight, center=center,
                                     method="N1")
            if rg is None:
                rg = rg_new
            else:
                rg = min(rg_new, rg)
        return rg


def collective_overlap(particle, other, a, side, normalize=True):
    """Compute collective overlap between two lists of particles."""
    # These optimizations via numpy are necessary. Computing O(N^2)
    # quantities in pure python takes forever.
    x = numpy.array([p.position for p in particle])
    y = numpy.array([p.position for p in other])
    dr = numpy.ndarray(x.shape[0])
    rij = numpy.asarray(y)
    q = 0
    for i in range(y.shape[0]):
        dr = x[:, :] - y[i, :]
        dr = dr - numpy.rint(dr / side) * side
        dr = numpy.sum(dr**2, axis=1)
        q += (dr < a**2).sum()
    if normalize:
        q /= float(len(particle))
    return q


def self_overlap(particle, other, a, normalize=True):
    """Compute self overlap between two lists of unfolded particles."""
    # This quantity is O(N) and therefore we use pure python
    rij = [sum(p.distance(o)**2) for p, o in zip(particle, other)]
    q = (numpy.array(rij) < a**2).sum()
    if normalize:
        q /= float(len(particle))
    return q


def show_3dmol(particle, cell=None, radius=1.0, palette=None):
    """
    Render particles in cell using 3dmol http://3dmol.csb.pitt.edu/
    """
    import py3Dmol
    if palette is None:
        palette = ["#50514f", "#f25f5c", "#ffe066", "#247ba0", "#70c1b3",
                   "#0cce6b", "#c200fb", "#e2a0ff", "#6622cc", "#119822"]
    colors = {}
    for i, s in enumerate(distinct_species(particle)):
        colors[s] = palette[i]
    view = py3Dmol.view()
    view.setBackgroundColor('white')
    for p in particle:
        view.addSphere({'center': {'x': p.position[0], 'y': p.position[1], 'z': p.position[2]},
                        'radius': radius * p.radius, 'color': colors[p.species]})
    if cell is not None:
        view.addBox({'center': {'x': cell.center[0], 'y': cell.center[1], 'z': cell.center[2]},
                     'dimensions': {'w': cell.side[0], 'h': cell.side[1], 'd': cell.side[2]},
                     'wireframe': True, 'color': "#000000"})
    return view


def show_matplotlib(particle, cell, outfile=None, linewidth=3, alpha=0.3, show=False):
    """
    Make a snapshot of the `particle`s in the `cell` and save the
    image in `outfile`. The image is returned for further
    customization.
    """
    import matplotlib.pyplot as plt
    from .particle import distinct_species

    color_db = ['b', 'r', 'g', 'y']
    species = distinct_species(particle)
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlim((-cell.side[0]/2, cell.side[0]/2))
    ax.set_ylim((-cell.side[1]/2, cell.side[1]/2))
    for p in particle:
        c = plt.Circle(p.position[: 2], p.radius,
                       facecolor=color_db[species.index(p.species)],
                       edgecolor='black', alpha=alpha, linewidth=linewidth)
        ax.add_artist(c)
    if outfile is not None:
        fig.savefig(outfile, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def show_ovito(particle, cell, outfile=None, radius=0.35,
               viewport=None, callback=None, tmpdir=None,
               camera_dir=(0, 1, 0), camera_pos=(0, -10, 0),
               size=(640, 480), zoom=True, perspective=False):
    """
    Render particle in cell using ovito
    """
    import os
    try:
        from ovito.io import import_file
    except ImportError:
        _log.warning('install ovito to display the particles')
        return
    from ovito.vis import Viewport, TachyonRenderer
    from ovito.vis import ParticlesVis
    import tempfile
    from atooms.core.utils import mkdir

    # Make sure dirname exists
    if outfile is not None:
        mkdir(os.path.dirname(outfile))

    # Get a temporary file to write the sample
    fh = tempfile.NamedTemporaryFile('w', dir=tmpdir, suffix='.xyz', delete=False)
    tmp_file = fh.name

    # Self-contained EXYZ dump (it is not clean to use trajectories here)
    fh.write('{}\n'.format(len(particle)))
    fh.write('Properties=species:S:1:pos:R:3 Lattice="{},0.,0.,0.,{},0.,0.,0.,{}"\n'.format(*cell.side))
    for p in particle:
        fh.write('{} {} {} {}\n'.format(p.species, *p.position))
    fh.close()

    # Ovito stuff. Can be customized by client code.
    pipeline = import_file(tmp_file)
    # Ovito seems to ignore the lattice info of exyz file
    # so we forcibly set the cell info here
    pipeline.source.data.cell_[0, 0] = cell.side[0]
    pipeline.source.data.cell_[1, 1] = cell.side[1]
    pipeline.source.data.cell_[2, 2] = cell.side[2]
    pipeline.source.data.cell_[:, 3] = -cell.side/2
    # Scale radius by default
    vis_element = pipeline.source.data.particles.vis
    vis_element.radius = radius
    # Apply client code callback
    if callback:
        callback(pipeline)
    pipeline.add_to_scene()

    # Define viewport
    if viewport:
        vp = vieport
    else:
        if perspective:
            vp = Viewport(type=Viewport.Type.Perspective, camera_dir=camera_dir, camera_pos=camera_pos)
        else:
            vp = Viewport(type=Viewport.Type.Ortho, camera_dir=camera_dir, camera_pos=camera_pos)

    # Render
    if zoom:
        vp.zoom_all()
    if outfile is None:
        outfile = tmp_file + '.png'

    vp.render_image(filename=outfile,
                    size=size,
                    renderer=TachyonRenderer())

    # Scene is a singleton, so we must clear it
    pipeline.remove_from_scene()

    from atooms.core.utils import rmf
    rmf(tmp_file)

    # Try to display the image (e.g. in a jupyter notebook)
    try:
        from IPython.display import Image
        return Image(outfile)
    except ImportError:
        return outfile

show = show_matplotlib


def decimate(particle, N):
    """
    Return a decimated list of N particles keep the same chemical
    composition as in the input list `particle`.
    """
    import random
    if N > len(particle):
        raise ValueError('cannot increase number of particles')
    x = composition(particle)
    scale = float(N) / len(particle)
    idx = {}
    for species in x:
        idx[species] = random.sample(range(x[species]), int(x[species] * scale))

    # if sum([len(idx[species]) for species in x]) != N:
    #     print sum([len(idx[species]) for species in x]), N,  [len(idx[species]) for species in x], x
    #     raise ValueError('cannot preserve composition')

    # Go through the particles once and map sequential indices to the subset of indices of a given species
    from collections import defaultdict
    idx_species = defaultdict(list)
    for i, p in enumerate(particle):
        idx_species[p.species].append(i)

    # Collect selected particles
    pnew = []
    for species in sorted(x):
        for i in idx[species]:
            ii = idx_species[species][i]
            pnew.append(particle[ii])

    return pnew
