# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Cut-off schemes for pair potentials.

The potential is truncated and optionally smoothed at the cut-off
distance `radius`.
"""

_db = {'c': 'cut',
        'cs': 'cut and shifted',
        'cspl': 'cut and cubic splined',
        'qs': 'cut and quadratically shifted'}

class CutOff(object):

    def __init__(self, scheme, radius):
        """
        Available schemes:

        - cut: `c` or `cut`
        - cut and shifted: `cs` or `CS`
        """
        self.scheme = scheme
        self.radius = radius
        self.rcutsq = radius**2
        self.radius_mid = radius
        self.radius_mid_sq = radius**2

    def __str__(self):
        return _db[self.scheme]

    @property
    def effective_radius(self):
        return self.radius

    def is_zero(self, rsquare):
        """Return true if `rsquare` is beyond cutoff `radius`"""
        return rsquare > self.radius**2

    def tailor(self, rsquare, u):
        """
        Adjust the cut off to a potential.

        `u` is tuple with the potential value and its derivatives
        evaluated at the cut off distance.
        """
        if self.scheme in ['cs', 'CS']:
            self.vcut = u[0]

        elif self.scheme in ['c', 'cut']:
            pass

        elif self.scheme in ['qs', 'QS']:
            self.Acut =   u[1] / 2
            self.Bcut = - u[1] * self.radius**2 / 2 - u[0]

        elif self.scheme in ['cspl', 'CSPL']:
            # the potential is cubic-splined between rcut1 < r < rcut
            # cubic-splined cutoff : u(r) = A*(B-r)^3 + C for rcut1 < r < rcut=B
            # where
            #   A = - u2^2/(12*u1)
            #   B = rcut1 - 2*(u1/u2) = rcut
            #   C = u0 + (2*u1^2)/(3*u2)
            # with u0=u(rcut1), u1=u'(rcut1), u2=u''(rcut1)
            u2 = u[2]*self.radius_mid**2 - u[1]  # second order derivative
            u1 = - u[1] * self.radius_mid  # first order derivative
            u0 = u[0]
            self.Acut = - u2**2 / (12.0*u1)
            self.Ccut = - u0 + (2.0*u1**2) / (3.0*u2)
            self.radius = self.radius_mid - (2.0*u1/u2)
            self.radius_sq = self.radius**2
        else:
            raise NotImplementedError()

    def smooth(self, rsquare, u):
        """
        Smooth the potential.

        This must be called by `atooms.interaction.PairPotential` in
        `atooms.interaction.PairPotential.compute`.
        """
        u_new = list(u)
        if self.scheme in ['cs', 'CS']:
            u_new[0] = u[0] - self.vcut

        elif self.scheme in ['c', 'cut']:
            pass

        elif self.scheme in ['qs', 'QS']:
            u_new[0] = u[0] + self.Acut * rsquare + self.Bcut
            u_new[1] = u[1] - self.Acut * 2

        elif self.scheme in ['cspl', 'CSPL']:
            # cubic splined between rcut1 and rcut :
            # function is overwritten when between rcut1 and rcut
            if rsquare > self.radius_mid_sq:
                rij = rsquare**0.5
                dr  = self.radius - rij
                u_new[0] = self.Acut * dr**3
                u_new[1] = 3 * self.Acut * dr**2 / rij
            else:
                u_new[0] = u[0] + self.Ccut

        else:
            raise NotImplementedError()
        return u_new
