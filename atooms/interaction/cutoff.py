# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Cut-off schemes for pair potentials.

The potential is truncated and optionally smoothed at the cut-off
distance `radius`.
"""

class CutOff(object):
    
    def __init__(self, name, radius):
        """
        Available schemes are
        
        - cut: `c` or `cut` 
        - cut and shifted: `cs` or `CS`
        """
        # TODO: call it scheme again, not name
        self.name = name
        self.radius = radius
        self.rcutsq = radius**2
    
    def __str__(self):
        return self.name

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
        if self.name in ['cs', 'CS']:
            self.vcut = u[0]
        elif self.name in ['c', 'cut']:
            pass
        else:
            raise NotImplementedError()

    def smooth(self, rsquare, u):
        """
        Smooth the potential.

        This must be called by `atooms.interaction.PairPotential` in
        `atooms.interaction.PairPotential.compute`.
        """
        u_new = list(u)
        if self.name in ['cs', 'CS']:
            u_new[0] = u[0] - self.vcut
        elif self.name in ['c', 'cut']:
            pass
        else:
            raise NotImplementedError()
        return u_new
