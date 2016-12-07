# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

from potential import PairPotentialBase

class LennardJones(PairPotentialBase):

    def _compute(self, rsq):
        if 'eps' in self.params:
            self.params['epsilon'] = self.params['eps']
        sigsq = self.params['sigma']**2
        u = 4 * self.params['epsilon'] * ((sigsq/rsq)**6 - (sigsq/rsq)**3)
        a = 24 * self.params['epsilon'] * (2*(sigsq/rsq)**6 - (sigsq/rsq)**3) / rsq
        return u, a
