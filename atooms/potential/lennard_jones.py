# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

from potential import PairPotentialBase

class LennardJones(PairPotentialBase):

    def compute(self, rsq):
        sigsq = self.params['sigma']**2
        u = 4 * self.params['eps'] * ((sigsq/rsq)**6 - (sigsq/rsq)**3)
        a = 24 * self.params['eps'] * (2*(sigsq/rsq)**6 - (sigsq/rsq)**3) / rsq
        return u, a
            
# from cutoff import CutOff

# p = LennardJones('LJ', {'eps':1.0, 'sigma':1.0}, (1, 1), cutoff=CutOff('cs', 2.5), npoints=100)
# rsq, u0, u1 = p.tabulate()
# for x, y, z in zip(rsq, u0, u1):
#     print x, y, z
