# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich


from potential import PairPotentialBase

class LennardJones(PairPotentialBase):

    def compute(self, rsquare):
        return 4 * self.params["eps"] * ((self.params["sigma"]**2/rsquare)**6 -
                                         (self.params["sigma"]**2/rsquare)**3 )
