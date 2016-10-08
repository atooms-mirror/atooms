# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich


class CutOff(object):
    
    def __init__(self, name, radius):
        self.name = name # this is what was called scheme
        self.radius = radius
        self.rcutsq = radius**2
        
    @property
    def effective_radius(self):
        return self.radius

    def is_zero(self, rsquare):
        return rsquare > self.radius**2

    def tailor(self, rsquare):
        raise NotImplementedError()

    def smooth(self, rsquare):
        raise NotImplementedError()
