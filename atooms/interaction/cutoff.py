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

    def tailor(self, rsquare, u):
        if self.name in ['cs', 'CS']:
            self.vcut = u[0]
        elif self.name in ['c', 'cut']:
            pass
        else:
            raise NotImplementedError()

    def smooth(self, rsquare, u):
        if self.name in ['cs', 'CS']:
            u[0] = u[0] - self.vcut
        elif self.name in ['c', 'cut']:
            pass
        else:
            raise NotImplementedError()
        return u
