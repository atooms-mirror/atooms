# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Interactions between particles.

Individual particles interact via a `Potential`. A simple example is a
two-body potential that depends only on the scalar distance between
two particles. The potentials can be cut off and smoothed by a
`CutOff`.

`Interaction` accounts for the total interaction of all the
particles in a system (order N^2 calculation). Valid interaction
backends must implement the interface of `Interaction`.
"""

from .interaction import Interaction
from .potential import PairPotential
from .cutoff import CutOff
