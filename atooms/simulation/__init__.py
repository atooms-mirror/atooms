# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Simulation framework for particulate systems. 

`atooms` provides a generic simulation interface that abstracts out
most of the common parts of particle-based simulations. It uses
callbacks to analyze and process simulation data on the fly.
"""

import logging
from atooms.utils import NullHandler
logging.getLogger(__name__).addHandler(NullHandler())

from .base import *
