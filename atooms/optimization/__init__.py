# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Optimization framework for systems of interacting particles.

`atooms` provides a generic optimization interface that abstracts out
most of the common parts of particle-based optimizations.
"""

import logging
from .core import Optimization
from atooms.core.utils import NullHandler
logging.getLogger(__name__).addHandler(NullHandler())
