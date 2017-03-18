# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Simulation package.
"""

import logging
from atooms.utils import NullHandler
logging.getLogger(__name__).addHandler(NullHandler())

from .base import *
