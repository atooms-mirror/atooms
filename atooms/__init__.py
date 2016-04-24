# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

"""
atooms - Atomistic object-oriented modeling and simulations
"""

__author__ = "Daniele Coslovich <daniele.coslovich@um2.fr>"

try:
    from _version import __version__, __date__
except ImportError:
    __version__ = "unknown"
    __date__ = "unknown"

ndim = 3
