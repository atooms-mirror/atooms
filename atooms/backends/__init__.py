# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Simulation backends.

A backend is an adapter to an actual implementation of a
simulation. The simulation backend interface must provide, at a
minimum, a `run()` method that performs the simulation. See the
`DryRun` backend for the interface that a valid backend should
implement.
"""
