# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Base optimization class.

`atooms` provides a generic optimization interface that abstracts out
most of the common parts of particle-based optimizations.

The actual optimization code is wrapped by a simulation "backend" that
exposes a minimal but coherent interface.

`Optimization` uses callbacks to analyze and process optimization data
on the fly. The interval in steps at which callbacks are executed is
controlled by a `Scheduler` instance.
"""

import sys
import logging

from atooms.simulation import Simulation
from atooms.simulation.observers import SimulationEnd

_log = logging.getLogger(__name__)


class OptimizationEnd(SimulationEnd):
    """Raised when an targeter reaches its target."""
    pass

def _target_backend_tolerance(sim, value):
    if not sim.backend.reached_steps:
        raise OptimizationEnd('reached backend tolerance')

def _target_force_norm_square(sim, value):
    if sim.system.force_norm_square(per_particle=True) < value:
        raise OptimizationEnd('reached target W %d' % value)


class Optimization(Simulation):

    """Optimization base class."""

    def __init__(self, backend, output_path=None, tolerance=None,
                 steps=sys.maxint):
        """
        Perform an optimization using the specified `backend` and
        optionally write output to `output_path`. This can be a file
        or directory path.

        Paths: to define output paths we rely on `output_path`, all
        other paths are defined based on it and on its base_path.
        """
        Simulation.__init__(self, backend, output_path=output_path)
        self.steps = steps
        self.tolerance = tolerance
        self._check_interval = 1000

    def __str__(self):
        return 'atooms optimization via %s' % self.backend

    def write_checkpoint(self):
        """Write checkpoint to allow restarting a simulation."""
        pass

    def read_checkpoint(self):
        """Read the checkpoint to restart a simulation."""
        pass

    def run(self, steps=None):
        if self.tolerance is not None:
            # We need to add this check only if targeters are present
            self.add(_target_backend_tolerance, self._check_interval, 0.0)
            self.add(_target_force_norm_square, self._check_interval, self.tolerance)
        Simulation.run(self, steps)
        
