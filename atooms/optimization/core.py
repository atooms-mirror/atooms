# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Base optimization class.

`atooms` provides a generic optimization interface that abstracts out
most of the common parts of particle-based optimizations.

The actual optimization code is wrapped by a simulation "backend" that
exposes a minimal but coherent interface.
"""

import os
import time
import datetime
import logging

import atooms.core.progress
from atooms.core import __version__
from atooms.core.utils import mkdir, barrier

_log = logging.getLogger(__name__)


# Default exceptions

class OptimizationEnd(Exception):
    """Raised when an targeter reaches its target."""
    pass

class OptimizationKill(Exception):
    """Raised when a simulation is terminated by SIGTERM."""
    pass


def _report(info, file_handle=None, log_echo=True):
    """
    Log `info` string to default logger at level info.

    Optionally write `info` to `file_handle` if the latter is
    given. Logging is disabled via `log_echo` is False.
    """
    if info is None:
        return

    if log_echo:
        for line in info.split('\n'):
            _log.info(line.strip())

    if file_handle is not None:
        file_handle.write(info)


class Optimization(object):

    """Optimization base class."""

    # A stripped down version of Simulation, removed
    # - checkpoint
    # - callbacks
    # - steps
    # - speedometer

    version = __version__

    def __init__(self, backend, output_path=None, method=None,
                 tolerance=1e-10, max_iterations=None):
        """
        Perform an optimization using the specified `backend` and
        optionally write output to `output_path`. This can be a file
        or directory path.

        Paths: to define output paths we rely on `output_path`, all
        other paths are defined based on it and on its base_path.
        """
        self.backend = backend
        self.output_path = output_path
        self.method = method
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.current_step = 0
        self.initial_step = 0
        self._start_time = time.time()

        # We expect subclasses to keep a ref to the trajectory object
        # self.trajectory used to store configurations
        if hasattr(self.backend, 'trajectory'):
            self.trajectory = self.backend.trajectory
        else:
            self.trajectory = None
        # Make sure the dirname of output_path exists. For instance,
        # if output_path is data/trajectory.xyz, then data/ should
        # exist. This creates the data/ folder and its parents folders.
        if self.output_path is not None:
            mkdir(os.path.dirname(self.output_path))

    @property
    def system(self):
        # Note that setting system as a reference in the instance, like
        #   self.system = self.backend.system
        # is unsafe because this won't follow the backend's system when the
        # latter is reassigned as in
        #   self.backend.system = None
        # So we defined it as a property.
        return self.backend.system

    @system.setter
    def system(self, value):
        self.backend.system = value

    def __str__(self):
        return 'atooms optimization via %s' % self.backend

    @property
    def base_path(self):
        return os.path.splitext(self.output_path)[0]

    def _elapsed_wall_time(self):
        """Elapsed wall time in seconds."""
        return time.time() - self._start_time

    def wall_time(self, per_step=False, per_particle=False):
        """
        Elapsed wall time in seconds.

        Optionally normalized per step and or per particle. It can be
        subclassed by more complex simulation classes.
        """
        norm = 1.0
        # Normalize per particle
        if per_particle:
            if len(self.system.particle) > 0:
                norm *= len(self.system.particle)
            else:
                return float('nan')
        # Normalize per step
        if per_step:
            if self.current_step - self.initial_step > 0:
                norm *= (self.current_step - self.initial_step)
            else:
                return float('nan')
        return self._elapsed_wall_time() / norm

    def run(self):
        """Run the optimization."""
        # Report
        _report(self._info_start())
        _report(self._info_backend())
        if hasattr(self.system, 'report'):
            _report(self.system.report())
        if hasattr(self.backend, 'report'):
            _report(self.backend.report())

        barrier()
        self.initial_step = self.current_step
        self._start_time = time.time()

        _log.info('starting at step: %d', self.current_step)
        _log.info('')

        import signal
        import sys

        def signal_term_handler(signal, frame):
            raise OptimizationKill('simulation terminated')

        signal.signal(signal.SIGTERM, signal_term_handler)

        try:
            self.backend.run()
            raise OptimizationEnd

        except OptimizationEnd as end:
            _log.info('simulation ended successfully: %s', end)
            _report(self._info_end())

        except KeyboardInterrupt:
            pass

        except OptimizationKill:
            _log.info('simulation terminated')

        except:
            _log.error('simulation failed')
            raise

    def _info_start(self):
        now = datetime.datetime.now().strftime('%Y-%m-%d at %H:%M')
        txt = """\

        {}

        version: {}
        atooms version: {}
        optimization started on: {}
        output path: {}\
        """.format(self, self.version, __version__, now, self.output_path)
        return txt

    def _info_backend(self):
        """Subclasses may want to override this method."""
        txt = 'backend: {}\n'.format(self.backend)
        if hasattr(self.backend, 'version'):
            txt += 'backend version: %s\n' % self.backend.version
        return txt

    def _info_observers(self):
        txt = []
        for f in self._callback:
            params = self._cbk_params[f]
            s = params['scheduler']
            if 'target' in _callable_name(f):
                args = params['args']
                txt.append('target %s: %s' % (_callable_name(f), args[0]))
            else:
                txt.append('writer %s: interval=%s calls=%s' %
                           (_callable_name(f), s.interval, s.calls))
        return '\n'.join(txt) + '\n'

    def _info_end(self):
        now = datetime.datetime.now().strftime('%Y-%m-%d at %H:%M')
        txt = """
        wall time [s]: {:.2f}
        average TSP [s/step/particle]: {:.2e}
        optimization ended on: {}\
        """.format(self.wall_time(), self.wall_time(per_step=True,
                                                    per_particle=True),
                   now)
        return txt
