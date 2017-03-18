# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Base simulation class with callback logic.

`atooms` provides a generic simulation interface that abstracts out
most of the common parts of particle-based simulations.

`Simulation` uses callbacks to analyze and process simulation data on
the fly. The module `atooms.simulation.observers` provides basic
callbacks to write data to disk, e.g. trajejectory files, and to stop
the simulation when certain targets are reached, e.g. mean squared
displacement larger than a threshold.

The interval in steps at which callbacks are executed is controlled by
a `Scheduler` instance.

The actual simulation code is wrapped by a simulation "backend" that
exposes a minimal but coherent interface.
"""

import os
import time
import datetime
import logging

from atooms.core import __version__, __commit__, __date__
from atooms.utils import mkdir, barrier
from .backends import DryRunBackend
from .observers import TargetSteps, Target, Speedometer, Scheduler
from .observers import SimulationEnd

log = logging.getLogger(__name__)


class Simulation(object):

    """Simulation base class."""

    def __init__(self, backend=None, output_path=None, steps=0,
                 checkpoint_interval=0, enable_speedometer=False,
                 restart=False):
        """
        Perform a simulation using the specified `backend` and optionally
        write output to `output_path`.

        Paths. To define output paths we rely on output_path, all
        other paths are defined based on it and on its
        base_path. Paths can then be defined locally by writers. Some
        glue is added in run_pre() to allow writers to cleanup their
        files.
        """
        self.backend = backend
        self.restart = restart
        self.output_path = output_path # can be None, file or directory
        self.max_steps = steps
        self.enable_speedometer = enable_speedometer
        self._checkpoint_scheduler = Scheduler(checkpoint_interval)
        self._targeter_steps = TargetSteps(self.max_steps)

        # Make sure the dirname of output_path exists. For instance,
        # if output_path is data/trajectory.xyz, then data/ should
        # exist. This creates the data/ folder and its parents folders.
        if self.output_path is not None:
            mkdir(os.path.dirname(self.output_path))
        if self.backend is None:
            self.backend = DryRunBackend()

        # Internal variables
        self._callback = []
        self.steps = 0
        self.initial_steps = 0
        self.start_time = time.time()
        # We expect subclasses to keep a ref to the trajectory object self.trajectory
        # used to store configurations, although this is not used in base class
        self.system = self.backend.system
        self.trajectory = self.backend.trajectory

        self.speedometer = None
        if enable_speedometer:
            self.speedometer = Speedometer()
            self.add(self.speedometer, Scheduler(None, calls=20, target=self.max_steps))

    def __str__(self):
        return 'atooms simulation via %s backend' % self.backend

    @property
    def base_path(self):
        # TODO: this if is not needed. If output_path is None, writers should never be called and simply disabled.
        if self.output_path is None:
            return None
        else:
            return os.path.splitext(self.output_path)[0]

    def add(self, callback, scheduler):
        """
        Register an observer `callback` to be called along with a
        `scheduler`.
        """
        # If the callback is already there we replace it
        # This allows to update targets / schedules on the way
        if callback in self._callback:
            self._callback.remove(callback)

        # Keep targeters last
        if not isinstance(callback, Target):
            callback.scheduler = scheduler
            self._callback.insert(0, callback)
        else:
            callback.scheduler = scheduler
            self._callback.append(callback)

    def remove(self, callback):
        """Remove the observer `callback`."""
        if callback in self._callback:
            self._callback.remove(callback)
        else:
            log.debug('attempt to remove inexistent callback %s (dont worry)', callback)

    def notify(self, observers):
        for o in observers:
            log.debug('notify %s at step %d', o, self.steps)
            o(self)

    @property
    def _targeters(self):
        return [o for o in self._callback if isinstance(o, Target)]

    @property
    def _non_targeters(self):
        return [o for o in self._callback if not isinstance(o, Target)]

    @property
    def _speedometers(self):
        return [o for o in self._callback if isinstance(o, Speedometer)]

    def write_checkpoint(self):
        # Tolerate missing implementation
        try:
            self.backend.write_checkpoint()
        except AttributeError:
            pass

    @property
    def rmsd(self):
        try:
            return self.backend.rmsd
        except AttributeError:
            return 0.0

    def elapsed_wall_time(self):
        return time.time() - self.start_time

    def wall_time_per_step(self):
        """
        Wall time per step in seconds.

        It can be subclassed by more complex simulation classes.
        """
        return self.elapsed_wall_time() / (self.steps-self.initial_steps)

    def wall_time_per_step_particle(self):
        """Wall time per step and particle in seconds."""
        try:
            # Be tolerant if there is no reference to system
            return self.wall_time_per_step() / len(self.system.particle)
        except AttributeError:
            return 0.0

    # Our template consists of two steps: run_pre() and run_until()
    # Typically a backend will implement the until method.
    # It is recommended to *extend* (not override) the base run_pre() in subclasses
    # TODO: when should checkpoint be read? The logic must be set here
    # Having a read_checkpoint() stub method would be OK here.
    def run_pre(self):
        """
        Preliminary step before run_until() to deal with restart
        conditions.
        """
        if self.output_path is not None:
            self.backend.output_path = self.output_path
            if not self.restart:
                # Clean up the trajectory folder and files.
                # Callbacks may implement their clean() methods
                for cbk in self._callback:
                    try:
                        cbk.clear(self)
                    except AttributeError:
                        pass

        self.backend.run_pre(self.restart)
        # If backend has reset the step because of restart, we update
        # it. Note that subclasses may overwrite this, because of
        # their own restart handling.
        if self.restart:
            self.steps = self.backend.steps
        barrier()

    def run_until(self, steps):
        """Run the simulation up to `steps`.

        Subclasses must set steps.
        """
        self.backend.run_until(steps)
        self.backend.steps = steps
        self.steps = steps

    def run(self, steps=None):
        """Run the simulation."""
        # If we are restaring we do not allow changing target steps on the fly.
        # because it might have side effects like non constant writing interval.
        if not self.restart or self.steps == 0:
            if steps is not None:
                self.max_steps = steps
            self.steps = 0
            self.backend.steps = 0

        # Targeter for max steps. Note that this will the replace an existing one.
        self._targeter_steps.target = self.max_steps
        self.add(self._targeter_steps, Scheduler(self.max_steps))

        self.report_header()
        self.run_pre()
        self.initial_steps = self.steps
        self.report()
        # Reinitialize speedometers
        for s in self._speedometers:
            s._init = False

        try:
            # Before entering the simulation, check if we can quit right away
            self.notify(self._targeters)
            # Then notify non targeters unless we are restarting
            if self.steps == 0:
                self.notify(self._non_targeters)
            else:
                self.notify(self._speedometers)
            log.info('starting at step: %d', self.steps)
            log.info('')
            while True:
                # Run simulation until any of the observers need to be called
                all_steps = [c.scheduler.next(self.steps) for c in self._callback]
                next_checkpoint = self._checkpoint_scheduler.next(self.steps)
                next_step = min(all_steps + [next_checkpoint])
                self.run_until(next_step)

                # Find observers indexes corresponding to minimum step
                # then get all corresponding observers
                next_step_ids = [i for i, step in enumerate(all_steps) if step == next_step]
                next_observers = [self._callback[i] for i in next_step_ids]

                # Observers should be sorted such that targeters are
                # last to avoid cropping output files
                self.notify(next_observers)
                if self.steps == next_checkpoint:
                    self.write_checkpoint()

        except SimulationEnd:
            # Checkpoint configuration at last step
            self.write_checkpoint()
            # We ignore errors due to performed steps being zero
            try:
                self._report_end()
            except:
                pass

        except KeyboardInterrupt:
            pass

        except:
            log.error('simulation failed')
            raise

        finally:
            log.info('goodbye')

    def report_header(self):
        txt = '%s' % self
        log.info('')
        log.info(txt)
        log.info('')
        log.info('atooms version: %s+%s (%s)', __version__, __commit__, __date__.split()[0])
        try:
            log.info('backend version: %s', self.backend.version)
        except:
            pass
        log.info('simulation starts on: %s', datetime.datetime.now().strftime('%Y-%m-%d at %H:%M'))
        log.info('output path: %s', self.output_path)

    def report(self):
        self._report()
        self._report_observers()

    def _report(self):
        """Implemented by subclasses"""
        pass

    def _report_observers(self):
        for f in self._callback:
            s = f.scheduler
            if isinstance(f, Target):
                log.info('target %s: %s', f, f.target)
            else:
                log.info('writer %s: interval=%s calls=%s', f, s.interval, s.calls)

    def _report_end(self):
        log.info('simulation ended on: %s', datetime.datetime.now().strftime('%Y-%m-%d at %H:%M'))
        log.info('final steps: %d', self.steps)
        log.info('final rmsd: %.2f', self.rmsd)
        log.info('wall time [s]: %.1f', self.elapsed_wall_time())
        log.info('average TSP [s/step/particle]: %.2e', self.wall_time_per_step_particle())
