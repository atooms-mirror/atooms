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
from .observers import target_steps, Speedometer, Scheduler, SimulationEnd

log = logging.getLogger(__name__)


class Simulation(object):

    """Simulation base class."""

    def __init__(self, backend, output_path=None, steps=0,
                 checkpoint_interval=0, enable_speedometer=False,
                 restart=False):
        """
        Perform a simulation using the specified `backend` and optionally
        write output to `output_path`. This can be a file or directory path.

        Paths. To define output paths we rely on output_path, all
        other paths are defined based on it and on its
        base_path. Paths can then be defined locally by writers. Some
        glue is added in run_pre() to allow writers to cleanup their
        files.
        """
        self.backend = backend
        self.restart = restart
        self.output_path = output_path
        # TODO: use steps for the target steps, current_step instead of steps
        self.max_steps = steps
        self.steps = 0
        self.initial_step = 0
        # We expect subclasses to keep a ref to the trajectory object
        # self.trajectory used to store configurations
        self.trajectory = self.backend.trajectory

        # Make sure the dirname of output_path exists. For instance,
        # if output_path is data/trajectory.xyz, then data/ should
        # exist. This creates the data/ folder and its parents folders.
        if self.output_path is not None:
            mkdir(os.path.dirname(self.output_path))

        # Internal variables
        self._callback = []
        self._start_time = time.time()
        self._speedometer = None
        self._checkpoint_scheduler = Scheduler(checkpoint_interval)
        self._targeter_steps = target_steps
        self._cbk_params = {}  # hold scheduler and parameters of callbacks
        if enable_speedometer:
            self._speedometer = Speedometer()
            self.add(self._speedometer, Scheduler(None, calls=20,
                                                  target=self.max_steps))

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
        return 'atooms simulation via %s backend' % self.backend

    @property
    def base_path(self):
        return os.path.splitext(self.output_path)[0]

    def add(self, callback, scheduler, *args, **kwargs):
        """
        ASS an observer `callback` to be called along with a `scheduler`.

        `scheduler` and `callback` must be callables accepting a
        Simulation instance as unique argument. `scheduler` must
        return the next step at which the observer has to be notified.

        An integer value is allowed for `scheduler`. In this case, a
        scheduler with fixed interval is generated internally and the
        observer is notified every `scheduler` steps.
        """
        # If the callback is already there we replace it
        # This allows to update targets / schedules on the way
        # TODO: this way we cannot the same observer with different schedules
        if callback in self._callback:
            self._callback.remove(callback)

        # Accept an integer interval
        if type(scheduler) is int:
            scheduler = Scheduler(scheduler)

        # Store scheduler, callback and its arguments
        # in a separate dict (NOT in the function object itself!)
        self._cbk_params[callback] = {'scheduler': scheduler, 
                                      'args': args, 
                                      'kwargs': kwargs}

        # Keep targeters last
        if 'target' not in callback.__name__.lower():
            self._callback.insert(0, callback)
        else:
            self._callback.append(callback)

    def remove(self, callback):
        """Remove the observer `callback`."""
        if callback in self._callback:
            self._callback.remove(callback)
            self._cbk_params.pop(callback)
        else:
            log.debug('attempt to remove inexistent callback %s (dont worry)', callback)

    def _notify(self, observers):
        for observer in observers:
            log.debug('notify %s at step %d', observer, self.steps)
            args = self._cbk_params[observer]['args']
            kwargs = self._cbk_params[observer]['kwargs']
            observer(self, *args, **kwargs)

    @property
    def _targeters(self):
        return [o for o in self._callback if 'target' in o.__name__.lower()]

    @property
    def _non_targeters(self):
        return [o for o in self._callback if 'target' not in o.__name__.lower()]

    @property
    def _speedometers(self):
        return [o for o in self._callback if isinstance(o, Speedometer)]

    def write_checkpoint(self):
        if self.output_path is not None:
            with open(self.output_path + '.chk.step', 'w') as fh:
                fh.write('%d' % self.steps)
        try:
            self.backend.write_checkpoint()
        except AttributeError:
            # Tolerate missing checkpoint implementation
            pass

    def read_checkpoint(self):
        if self.output_path is not None:
            with open(self.output_path + '.chk.step') as fh:
                self.steps = int(fh.read())
        try:
            self.backend.read_checkpoint()
        except AttributeError:
            # Tolerate missing checkpoint implementation
            pass

    @property
    def rmsd(self):
        try:
            return self.backend.rmsd
        except AttributeError:
            return 0.0

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
        if len(self.system.particle) > 0:
            norm *= len(self.system.particle)
        else:
            return float('nan')
        if self.steps - self.initial_step > 0:
            norm *= (self.steps - self.initial_step)
        else:
            return float('nan')
        return self._elapsed_wall_time() / norm

    # Our template consists of two steps: run_pre() and run_until()
    # Typically a backend will implement the until method.
    # It is recommended to *extend* (not override) the base run_pre() in subclasses
    # TODO: when should checkpoint be read? The logic must be set here
    # Having a read_checkpoint() stub method would be OK here.
    def initialize(self):
        """
        Preliminary step before run_until() to deal with restart
        conditions.
        """
        self._start_time = time.time()
        # TODO: This way the backend inherits the output path and no need to set it there. We could do it the other way round?
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

            if self.restart:
                self.read_checkpoint()

        barrier()

    def run_until(self, steps):
        """
        Run the simulation up to `steps`.

        Subclasses must set steps.
        """
        self.backend.run(steps - self.steps)
        #self.backend.steps = steps
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

        # Targeter for max steps. This will the replace an existing one.
        self.add(self._targeter_steps, Scheduler(self.max_steps), self.max_steps)
        self._report_header()
        self.initialize()
        self.initial_step = self.steps
        self._report()
        self._report_observers()
        # Reinitialize speedometers
        for s in self._speedometers:
            s._init = False

        try:
            # Before entering the simulation, check if we can quit right away
            self._notify(self._targeters)
            # Then notify non targeters unless we are restarting
            if self.steps == 0:
                self._notify(self._non_targeters)
            else:
                self._notify(self._speedometers)
            log.info('starting at step: %d', self.steps)
            log.info('')
            while True:
                # Run simulation until any of the observers need to be called
                all_steps = [self._cbk_params[c]['scheduler'](self) for c in self._callback]
                next_checkpoint = self._checkpoint_scheduler(self)
                next_step = min(all_steps + [next_checkpoint])
                self.run_until(next_step)

                # Find observers indexes corresponding to minimum step
                # then get all corresponding observers
                next_step_ids = [i for i, step in enumerate(all_steps) if step == next_step]
                next_observers = [self._callback[i] for i in next_step_ids]

                # Observers should be sorted such that targeters are
                # last to avoid cropping output files
                self._notify(next_observers)
                if self.steps == next_checkpoint:
                    self.write_checkpoint()
        except SimulationEnd:
            # Checkpoint configuration at last step
            self.write_checkpoint()
            self._report_end()
        except KeyboardInterrupt:
            pass
        except:
            log.error('simulation failed')
            raise
        finally:
            log.info('goodbye')

    def _report_header(self):
        txt = '%s' % self
        log.info('')
        log.info(txt)
        log.info('')
        log.info('atooms version: %s+%s (%s)', __version__, __commit__, __date__)
        try:
            log.info('backend version: %s', self.backend.version)
        except AttributeError:
            pass
        log.info('simulation starts on: %s', datetime.datetime.now().strftime('%Y-%m-%d at %H:%M'))
        log.info('output path: %s', self.output_path)

    def _report(self):
        """Implemented by subclasses"""
        pass

    def _report_observers(self):
        for f in self._callback:
            params = self._cbk_params[f]
            s = params['scheduler']
            if 'target' in f.__name__.lower():
                args = params['args']
                log.info('target %s: %s', f.__name__, args[0])
            else:
                log.info('writer %s: interval=%s calls=%s', f.__name__, s.interval, s.calls)

    def _report_end(self):
        now = datetime.datetime.now().strftime('%Y-%m-%d at %H:%M')
        log.info('simulation ended on: %s', now)
        log.info('final steps: %d', self.steps)
        log.info('final rmsd: %.2f', self.rmsd)
        log.info('wall time [s]: %.1f', self.wall_time())
        log.info('average TSP [s/step/particle]: %.2e', 
                 self.wall_time(per_step=True, per_particle=True))
