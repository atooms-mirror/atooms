# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Simulation base class with callback logic"""

import os
import time
import datetime
import logging

from atooms.core import __version__, __date__
from atooms.utils import mkdir, barrier
from atooms.backends.dryrun import DryRunBackend

from .observers import *

log = logging.getLogger(__name__)

class Simulation(object):

    """Simulation base class."""

    # Trick to allow subclass to use custom observers for target and
    # writer without overriding setup(): have class variables to point
    # to the default observer classes that may be switched to custom
    # ones in subclasses
    _TARGET_STEPS = TargetSteps
    _TARGET_RMSD = TargetRMSD
    _WRITER_THERMO = WriterThermo
    _WRITER_CONFIG = WriterConfig
    _WRITER_CHECKPOINT = WriterCheckpoint

    def __init__(self, backend=DryRunBackend(),
                 output_path=None,
                 steps=0, rmsd=None,
                 target_steps=0, target_rmsd=None,
                 thermo_interval=0, thermo_number=0,
                 config_interval=0, config_number=0,
                 checkpoint_interval=0, checkpoint_number=0,
                 enable_speedometer=True,
                 restart=False):
        """
        Perform a simulation using the specified *backend* and writing
        output to *output_path*.

        Paths. To define output paths we rely on output_path, all
        other paths are defined based on it and on its
        base_path. Paths can then be defined locally by writers. Some
        glue is added in run_pre() to allow writers to cleanup their
        files.
        """
        self.backend = backend
        self.restart = restart
        self.output_path = output_path # can be None, file or directory
        self.target_steps = target_steps
        self.target_rmsd = target_rmsd
        self.thermo_interval = thermo_interval
        self.thermo_number = thermo_number
        self.config_interval = config_interval
        self.config_number = config_number
        self.checkpoint_interval = checkpoint_interval
        self.checkpoint_number = checkpoint_number
        self.enable_speedometer = enable_speedometer
        # Convenience shortcuts (might be dropped in the future)
        if steps > 0:
            self.target_steps = steps
        if rmsd is not None:
            self.target_rmsd = rmsd

        # Internal variables
        self._callback = []
        self.steps = 0
        self.initial_steps = 0
        self.start_time = time.time()
        # We expect subclasses to keep a ref to the trajectory object self.trajectory
        # used to store configurations, although this is not used in base class
        self.system = self.backend.system
        self.trajectory = self.backend.trajectory

        # Make sure the dirname of output_path exists. For instance,
        # if output_path is data/trajectory.xyz, then data/ should
        # exist. This creates the data/ folder and its parents folders.
        if self.output_path is not None:
            mkdir(os.path.dirname(self.output_path))

        # Setup writer callbacks
        self.targeter_rmsd_period = 10000
        self.targeter_steps = self._TARGET_STEPS(self.target_steps)
        self.targeter_rmsd = self._TARGET_RMSD(self.target_rmsd)
        self.writer_thermo = self._WRITER_THERMO()
        self.writer_config = self._WRITER_CONFIG()
        self.writer_checkpoint = self._WRITER_CHECKPOINT()
        self.speedometer = Speedometer()

    def __str__(self):
        return 'ATOOMS simulation (backend: %s)' % self.backend
    
    @property
    def base_path(self):
        # TODO: this if is not needed. If output_path is None, writers should never be called and simply disabled.
        if self.output_path is None:
            return None
        else:
            return os.path.splitext(self.output_path)[0]

    def _setup(self):
        """Add all internal observers to callbacks"""
        # This is called in run() if we are not restarting.
        # First make sure there are no default callbacks
        for c in [self.targeter_steps, self.targeter_rmsd,
                  self.writer_thermo, self.writer_config, self.writer_checkpoint,
                  self.speedometer]:
            self.remove(c)

        # Target steps are checked only at the end of the simulation
        self.targeter_steps = self._TARGET_STEPS(self.target_steps)
        self.add(self.targeter_steps, Scheduler(max(1, self.target_steps)))

        # TODO: implement dynamic scheduling or fail when interval is None and targeting is rmsd
        if self.target_rmsd is not None:
            self.targeter_rmsd = self._TARGET_RMSD(self.target_rmsd)
            self.add(self.targeter_rmsd, Scheduler(self.targeter_rmsd_period))

        # Setup schedulers
        if self.enable_speedometer:
            self.speedometer = Speedometer()
            self.add(self.speedometer, Scheduler(None, calls=20, target=self.target_steps))
        # TODO: if we are restarting and nsteps is increased, this will change the interval between cfgs
        if self.output_path is not None:
            self.add(self.writer_thermo, Scheduler(self.thermo_interval, self.thermo_number, self.target_steps))
            self.add(self.writer_config, Scheduler(self.config_interval, self.config_number, self.target_steps))
            self.add(self.writer_checkpoint, Scheduler(self.checkpoint_interval, self.checkpoint_number, self.target_steps))

    def add(self, callback, scheduler):
        """Register an observer (callback) to be called along with a scheduler"""
        # There are certainly more elegant ways of sorting Writers < Checkpoint < Target but anyway...
        # Keep targeters last
        if not isinstance(callback, Target):
            callback.scheduler = scheduler
            self._callback.insert(0, callback)
        else:
            callback.scheduler = scheduler
            self._callback.append(callback)

        # Enforce checkpoint is last among non_targeters
        try:
            idx = self._callback.index(self.writer_checkpoint)
            cnt = len(self._non_targeters)
            self._callback.insert(cnt-1, self._callback.pop(idx))
        except ValueError:
            pass

    def remove(self, callback):
        """Remove an observer (callback)"""
        if callback in self._callback:
            self._callback.remove(callback)
        else:
            log.debug('attempt to remove inexistent callback %s (dont worry)', callback)

    def notify(self, observers):
        for o in observers:
            log.debug('notify %s', o)
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

    @property
    def rmsd(self):
        try:
            return self.backend.rmsd
        except AttributeError:
            return 0.0

    def write_checkpoint(self):
        # Tolerate missing implementation
        # TODO: we should really check for write_checkpoint
        try:
            self.backend.write_checkpoint()
        except AttributeError:
            pass

    def elapsed_wall_time(self):
        return time.time() - self.start_time

    def wall_time_per_step(self):
        """
        Wall time per step in seconds.
        Can be conventiently subclassed by more complex simulation classes.
        """
        return self.elapsed_wall_time() / (self.steps-self.initial_steps)

    def wall_time_per_step_particle(self):
        """
        Wall time per step and particle in seconds.  Can be conventiently
        subclassed by more complex simulation classes.
        """
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
        """Deal with restart conditions before run we call_until()"""
        log.debug('calling backend pre at steps %d', self.steps)
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
        # /Design/: it is run_until responsability to set steps: self.steps = n
        # Bear it in mind when subclassing
        self.backend.run_until(steps)
        self.backend.steps = steps
        self.steps = steps

    def run(self, steps=None, rmsd=None):
        # If we are restaring we do not allow changing steps on the fly
        # Changing target_steps when restarting might have side effects like non constant writing interval.
        # TODO: spaghetti code
        if not self.restart or self.steps == 0:
            if steps is not None:
                self.target_steps = steps
                self.target_rmsd = rmsd
            if rmsd is not None:
                self.target_rmsd = rmsd
            self.steps = 0
            self.backend.steps = 0
            self._setup()
        self.run_pre()
        self.initial_steps = self.steps
        self.report()
        # Reinitialize speedometers
        for s in self._speedometers:
            s._init = False

        try:
            # Before entering the simulation, check if we can quit right away
            # Then notify non targeters unless we are restarting
            self.notify(self._targeters)
            if self.steps == 0:
                self.notify(self._non_targeters)
            else:
                self.notify(self._speedometers)
            log.info('starting at step: %d', self.steps)
            log.info('')
            while True:
                # Run simulation until any of the observers need to be called
                all_steps = [c.scheduler.next(self.steps) for c in self._callback]
                next_step = min(all_steps)
                # Find observers indexes corresponding to minimum step
                # then get all corresponding observers
                next_step_ids = [i for i, step in enumerate(all_steps) if step == next_step]
                next_obs = [self._callback[i] for i in next_step_ids]
                self.run_until(next_step)
                # Observers are sorted such that targeters are last
                # and checkpoint is last among writers
                self.notify(next_obs)

        except SimulationEnd as err:

            # Checkpoint configuration at last step
            self.writer_checkpoint(self)
            log.info(err.message)
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

    def report(self):
        txt = '%s' % self
        log.info('')
        log.info(txt)
        log.info('')
        log.info('atooms version: %s (%s)', __version__, __date__.split()[0])
        log.info('simulation starts on: %s', datetime.datetime.now().strftime('%Y-%m-%d at %H:%M'))
        log.info('output path: %s', self.output_path)
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


