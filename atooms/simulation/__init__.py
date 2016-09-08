# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Simulation base class with callback logic"""

import sys
import os
import time
import datetime
import logging
import copy

from atooms import __version__, __date__
from atooms.utils import mkdir
from atooms.utils import rank, size, barrier
from atooms.backends.dryrun import DryRunBackend

# Logging facilities

class ParallelFilter(logging.Filter):
    def filter(self, rec):
        if hasattr(rec, 'rank'):
            if rec.rank == 'all':
                return True
            else:
                return rank == rec.rank
        else:
            return rank == 0

class MyFormatter(logging.Formatter):
    def format(self, record):
        if record.levelname in ['WARNING', 'ERROR']:
            return '# %s: %s' % (record.levelname, record.msg)
        else:
            return '# %s' % record.msg

log = logging.getLogger('atooms')
formatter = MyFormatter()
handler = logging.StreamHandler(stream=sys.stdout)
handler.setFormatter(formatter)
log.addFilter(ParallelFilter())
log.addHandler(handler)
log.setLevel(logging.INFO)

# Default exceptions

class SimulationEnd(Exception):
    pass

class WallTimeLimit(Exception):
    pass

class SchedulerError(Exception):
    pass

class ParameterError(Exception):
    pass

# Default observers

# Different approaches are possible:
# 1. use callable classes passing args to __init__ and interval to add()
#    Exemple of client code:
#      sim.add(simulation.WriterThermo(), interval=100)
#      sim.add(simulation.TargetRMSD(5.0))

# 2. use functions passing interval and args as kwargs 
#    args are then passed to the function upon calling
#    In this case, we should differentiate types of callbacks via different add() methods
#    Exemple of client code:
#      sim.add(simulation.writer_thermo, interval=100)
#      sim.add(simulation.target_rmsd, rmsd=5.0)

# At least for default observers, it should be possible to use shortcuts
#   sim.writer_thermo_interval = 100
#   sim.target_rmsd = 5.0
# which will then add / modify the callbacks

# We identify callbacks by some types
# * target : these callbacks raise a SimualtionEnd when it's over
# * writer : these callbacks dump useful stuff to file
# and of course general purpose callback can be passed to do whatever

# It os the backend's responsibility to implement specific Writer observers.

# TODO: several related todos on where to store output file paths. It
# would make sense to keep them in Writers and delegate to them
# (instead of having them in simulation objects). This would require being able
# to access writers more directly.

class WriterCheckpoint(object):

    def __str__(self):
        return 'checkpoint'

    def __call__(self, e):
        e.write_checkpoint()

class WriterConfig(object):

    def __str__(self):
        return 'config'

    def __call__(self, e):
        if e.output_path is None:
            log.warning('config writing is request but output path is None')
            return
        f = os.path.join(e.output_path, 'trajectory.' + e.trajectory.suffix)
        with e.trajectory(f, 'a') as t:
            t.write(e.system, e.steps)

class WriterThermo(object):

    def __str__(self):
        return 'thermo'

    def __call__(self, e):
        if e.output_path is None:
            log.warning('thermo writing is request but output path is None')
            return
        f = os.path.join(e.output_path, 'trajectory.thermo')
        with open(f, 'a') as fh:
            fh.write('%d %g %g\n' % (e.steps, e.system.potential_energy(), e.rmsd))

def sec2time(t):
    eta_d = t / (24.0*3600)
    eta_h = (eta_d - int(eta_d)) * 24
    eta_m = (eta_h - int(eta_h)) * 60.0
    eta_s = (eta_m - int(eta_m)) * 60.0
    return '%dd:%02dh:%02dm:%02ds' % (eta_d, eta_h, eta_m, eta_s)


class Speedometer(object):

    def __init__(self):
        self._init = False

    def __str__(self):
        return 'speedometer'

    def __call__(self, e):
        if not self._init:
            # We could store all this in __init__() but this
            # way we allow targeters added to simulation via add()
            for c in e._callback:
                if isinstance(c, Target):
                    self.name_target = c.name
                    self.x_target = c.target
                    self.t_last = time.time()
                    # TODO: this assumes that targeters all get their target as attributes of simulation.
                    # We should fail or ask the targeter a cached value
                    self.x_last = float(getattr(e, self.name_target))
                    self._init = True
                    return

        if self.x_target > 0:
            t_now = time.time()
            x_now = float(getattr(e, self.name_target))
            # Get the speed at which the simulation advances
            speed = (x_now - self.x_last) / (t_now - self.t_last)
            # Report fraction of target achieved and ETA
            frac = float(x_now) / self.x_target
            try:
                eta = (self.x_target-x_now) / speed
                delta = sec2time(eta)
                d_now = datetime.datetime.now()
                d_delta = datetime.timedelta(seconds=eta)
                d_eta = d_now + d_delta
                log.info('%s: %d%% %s/%s estimated end: %s rate: %.2e TSP: %.2e' % \
                         (self.name_target, int(frac * 100), \
                          getattr(e, self.name_target), \
                          self.x_target,
                          d_eta.strftime('%Y-%m-%d %H:%M'),
                          speed, e.wall_time_per_step_particle()))
            except ZeroDivisionError:
                print x_now, self.x_last
                raise

        self.t_last = t_now
        self.x_last = x_now


class Target(object):

    def __init__(self, name, target):
        self.name = name
        self.target = target

    def __call__(self, sim):
        x = float(getattr(sim, self.name))        
        if self.target > 0:
            frac = float(x) / self.target
            log.debug('targeting %s now at %g [%d]' % (self.name, x, int(frac * 100)))
        if x >= self.target:
            raise SimulationEnd('achieved target %s: %s' % (self.name, self.target))

    def fraction(self, sim):
        """Fraction of target value already achieved"""
        return float(getattr(sim, self.name)) / self.target

    def __str__(self):
        return self.name
        
class TargetSteps(Target):

    # Note: this class is there just as an insane proof of principle
    # Steps targeting can/should be implemented by checking a simple int variable

    def __init__(self, target):
        Target.__init__(self, 'steps', target)
        
class TargetRMSD(Target):

    def __init__(self, target):
        Target.__init__(self, 'rmsd', target)

class TargetWallTime(Target):

    def __init__(self, wall_time):
        self.wtime_limit = wall_time

    def __call__(self, sim):
        if sim.elapsed_wall_time() > self.wtime_limit:
            raise WallTimeLimit('target wall time reached')
        else:
            t = sim.elapsed_wall_time()
            dt = self.wtime_limit - t
            log.debug('elapsed time %g, reamining time %g' % (t, dt))

class UserStop(object):
    """Allows a user to stop the simulation smoothly by touching a STOP
    file in the output root directory.
    Currently the file is not deleted to allow parallel jobs to all exit.
    """
    def __call__(self, e):
        # To make it work in parallel we should broadcast and then rm 
        # or subclass userstop in classes that use parallel execution
        if e.output_path is not None and self.STORAGE == 'directory':
            log.debug('User Stop %s/STOP' % e.output_path)
            # TODO: support files as well
            if os.path.exists('%s/STOP' % e.output_path):
                raise SimulationEnd('user has stopped the simulation')
        else:
            raise RuntimeError('USerStop wont work atm with file storage')

class Scheduler(object):

    #TODO: interval can be a function to allow non linear sampling
    #TODO: base scheduler plus derived scheduler for fixed ncalls 

    """Scheduler to call observer during the simulation"""

    def __init__(self, interval, calls=None, target=None):
        self.interval = interval
        self.calls = calls
        self.target = target
        # Main switch
        if interval>0:
            # Fixed interval.
            self.interval = interval
        else:
            if calls>0:
                # Fixed number of calls.
                if self.target is not None:
                    # If both calls and target are not None, we determine interval
                    self.interval = max(1, self.target / self.calls)
                else:
                    # Dynamic scheduling
                    raise SchedulerError('dynamic scheduling not implemented')
            # else:
            #     # Scheduler is disabled
            #     self.interval = 0 #sys.maxint

    def next(self, this):
        if self.interval>0:
            return (this / self.interval + 1) * self.interval
        else:
            return sys.maxint

class OnetimeScheduler(object):

    """Scheduler to call observer during the simulation"""

    def __init__(self, interval, sim):
        self.sim = sim
        self.interval = interval
        self.calls = None
        self.target = None

    def next(self, this):
        log.debug('one time initial steps %d this %d' % (self.sim.initial_steps, this))
        if (this-self.sim.initial_steps) / self.interval == 0:
            return self.interval
        else:
            return sys.maxint


class Simulation(object):

    """Simulation base class."""

    # TODO: write initial configuration as well

    # Comvoluted trick to allow subclass to use custom observers for
    # target and writer without overriding setup(): have class
    # variables to point to the default observer classes that may be
    # switched to custom ones in subclasses
    _TARGET_STEPS = TargetSteps
    _TARGET_RMSD = TargetRMSD
    _WRITER_THERMO = WriterThermo
    _WRITER_CONFIG = WriterConfig    
    _WRITER_CHECKPOINT = WriterCheckpoint

    # Storage class attribute: backends can either dump to a directory
    # or following a file suffix logic. Meant to be subclassed.
    STORAGE = 'directory'

    def __init__(self, backend=DryRunBackend(),
                 output_path=None, 
                 steps=0, rmsd=None,
                 target_steps=0, target_rmsd=None,
                 thermo_interval=0, thermo_number=0, 
                 config_interval=0, config_number=0,
                 checkpoint_interval=0, checkpoint_number=0,
                 enable_speedometer=True,
                 restart=False):
        """We expect input and output paths as input.
        Alternatively, input might be a system (or trajectory?) instance.
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
        if steps>0:
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

        # Setup writer callbacks
        # TODO: if output_path is None we should disable writers
        self.targeter_rmsd_period = 10000
        self.targeter_steps = self._TARGET_STEPS(self.target_steps)
        self.targeter_rmsd = self._TARGET_RMSD(self.target_rmsd)
        self.writer_thermo = self._WRITER_THERMO()
        self.writer_config = self._WRITER_CONFIG()
        self.writer_checkpoint = self._WRITER_CHECKPOINT()
        self.speedometer = Speedometer()

    def __str__(self):
        return 'ATOOMS simulation (backend: %s)' % self.backend
        
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
            log.debug('attempt to remove inexistent callback %s (dont worry)' % callback)

    @property
    def _targeters(self):
        return [o for o in self._callback if isinstance(o, Target)]

    @property
    def _non_targeters(self):
        return [o for o in self._callback if not isinstance(o, Target)]

    @property
    def _speedometers(self):
        return [o for o in self._callback if isinstance(o, Speedometer)]

    def notify(self, observers):
        for o in observers:
            log.debug('notify %s' % o)
            o(self)
    
    @property
    def rmsd(self):
        # TODO: provide rmsd by species 07.12.2014
        try:
            return self.backend.rmsd
        except AttributeError:
            log.warning('rmsd has not been subclassed')
            return 0.0

    def write_checkpoint(self):
        self.backend.write_checkpoint()

    def elapsed_wall_time(self):
        return time.time() - self.start_time

    def wall_time_per_step(self):
        """Return the wall time in seconds per step.
        Can be conventiently subclassed by more complex simulation classes.
        """
        return self.elapsed_wall_time() / (self.steps-self.initial_steps)

    def wall_time_per_step_particle(self):
        """Return the wall time in seconds per step and particle.
        Can be conventiently subclassed by more complex simulation classes.
        """
        try:
            # Be tolerant if there is no reference to system
            return self.wall_time_per_step() / len(self.system.particle)
        except:
            return 0.0

    # Our template consists of two steps: run_pre() and run_until()
    # Typically a backend will implement the until method.
    # It is recommended to *extend* (not override) the base run_pre() in subclasses
    # TODO: when should checkpoint be read? The logic must be set here
    # Having a read_checkpoint() stub method would be OK here.
    def run_pre(self):
        """Deal with restart conditions before run we call_until()"""
        # TODO: moved from init to here, is it ok?
        log.debug('calling backend pre at steps %d' % self.steps)
        if self.output_path is not None and self.STORAGE == 'directory':
            mkdir(self.output_path)
        if self.output_path is not None:
            self.backend.output_path = self.output_path
        self.backend.run_pre(self.restart)
        barrier()
            
    def run_until(self, n):
        # /Design/: it is run_until responsability to set steps: self.steps = n
        # Bear it in mind when subclassing 
        self.backend.run_until(n)
        self.backend.steps = n
        self.steps = n

    def run(self, steps=None, rmsd=None):
        # If we are restaring we do not allow changing steps on the fly
        if not self.restart or self.steps==0:
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
            log.info('starting at step: %d' % self.steps)
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

        except SimulationEnd as s:
            # Checkpoint configuration at last step
            self.writer_checkpoint(self)
            log.info('%s' % s.message)
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
        nch = len(txt)
        log.info('')
        log.info(txt)
        log.info('')
        log.info('atooms version: %s (%s)' % (__version__, __date__.split()[0]))
        log.info('simulation starts on: %s' % datetime.datetime.now().strftime('%Y-%m-%d at %H:%M'))
        self._report()
        self._report_observers()
        
    def _report(self):
        """Implemented by subclassed"""
        pass

    def _report_observers(self):
        for f in self._callback:
            s = f.scheduler
            if isinstance(f, Target):
                log.info('target %s: %s' % (f, f.target)) #, s.interval, s.calls))
            else:
                log.info('writer %s: interval=%s calls=%s' % (f, s.interval, s.calls))

    def _report_end(self):
        log.info('simulation ended on: %s' % datetime.datetime.now().strftime('%Y-%m-%d at %H:%M'))
        log.info('final steps: %d' % self.steps)
        log.info('final rmsd: %.2f' % self.rmsd)
        log.info('wall time [s]: %.1f' % self.elapsed_wall_time())
        log.info('average TSP [s/step/particle]: %.2e' % (self.wall_time_per_step_particle()))


