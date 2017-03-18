# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Default observers and schedulers for simulations.

To add a callback `func` to a `Simulation` instance `sim` and have it called every 100 steps

    #!python
    sim.add(func, Scheduler(100))
"""

import sys
import os
import shutil
import time
import datetime
import logging

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

log = logging.getLogger(__name__)

# Default exceptions

class SimulationEnd(Exception):
    pass

class WallTimeLimit(Exception):
    pass

class SchedulerError(Exception):
    pass


# Writers

class WriterConfig(object):

    """
    Callable class that writes configurations to a trajectory file.

    The trajectory format is taken from the passed Simulation instance
    that calls the callbacks.
    """

    def __str__(self):
        return 'config'

    def __call__(self, sim):
        with sim.trajectory(sim.output_path, 'a') as t:
            t.write(sim.system, sim.steps)

    def clear(self, sim):
        # TODO: refactor as rm()
        if os.path.isdir(sim.output_path):
            shutil.rmtree(sim.output_path)
        elif os.path.isfile(sim.output_path):
            os.remove(sim.output_path)

class WriterThermo(object):

    """Callable class that writes thermodynamic data to disk."""

    def __str__(self):
        return 'thermo'

    def __call__(self, sim):
        f = sim.base_path + '.thermo'
        with open(f, 'a') as fh:
            fh.write('%d %g %g\n' % (sim.steps, sim.system.potential_energy(), sim.rmsd))

    def clear(self, sim):
        f = sim.base_path + '.thermo'
        if os.path.exists(f):
            os.remove(f)

def sec2time(time_interval):
    """
    Convert a time interval in seconds to (day, hours, minutes,
    seconds) format.
    """
    eta_d = time_interval / (24.0*3600)
    eta_h = (eta_d - int(eta_d)) * 24
    eta_m = (eta_h - int(eta_h)) * 60.0
    eta_s = (eta_m - int(eta_m)) * 60.0
    retur
    n '%dd:%02dh:%02dm:%02ds' % (eta_d, eta_h, eta_m, eta_s)


class Speedometer(object):

    """Display speed of simulation and remaining time to reach target."""

    def __init__(self):
        self._init = False

    def __str__(self):
        return 'speedometer'

    def __call__(self, sim):
        if not self._init:
            # We could store all this in __init__() but this
            # way we allow targeters added to simulation via add()
            for c in sim._callback:
                if isinstance(c, Target):
                    self.name_target = c.name
                    self.x_target = c.target
                    self.t_last = time.time()
                    # TODO: this assumes that targeters all get their target as attributes of simulation.
                    # We should fail or ask the targeter a cached value
                    self.x_last = float(getattr(sim, self.name_target))
                    self._init = True
                    return

        if self.x_target > 0:
            t_now = time.time()
            x_now = float(getattr(sim, self.name_target))
            # Get the speed at which the simulation advances
            speed = (x_now - self.x_last) / (t_now - self.t_last)
            # Report fraction of target achieved and ETA
            frac = float(x_now) / self.x_target
            try:
                eta = (self.x_target-x_now) / speed
                d_now = datetime.datetime.now()
                d_delta = datetime.timedelta(seconds=eta)
                d_eta = d_now + d_delta
                log.info('%s: %d%% %s/%s estimated end: %s rate: %.2e TSP: %.2e', \
                         self.name_target, int(frac * 100), \
                         getattr(sim, self.name_target), \
                         self.x_target,
                         d_eta.strftime('%Y-%m-%d %H:%M'),
                         speed, sim.wall_time_per_step_particle())
            except ZeroDivisionError:
                print x_now, self.x_last
                raise

        self.t_last = t_now
        self.x_last = x_now


class Target(object):

    """Base targeter class."""

    def __init__(self, name, target):
        self.name = name
        self.target = target

    def __call__(self, sim):
        x = float(getattr(sim, self.name))
        if self.target > 0:
            frac = float(x) / self.target
            log.debug('targeting %s now at %g [%d]', self.name, x, int(frac * 100))
        if x >= self.target:
            raise SimulationEnd('achieved target %s: %s', self.name, self.target)

    def fraction(self, sim):
        """Fraction of target value already achieved"""
        return float(getattr(sim, self.name)) / self.target

    def __str__(self):
        return self.name

class TargetSteps(Target):

    """
    Targeter a fixed number of steps.

    Note: this class is here as an insane proof of principle. Steps
    targeting can (should?) be implemented in `Simulation` by checking
    a simple integer variable.
    """

    def __init__(self, target):
        Target.__init__(self, 'steps', target)

class TargetRMSD(Target):

    """Target a value of the total root mean squared displacement."""

    def __init__(self, target):
        Target.__init__(self, 'rmsd', target)

class TargetWallTime(Target):

    """
    Target a value of the elapsed wall time from the beginning of the
    simulation.

    Useful to self restarting jobs in a queining system with time
    limits.
    """

    def __init__(self, wall_time):
        self.wtime_limit = wall_time

    def __call__(self, sim):
        if sim.elapsed_wall_time() > self.wtime_limit:
            raise WallTimeLimit('target wall time reached')
        else:
            t = sim.elapsed_wall_time()
            dt = self.wtime_limit - t
            log.debug('elapsed time %g, reamining time %g', t, dt)

class UserStop(object):
    """Allows a user to stop the simulation smoothly by touching a STOP
    file in the output root directory.
    Currently the file is not deleted to allow parallel jobs to all exit.
    """
    def __call__(self, sim):
        # To make it work in parallel we should broadcast and then rm
        # or subclass userstop in classes that use parallel execution
        if sim.output_path is not None and sim.storage == 'directory':
            log.debug('User Stop %s/STOP', sim.output_path)
            # TODO: support files as well
            if os.path.exists('%s/STOP' % sim.output_path):
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

        if interval > 0:
            # Fixed interval.
            self.interval = interval
        else:
            if calls > 0:
                # Fixed number of calls.
                if self.target is not None:
                    # If both calls and target are not None, we determine interval
                    self.interval = max(1, self.target / self.calls)
                else:
                    # Dynamic scheduling
                    raise SchedulerError('dynamic scheduling not implemented')

    def next(self, this):
        if self.interval > 0:
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
        log.debug('one time initial steps %d this %d', self.sim.initial_steps, this)
        if (this-self.sim.initial_steps) / self.interval == 0:
            return self.interval
        else:
            return sys.maxint
