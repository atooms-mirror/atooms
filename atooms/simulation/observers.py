# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""
Scheduler and callbacks (aka observers) to be called during a simulation.

To add a callback `func` to a `Simulation` instance `sim` and have it
called every 100 steps

    #!python
    sim.add(func, Scheduler(100))

There are two ways to set them up and add them to a simulation:

1. callable classes

   Example of client code:
     sim.add(simulation.WriterThermo(), 100)
     sim.add(simulation.TargetRMSD(5.0))

2. use functions passing optional *args and **kwargs.

   Example of client code:
     sim.add(simulation.writer_thermo, 100)
     sim.add(simulation.target_rmsd, 100, rmsd=5.0)

To differentiate different types of callbacks, we following a naming
convention (for classes and function). If they contain

- target : these callbacks raise a SimulationEnd when it's over
- writer : these callbacks dump useful stuff to file

Of course, general purpose callback can be passed to do whatever.
"""

import sys
import os
import shutil
import time
import datetime
import logging

__all__ = ['SimulationEnd', 'WallTimeLimit', 'Scheduler',
           'write_config', 'write_thermo', 'write', 'target',
           'target_rmsd', 'target_steps', 'target_walltime',
           'user_stop', 'Speedometer']

_log = logging.getLogger(__name__)


# Helper functions

def _sec2time(time_interval):
    """
    Convert a time interval in seconds to (day, hours, minutes,
    seconds) format.
    """
    eta_d = time_interval / (24.0 * 3600)
    eta_h = (eta_d - int(eta_d)) * 24
    eta_m = (eta_h - int(eta_h)) * 60.0
    eta_s = (eta_m - int(eta_m)) * 60.0
    return '%dd:%02dh:%02dm:%02ds' % (eta_d, eta_h, eta_m, eta_s)


# Default exceptions

class SimulationEnd(Exception):
    """Raised when an targeter reaches its target."""
    pass


class WallTimeLimit(Exception):
    """Raised when the wall time limit is reached."""
    pass


# Scheduler classes

class Scheduler(object):

    """
    Schedule observer calls during the simulation.

    This is nothing but a callable that takes a simulation instance
    and returns the next step at which an observer has to to notified.
    """

    def __init__(self, interval=None, calls=None, steps=None,
                 block=None, seconds=None):
        """
        Only one of the arguments can be different from None.

        - `interval`: notify at a fixed steps interval
        - `calls`: fixed number of notification
        - `steps`: list of steps at which the observer will be notified
        - `block`: as steps, but will be called periodically
        - `seconds`: notify every `seconds`
        """
        self.interval = interval
        self.calls = calls
        self.steps = steps
        self.block = block
        self.seconds = seconds

    def __call__(self, sim):
        """
        Given a simulation instance `sim`, return the next step at which
        the observer will be called.
        """
        if self.interval is not None and self.interval > 0:
            # Regular interval
            return (sim.current_step / self.interval + 1) * self.interval
        elif self.calls is not None and self.interval > 0:
            # Fixed number of calls
            interval = int(sim.steps / self.calls)
            return (sim.current_step / interval + 1) * interval
        elif self.steps is not None:
            # List of selected steps
            inext = self.steps[0]
            for i, step in enumerate(self.steps[:-1]):
                if sim.steps >= step:
                    inext = self.steps[i+1]
                    break
            return inext
        elif self.block is not None:
            # like steps but with % on sim.steps
            pass
        elif self.seconds is not None:
            pass
        else:
            return sys.maxint


# Writer callbacks
# Callbacks as pure function to distinguish their role we adopt a naming convention:
# if the callback contains write (target) in its __name__ then it is a writer (targeter).

def write_config(sim):
    """
    Write configurations to a trajectory file.

    The trajectory format is taken from the passed Simulation
    instance.
    """
    if sim.current_step == 0:
        # TODO: folder-based trajectories should ensure that mode='w' clears up the folder
        # TODO: refactor as rm()
        if os.path.isdir(sim.output_path):
            shutil.rmtree(sim.output_path)
        elif os.path.isfile(sim.output_path):
            os.remove(sim.output_path)

    with sim.trajectory(sim.output_path, 'a') as t:
        t.write(sim.system, sim.current_step)

def write_thermo(sim):
    """Write basic thermodynamic data."""
    f = sim.output_path + '.thermo'
    if sim.current_step == 0:
        with open(f, 'w') as fh:
            fh.write('# columns:' + ', '.join(['steps',
                                               'temperature',
                                               'potential energy',
                                               'kinetic energy',
                                               'total energy',
                                               'rmsd']) + '\n')
    with open(f, 'a') as fh:
        fh.write('%d %g %g %g %g %g\n' % (sim.current_step,
                                          sim.system.temperature,
                                          sim.system.potential_energy(normed=True),
                                          sim.system.kinetic_energy(normed=True),
                                          sim.system.total_energy(normed=True),
                                          sim.rmsd))

def write(sim, name, attributes):
    """
    Write generic attributes of simulation and system to a file.

    `name` is a tag appended to `sim.base_path` to define the output
    file path.

    `attributes` must be a list of valid properties of the Simulation
    instance `sim` or of its System instance `sim.system`.
    """
    f = sim.output_path + '.' + name
    if sim.current_step == 0:
        with open(f, 'w') as fh:
            fh.write('# columns: %s\n' % ', '.join(attributes))
    else:
        # Extract the requested attributes
        values = []
        for attr in attributes:
            level = len(attr.split('.'))
            if level == 1:
                values.append(getattr(sim, attr))
            elif level == 2:
                system_attr = attr.split('.')[-1]
                if attr.startswith('system'):
                    values.append(getattr(sim.system, system_attr))
            else:
                raise ValueError('attribute is too deep')
        # Format output string
        fmt = ('%s ' * len(attributes)) + '\n'
        with open(f, 'a') as fh:
            fh.write(fmt % tuple(values))


# Target callbacks

def target(sim, attribute, value):
    """
    An observer that raises a `SimulationEnd` exception when a given
    target `value` of a property is reached during a simulation. The
    property is `attribute` and is assumed to be an attribute of
    simulation.
    """
    x = float(getattr(sim, attribute))
    if value > 0:
        frac = float(x) / value
        _log.debug('target %s now at %g [%d]', attribute, x, int(frac * 100))
    if x >= value:
        raise SimulationEnd('reached target %s: %s', attribute, value)
    return frac

def target_rmsd(sim, value):
    """Target the root mean squared displacement."""
    return target(sim, 'rmsd', value)

def target_steps(sim, value):
    """Target the number of steps."""
    if sim.current_step >= value:
        raise SimulationEnd('reached target steps %d' % value)

def target_walltime(sim, value):
    """
    Target a value of the elapsed wall time from the beginning of the
    simulation.

    Useful to self restarting jobs in a queining system with time
    limits.
    """
    wtime_limit = value
    if sim.elapsed_wall_time() > wtime_limit:
        raise WallTimeLimit('target wall time reached')
    else:
        t = sim.elapsed_wall_time()
        dt = wtime_limit - t
        _log.debug('elapsed time %g, reamining time %g', t, dt)

def user_stop(sim):
    """
    Allows a user to stop the simulation smoothly by touching a STOP
    file in the output root directory.  Currently the file is not
    deleted to allow parallel jobs to all exit.
    """
    # To make it work in parallel we should broadcast and then rm
    # or subclass userstop in classes that use parallel execution
    # TODO: support files as well
    if sim.output_path is not None:
        try:
            _log.debug('User Stop %s/STOP', sim.output_path)
            if os.path.exists('%s/STOP' % sim.output_path):
                raise SimulationEnd('user has stopped the simulation')
        except IOError:
            raise IOError('user_stop wont work atm with file storage')


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
                if 'target' in c.__name__.lower():
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
                eta = (self.x_target - x_now) / speed
                d_now = datetime.datetime.now()
                d_delta = datetime.timedelta(seconds=eta)
                d_eta = d_now + d_delta
                _log.info('%s: %d%% %s/%s estimated end: %s rate: %.2e TSP: %.2e',
                          self.name_target, int(frac * 100),
                          getattr(sim, self.name_target),
                          self.x_target,
                          d_eta.strftime('%Y-%m-%d %H:%M'),
                          speed, sim.wall_time_per_step_particle())
            except ZeroDivisionError:
                print x_now, self.x_last
                raise

        self.t_last = t_now
        self.x_last = x_now
