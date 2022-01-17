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
import time
import datetime
import logging
from atooms.core.utils import rmd, rmf

__all__ = ['SimulationEnd', 'WallTimeLimit', 'SimulationKill',
           'Scheduler', 'write_config', 'write_thermo', 'store',
           'write_trajectory', 'write', 'target', 'target_rmsd',
           'target_steps', 'target_walltime', 'user_stop',
           'target_user_stop', 'Speedometer', 'shell_stop',
           'target_shell_stop', 'target_python_stop']

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

class SimulationKill(Exception):
    """Raised when a simulation is terminated by SIGTERM."""
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

        # Normalize non-positive intervals and n. of calls
        if self.interval is not None and self.interval <= 0:
            self.interval = None
        if self.calls is not None and self.calls <= 0:
            self.calls = None

    def __call__(self, sim):
        """
        Given a simulation instance `sim`, return the next step at which
        the observer will be called.
        """
        if self.interval is not None and self.calls is None:
            # Regular interval
            return (sim.current_step // self.interval + 1) * self.interval

        elif self.calls is not None:
            # Fixed number of calls
            interval = max(1, sim.steps // self.calls)
            return (sim.current_step // interval + 1) * interval

        elif self.steps is not None:
            # List of selected steps
            inext = sys.maxsize
            for i, step in enumerate(self.steps):
                if step > sim.current_step:
                    inext = self.steps[i]
                    break
            return inext

        elif self.block is not None:
            # Periodic block of steps
            step_of_last_block = (sim.current_step // self.block[-1]) * self.block[-1]
            inext = sys.maxsize
            for i, step in enumerate(self.block):
                if step > sim.current_step % self.block[-1]:
                    inext = self.block[i] + step_of_last_block
                    break
            return inext

        elif self.seconds is not None:
            pass
        else:
            return sys.maxsize


# Writer callbacks
# Callbacks as pure function to distinguish their role we adopt a naming convention:
# if the callback contains write (target) in its __name__ then it is a writer (targeter).

def write_to_ram(sim, trajectory_ram):
    """
    Write configurations to a trajectory in ram.
    """
    # TODO: deprecate this since write_trajectory now covers it
    trajectory_ram.write(sim.system, sim.current_step)


def write_trajectory(sim, variables=None, precision=None, trajectory=None,
                     trajectory_class=None, fields=None):
    """
    Write trajectory frame from `sim` Simulation instance

    The trajectory format is taken from `sim.trajectory_class` and a
    local instance with that format is used appending the frames on
    successive calls to this function.

    If `trajectory` is a Trajectory instance, it is used instead.
    """
    # TODO: deprecate fields in favor of variables
    if fields is not None:
        variables = fields

    # Clear up everything
    if sim.current_step == 0 and trajectory is None and \
       trajectory_class is None:
        # TODO: bug here with trajectory the file is removed!
        # TODO: folder-based trajectories should ensure that mode='w' clears up the folder
        rmd(sim.output_path)
        rmf(sim.output_path)

    if trajectory is None and trajectory_class is None:
        th = sim.trajectory_class(sim.output_path, 'a')
    elif trajectory_class is not None:
        th = trajectory_class(sim.output_path, 'a')
    else:
        th = trajectory

    # Write stuff
    if hasattr(sim.backend, 'timestep'):
        th.timestep = sim.backend.timestep
    if precision is not None:
        th.precision = precision
    if variables is not None:
        th.variables = variables
    th.write(sim.system, sim.current_step)

    if trajectory is None:
        th.close()


# Deprecated alias
write_config = write_trajectory

def write_thermo(sim, fields=None, fmt=None, precision=6, functions=None):
    """
    Write thermodynamic properties to a file.

    By default, available fields are:
    - steps
    - temperature
    - potential energy per particle
    - kinetic energy per particle
    - total energy
    - pressure
    - rmsd

    The set of available `fields` can be augmented by passing an extra
    `functions` dictionary.

    The `fmt` dictionary can be used to provide custom formatting
    options for fields.

    The `precision` parameter controls the default precision of
    floating point fields.
    """
    # TODO: deprecate fields
    # By default write minimal info
    if fields is None:
        fields = ['steps',
                  'temperature',
                  'potential energy per particle',
                  'kinetic energy per particle',
                  'total energy per particle',
                  'rmsd']

    # Internal function database.
    # It can be augmented via functions parameter
    _db_func = {
        'steps': lambda x: x.current_step,
        'potential energy per particle': lambda x: x.system.potential_energy(True),
        'kinetic energy per particle': lambda x: x.system.kinetic_energy(True),
        'total energy per particle': lambda x: x.system.total_energy(True, cache=True),
        'temperature': lambda x: x.system.temperature,
        'density': lambda x: x.system.density,
        'pressure': lambda x: x.system.pressure,
        'rmsd': lambda x: x.rmsd,
    }

    # Update db with extra functions
    if functions is not None:
        _db_func.update(functions)

    # Internal database for formats
    _db_fmt = {}
    for key in _db_func:
        # Default to float formatting
        _db_fmt[key] = '{{:.{precision}{form}}}'.format(precision=precision, form='g')
    # Steps are integer
    _db_fmt['steps'] = '{:d}'

    # Update db with extra formats
    if fmt is not None:
        _db_fmt.update(fmt)

    # Header
    if sim.current_step == 0:
        with open(sim.output_path + '.thermo', 'w') as fh:
            txt = ', '.join(fields)
            fh.write('# columns: {}\n'.format(txt))

    # Line
    with open(sim.output_path + '.thermo', 'a') as fh:
        values = [_db_func[field](sim) for field in fields]
        result = ' '.join([_db_fmt[field].format(value) for value, field in zip(values, fields)])
        fh.write('{}\n'.format(result))


def _setup_callbacks(what):
    """Setup callbacks from `what` list, see `write` for definitions"""
    from operator import attrgetter

    # Default callbacks that take simulation as first argument
    _callbacks = {
        'steps': lambda x: x.current_step,
        'potential energy per particle': lambda x: x.system.potential_energy(True),
        'kinetic energy per particle': lambda x: x.system.kinetic_energy(True),
        'total energy per particle': lambda x: x.system.total_energy(True, cache=True),
        'temperature': lambda x: x.system.temperature,
        'density': lambda x: x.system.density,
        'pressure': lambda x: x.system.pressure,
        'rmsd': lambda x: x.rmsd,
    }

    names, callbacks = [], []
    for attribute in what:
        if attribute in _callbacks:
            # A predefined callback
            names.append(attribute)
            callbacks.append(_callbacks[attribute])
        elif isinstance(attribute, list) or isinstance(attribute, tuple):
            # A tuple or list (name, callback)
            assert len(attribute) == 2
            names.append(attribute[0])
            callbacks.append(attribute[1])
        else:
            # Generic simulation attribute
            names.append(attribute)
            callbacks.append(attrgetter(attribute))

    return names, callbacks


def write(sim, what, suffix=None, path=None):
    """
    Write generic attributes of simulation `sim` to a file.

    `suffix` is a tag appended to `sim.output_path` to define the output
    file path. If `path` is provided, however, it is used instead as
    output file path.

    `what` tells the function what to write and must be a list of either:

    - a string representing a valid property of the Simulation instance `sim`

    - a tuple `(name, callback)` where `name` is a descriptive name of
    the value returned by the `callback`, which is a function that
    takes `sim` as first argument

    - a string from the following list:
    steps
    temperature
    potential energy per particle
    kinetic energy per particle
    total energy
    pressure
    rmsd
    conserved energy
    """
    # TODO: add formats / precision
    assert suffix is not None or path is not None
    assert not (suffix is not None and path is not None)

    if path is None:
        path = sim.output_path + '.' + suffix

    # Define callbacks
    names, callbacks = _setup_callbacks(what)

    # Extract the requested attribute
    values = []
    for callback in callbacks:
        values.append(callback(sim))

    # Header
    if sim.current_step == 0:
        with open(path, 'w') as fh:
            fh.write('# columns: {}\n'.format(', '.join(names)))

    # Format output string
    fmt = ('{} ' * len(values)) + '\n'
    with open(path, 'a') as fh:
        fh.write(fmt.format(*values))


def store(sim, what, db):
    """
    Store generic attributes of simulation `sim` in the dictonary `db`.

    `what` tells the function what to write and must be a list of either:

    - a string representing a valid property of the Simulation instance `sim`

    - a tuple `(name, callback)` where `name` is a descriptive name of
    the value returned by the `callback`, which is a function that
    takes `sim` as first argument

    - a string from the following list:
    steps
    temperature
    potential energy per particle
    kinetic energy per particle
    total energy
    pressure
    rmsd
    conserved energy
    """
    # TODO: allow list?
    # If the dict is empty fill it
    if len(db) == 0:
        for attribute in what:
            db[attribute] = []

    # Define callbacks
    names, callbacks = _setup_callbacks(what)

    # Extract the requested attribute
    for name, callback in zip(names, callbacks):
        db[name].append(callback(sim))

    return db


# Target callbacks.
# They should return a fractional measure of completion

def target(sim, attribute, value):
    """
    An observer that raises a `SimulationEnd` exception when a given
    target `value` of a property is reached during a simulation. The
    property is `attribute` and is assumed to be an attribute of
    simulation.

    Return: the ratio between current and target values of the attribute.
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
    return float(sim.current_step) / value

def target_walltime(sim, value):
    """
    Target a value of the elapsed wall time in seconds from the
    beginning of the simulation.

    Useful to self restarting jobs in a queining system with time
    limits.
    """
    wtime_limit = value
    if sim.wall_time() > wtime_limit:
        raise SimulationEnd('target wall time reached')
    else:
        t = sim.wall_time()
        dt = wtime_limit - t
        _log.debug('elapsed time %g, reamining time %g', t, dt)


def target_python_stop(sim, condition):
    """
    Stop the simulation if `condition` is True.

    `condition` will interpolate attributes of the passed `sim`
    instance. For instance, the condition

        #!python
        {current_step} > 1000 and {rmsd} > 1.0

    will stop the simulation when the step is > 1000 and the rmsd > 1.
    """
    # We do nothing on the first step
    if sim.current_step == 0:
        return
    # Interpolate the command string
    cmd = condition.replace('{', '{0.')
    cmd = cmd.format(sim)
    if eval(cmd):
        raise SimulationEnd('condition "{}" satisfied'.format(condition))

def shell_stop(sim, cmd, exit_code=1):
    """
    Execute the shell command `cmd` and stop the simulation if the
    command returns an exit value equal to `exit_code`.

    `cmd` is actually a format string that may contain references to
    the passed `sim` instance. For instance, a valid command is

        #!python
        echo {sim.current_step} {sim.rmsd} >> {sim.output_path}.out

    which will append the step and rmsd to {sim.output_path}.out.
    """
    import subprocess
    # We do nothing on the first step
    if sim.current_step == 0:
        return
    try:
        # Interpolate the command string
        wrap_cmd = cmd.format(sim=sim)
        # Run the shell command
        output = subprocess.check_output(wrap_cmd, shell=True,
                                         stderr=subprocess.STDOUT, executable="/bin/bash")
        if len(output) > 0:
            _log.info('shell command "{}" returned: {}'.format(cmd, output.strip()))

    except subprocess.CalledProcessError as e:
        # We stop the simulation
        if e.returncode == exit_code:
            raise SimulationEnd('shell command "{}" returned "{}"'.format(cmd, e.output.strip()))
        else:
            _log.error('shell command {} failed with output {}'.format(cmd, e.output))
            raise

def user_stop(sim):
    """
    Allows a user to stop the simulation smoothly by touching a STOP
    file in the output root directory.  Currently the file is not
    deleted to allow parallel jobs to all exit.
    """
    # To make it work in parallel we should broadcast and then rm
    # or subclass userstop in classes that use parallel execution
    if sim.output_path is not None:
        if os.path.isdir(sim.output_path):
            dirpath = sim.output_path
        else:
            dirpath = os.path.dirname(sim.output_path)
        if os.path.exists('%s/STOP' % dirpath):
            raise SimulationEnd('user has stopped the simulation')

# Aliases


target_user_stop = user_stop
target_shell_stop = shell_stop


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
            for observer in sim._observer:
                c = observer['callback']
                if c is self:
                    continue
                if 'target' in c.__name__.lower():
                    self._observer = observer
                    args = observer['args']
                    kwargs = observer['kwargs']
                    self.x_last = c(sim, *args, **kwargs)
                    self.t_last = time.time()
                    self._init = True
                    return

        t_now = time.time()
        args = self._observer['args']
        kwargs = self._observer['kwargs']
        x_now = self._observer['callback'](sim, *args, **kwargs)
        # Get the speed at which the simulation advances
        speed = (x_now - self.x_last) / (t_now - self.t_last)
        # Report fraction of target achieved and ETA
        frac = float(x_now) / 1
        try:
            eta = (1.0 - x_now) / speed
            d_now = datetime.datetime.now()
            d_delta = datetime.timedelta(seconds=eta)
            d_eta = d_now + d_delta
            # self._callback.__name__,
            _log.info('%2d%% ETA: %s S/T: %.1f T/SP: %.2e',
                      int(frac * 100),
                      d_eta.strftime('%Y-%m-%d %H.%M'),
                      1./sim.wall_time(per_step=True),
                      sim.wall_time(per_step=True, per_particle=True))
        except ZeroDivisionError:
            print(x_now, self.x_last)
            raise

        self.t_last = t_now
        self.x_last = x_now
