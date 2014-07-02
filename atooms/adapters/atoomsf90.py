# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich


"""Atooms simulation package in fortran 90.

We rely on command line tools 
"""

import os
import copy
import subprocess
from atooms import trajectory
from atooms import simulation
from atooms import system
from atooms.system import thermostat

def _get_opts(opts):
    s = ''
    for i in opts:
        s += ' %s %s' % (i, opts[i])
    return s

class System(system.System):

    def __init__(self, filename, opts={}):
        if os.path.exists(filename):
            try:
                t = trajectory.TrajectoryHDF5(filename)
                s = t.read_initial_state()
                s = t.read_sample(t.samples[-1])
                t.close()
            except:
                print 'error with ', filename
                raise
            super(System, self).__init__(s.particle, s.cell, s.interaction, matrix=s.matrix, thermostat=s.thermostat)
        else:
            super(System, self).__init__()

        self.filename = filename

        # Some command line options may be needed to initialize the system.
        # For instance, if a thermostat is specified we must update the attribute in system
        # For the moment we only retain the temperature
        # TODO: actually this could go into a ThermostatAtooms
        self.opts = opts
        if '--thermostat' in self.opts:
            t = thermostat.Thermostat(self.opts['--thermostat'], self.opts['--thermostat-temperature'])
            if '--thermostat-mass' in self.opts:
                t.mass = self.opts['--thermostat-mass']
            if '--thermostat-period' in self.opts:
                t.collision_period = self.opts['--thermostat-period']
               
    def potential_energy(self):
        """ Full calculation of potential energy from file """
        cmd = '/home/coslo/codes/atooms/bin/energy.x -f -1 %s' % self.filename
        out = subprocess.check_output(cmd, shell=True).split()[2]
        return float(out)

    def force(self, seed=0):
        """ Full calculation of forces from file """
        cmd = '/home/coslo/codes/atooms/bin/force.x -s %d  -n -1 %s' % (seed, self.filename)
        out = subprocess.check_output(cmd, shell=True).split()
        return [float(o) for o in out]

    def potential_energy_tail(self):
        """ Return last line of thermo file """
        cmd = 'tail -1 %s.thermo' % self.filename
        out = subprocess.check_output(cmd, shell=True).split()[1]
        return float(out)

    def mean_square_displacement(self, reference):
        raise NotImplementedYet('rmsd missing from atooms')

# TODO: each simulation backend should have a storage attribute which can be either directory or file
# depending on how files will be arranged. Otherwise it should always decide internally what to do.
# Perhaps we can keep redundancy, but it seems that directories output is the most general approach
# However this would fail if many simulations are kept in the same directory.

class Simulation(simulation.Simulation):

    STORAGE = 'file'

    def __init__(self, file_output, file_input=None, opts={}):
        super(Simulation, self).__init__(file_output)
        self.file_output = file_output
        self.file_input = file_input
        if file_input is None:
            self.file_input = self.file_output
        self.opts = opts

        self.file_output_tmp = file_output + '.tmp'

        # TODO: should initial state be an input variable or just an entry in the opts dict?
        self.opts['--initial-state'] = file_input
        self.verbosity = 0
        self.system = System(self.file_output_tmp, self.opts) #, {k: self.opts[k] for k in ('--temperature',)})

    @property
    def trajectory(self):
        try:
            # TrajectoryAtooms will delay opening / writing until
            # we actually call write_sample. This behavior might be the best
            # in general (actually this is what atooms does, no?)
            #return TrajectoryAtooms(self.file_output, 'w')
            return trajectory.TrajectoryHDF5(self.file_output, 'w')
        except:
            print self.file_output
            raise

    def _update_thermostat(self, system):
        """
        Check whether system has modifications on the thermostat and update 
        the command line options accordingly
        """
        if system.thermostat:
            self.opts.update({'--temperature':system.thermostat.temperature})

    def run(self):
        if self.verbosity == 0:
            # By setting the write config period equal to nsteps
            # we'll only dump the first and last configurations.
            self.opts['-c'] = self.target_steps

        # Convert thermostat options from system
        self._update_thermostat(self.system)

        cmd = 'md.x -n %d' % self.target_steps
        cmd += _get_opts(self.opts)
        cmd += ' %s' % self.file_output_tmp
        out = subprocess.check_output(cmd, shell=True)

        # Update internal reference system
        self.system = copy.deepcopy(System(self.file_output_tmp))

        # Next time we run we'll restart from the current tmp config file
        if '--initial-state' in self.opts:
            del self.opts['--initial-state']
       
    # This won't work too well because trajectories will get opened
    # by the running md.x ... we could rewrite a lazy version of
    # Trajectory which only opens the file when requested, adjusting
    # the mode to the actual call
    # def __init__(self, trajectory, trajectory_input, opts):
    #     self.trajectory_input = trajectory_input
    #     self.trajectory = trajectory
    #     self.opts = opts
    #     self.opts['--initial-state'] = trajectory_input.filename
    #     self.verbosity = 0
