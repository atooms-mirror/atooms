# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Adapters for RUMD simulation package"""

import sys
import os
import logging
from atooms import simulation
from atooms.utils import mkdir
from rumd import *
from rumdSimulation import rumdSimulation

class WriterCheckpoint(object):

    def __call__(self, e):
        logging.debug('write checkpoint %d' % e.steps)
        e.write_checkpoint()

class WriterConfig(object):

    def __call__(self, e):
        logging.debug('write config %d' % e.steps)
        e.trajectory.write_sample(e.system, e.steps, e.steps)

class WriterThermo(object):

    def __call__(self, e):
        logging.debug('write thermo %d' % e.steps)
        with open(e.trajectory.filename + '.thermo', 'a') as fh:
            fh.write('%d %g %g\n' % (e.steps, e.system.potential_energy(), e.rmsd))

class Simulation(simulation.Simulation):

    _WRITER_CHECKPOINT = WriterCheckpoint
    _WRITER_CONFIG = WriterConfig
    _WRITER_THERMO = WriterThermo
    STORAGE = 'directory'

    def __init__(self, sim, dir_output):
        self._sim = sim
        self._sim.sample.SetOutputDirectory(dir_output)
        self._sim.SetMomentumResetInterval(10)
        simulation.Simulation.__init__(self, dir_output)
        self.__set_verbosity(0)
        self._initialize_output = True
        self.dir_output = dir_output
        self.trajectory = Trajectory(os.path.join(dir_output, 'config.xyz'), 'w')
        mkdir(self.dir_output)

    # Temporarily use a different method
    def __set_verbosity(self, value):
        if value == 0:
            self._sim.sample.SetVerbose(False)
            self._sim.SetVerbose(False)
            # self._sim.SetOutputScheduling("energies", "none")
            # self._sim.SetOutputScheduling("trajectory", "none")
            self._suppress_all_output = True
        else:
            self._sim.sample.SetVerbose(True)
            self._sim.SetVerbose(True)
            self._suppress_all_output = False
            # TODO: reset output parameters, or get rid of it
            
    def _get_system(self):
        return System(self._sim)

    def _set_system(self, value): 
        self._sim.sample = value.sample

    system = property(_get_system, _set_system, 'System')

    def write_checkpoint(self):
        f = self.trajectory.filename + '.chk'
        t = Trajectory(f, 'w')
        t.write_sample(self.system, None, None)
        t.close()
        with open(f + '.step', 'w') as fh:
            fh.write('%d' % self.steps)

#    def read_checkpoint(self, sim):
    def read_checkpoint(self):
        f = self.trajectory.filename + '.chk'
        self._sim.sample.ReadConf(f)
        with open(f + '.step') as fh:
            self.steps = int(fh.read())
        logging.info('restarting from %d' % self.steps)

    def __check_restart(self):
        self._ibl = None
        # TODO: this is never reset!!
        self._suppress_all_output = False
        #self._suppress_all_output = True
        # @thomas unfortunately RUMD does not seem to write the last restart 
        # when the simulation is over therefore the last block is always rerun
        # To work around it we check is final.xyz.gz exists (we write it in run_end)
        
        # This may also a problem for self restarting blocks...??
        # TODO: *** why and when does rumd moves to .bak???! ***
        if self.restart:
            
            self.restart = False

            # If we find our own checkpoint file, we ignore RUMD checkpoint
            f = self.trajectory.filename + '.chk'
            if os.path.exists(f):
                # TODO: this shouldn't be necessary because we should write checkpoint also at the end
                if os.path.exists(self.dir_output + '/final.xyz.gz'):
                    if os.path.getmtime(self.dir_output + '/final.xyz.gz') > \
                            os.path.getmtime(f):
                        # There exists a final configuration, exit immediately
                        # Update the sample from final configuartion
                        self._sim.sample.ReadConf(self.dir_output + '/final.xyz.gz')
                        raise simulation.SimulationEnd('already completed')

                self.read_checkpoint()

            else:
                # Use RUMD checkpoint
                if os.path.exists(self.dir_output + '/final.xyz.gz'):
                    if os.path.getmtime(self.dir_output + '/final.xyz.gz') > \
                            os.path.getmtime(self.dir_output + '/LastComplete_restart.txt'):
                        # There exists a final configuration, exit immediately
                        # Update the sample from final configuartion
                        self._sim.sample.ReadConf(self.dir_output + '/final.xyz.gz')
                        raise simulation.SimulationEnd('already completed')

                if os.path.exists(self.dir_output + '/LastComplete_restart.txt'):
                    with open(self.dir_output + '/LastComplete_restart.txt') as fh:
                        ibl, nstep = fh.readline().split()
                    self._ibl = int(ibl)
                    self.steps = int(nstep) * int(ibl)

                    # Delete old RUMD restart files
                    import glob
                    restart_files = glob.glob(self.dir_output + '/restart*')
                    restart_files.sort()
                    for f in restart_files[:-1]:
                        logging.debug('removing restart file %s' % f)
                        os.remove(f)
        
    def run_pre(self):
        super(Simulation, self).run_pre()
        if self._sim.blockSize is None:
            # We set RUMD block to infinity unless the user set it already
            self._sim.SetBlockSize(sys.maxint)
        logging.debug('RUMD block is %d' % self._sim.blockSize)
        self.__check_restart()
        self._initialize_output = True

    def run_until(self, n):

        # 1. suppress all RUMD output and use custom writers. PROS:
        #    this way running batches of simulations will work without
        #    rereading the restart file (everything stays in memory)
        #    CONS: we loose native RUMD output, log-lin and we have to
        #    pass through atooms trajectory (=> implement system interface)
        #    or thorugh RUMD WriteSample with names after time steps
        #    because WriteSample does not store step information!
        
        # 2. keep RUMD output. PROS: no need to recalculate things
        #    during the simulation, it should be more efficient from
        #    this point of view. CONS: we must read restart files at
        #    every batch. RUMD has a complicated output logic (blocks)

        # Large performance penalty when periods are too short in real time
                
        # TODO: why is it necessary to check this!? When running
        # successive blocks, there shouldn't be any need to restart
        # from the checkpoint because the configuration remains in
        # memory. Right, but rumd will delete the existing files and
        # start over again. This would only work if RUMD does not
        # write his files.
        #
        # *** It's a mess. If we use our own restart we must tell RUMD
        # to run only the difference n-self.steps. However, if we use
        # the native restart, we must keep n as is.

        self.__check_restart()
        if self._ibl:
            self._sim.Run(n, restartBlock=self._ibl,
                          suppressAllOutput=self._suppress_all_output)
        else:
            self._sim.Run(n-self.steps, restartBlock=self._ibl,
                          initializeOutput = self._initialize_output,
                          suppressAllOutput=self._suppress_all_output)
            #self._initialize_output = False
            self.steps = n

    def run_end(self):
        super(Simulation, self).run_end()
        # Make sure we write final.xyz.gz, this way we can avoid
        # restarting from the but to last block
        self._sim.WriteConf(self.dir_output + '/final.xyz.gz')

import numpy
from atooms.system.particle import Particle

class System(object):
    
    def __init__(self, sim):
        self._sim = sim
        self.sample = sim.sample

    def __copy__(self):
        # This is not really needed, it's just there for reference
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        result.sample = self.sample.Copy()
        return result
        
    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        result.__dict__.update(self.__dict__)
        # Use Copy() method of sample, 
        result.sample = self.sample.Copy()
        # We do not copy recursively, deepcopy fails when wrapping SWIG classes
        # from copy import deepcopy
        # for k, v in self.__dict__.items():
        #     setattr(result, k, deepcopy(v, memo))
        return result

    def potential_energy(self):
        self.sample.CalcF()
        return self.sample.GetPotentialEnergy()

    def kinetic_energy(self):
        return ekin

    def temperature(self):
        npart = self.sample.GetNumberOfParticles()
        vel = self.sample.GetVelocities()
        # TODO: missing masses in T !
        T = sum(sum(vel**2.0))/(3*npart - 3)
        return T

    def mean_square_displacement(self, reference):
        """ Compute the mean square displacement between actual sample and the reference sample """
        if reference.sample is self.sample:
            raise Exception('rmsd between two references of the same system does not make sense (use deepecopy?)')

        ndim = 3 # hard coded
        N = self.sample.GetNumberOfParticles()
        L = [self.sample.GetSimulationBox().GetLength(i) for i in range(ndim)]

        # Unfold positions using periodic image information
        ref = reference.sample.GetPositions() + reference.sample.GetImages() * L
        unf = self.sample.GetPositions() + self.sample.GetImages() * L

        return sum(sum((unf-ref)**2)) / N

    def _get_thermostat(self):
        return self.sample.GetIntegrator()

    def _set_thermostat(self, value): 
        self._sim.SetIntegrator(value)

    thermostat = property(_get_thermostat, _set_thermostat, 'Thermostat (actually an instance of RUMD Integrator')

    @property
    def particle(self):
        nmap = ['A', 'B', 'C', 'D']
        # TODO: get mass from rumd!
        #meta = self.sample.GetTrajectoryConfMetaData()
        n = self.sample.GetNumberOfParticles()
        pos = self.sample.GetPositions()
        vel = self.sample.GetVelocities()
        nsp = self.sample.GetNumberOfTypes()
        spe = numpy.ndarray(n, dtype=int)
        name = numpy.ndarray(n, dtype='|S1')
        mass = numpy.ndarray(n, dtype=float)
        ii = 0
        for i in range(nsp):
            ni = self.sample.GetNumberThisType(i)
            #mi = meta.GetMassOfType(i)
            mi = 1.0
            spe[ii:ii+ni] = i+1
            name[ii:ii+ni] = nmap[i]
            mass[ii:ii+ni] = mi
        p = [Particle(s, n, m, p, v) for s, n, m, p, v in zip(spe, name, mass, pos, vel)]
        return p

class Trajectory(object):

    def __init__(self, filename, mode='r'):
        self.filename = filename
        self.mode = 'w'

        # Remove the .gz extension from path. It will be added by rumd anyway
        base, ext = os.path.splitext(self.filename)
        if ext == '.gz':
            self.filename = base

    def write_initial_state(self, system):
        pass

    def write_sample(self, system, step, sample):
        f = self.filename 
        if step:
            f += '_%011d' % step
        system.sample.WriteConf(f, self.mode)

        # Next time we append
        # if self.mode == 'w':
        #     self.mode = 'a'

    # def read_sample(self, system, sample):
    #     return system.sample.ReadConf(self.filename)            

    def close(self):
        # Assuming something has been written, unzip the trajectory file
        if os.path.exists(self.filename + '.gz'):
            os.system("gunzip -f %s.gz" % self.filename)

    
