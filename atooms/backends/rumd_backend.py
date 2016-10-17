# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Adapters for RUMD simulation package"""

import sys
import os
import copy
from atooms import simulation
from atooms.simulation import log
from rumd import *
from rumdSimulation import rumdSimulation

class RumdBackend(object):

    # TODO: check the use case of repeated run() without restart

    def __init__(self, sim, output_path=None, fixcm_interval=0):
        self.steps = 0
        self.fixcm_interval = fixcm_interval
        self._sim = sim
        self._sim.sample.EnableBackup(False)
        self._sim.sample.SetVerbose(False)
        self._sim.SetVerbose(False)
        self._sim.SetMomentumResetInterval(self.fixcm_interval)
        self.output_path = output_path
        # TODO: need some switch to use or not RUMD checkpoint. If checkpoint_interval is set e.g. we suppress RUMD's one
        # Keep a reference of the Trajectory backend class
        self.trajectory = Trajectory
        # Copy of initial state (it is not always enough to do it in run_pre())
        self._initial_sample = self._sim.sample.Copy()
        # We set RUMD block to infinity unless the user set it already
        if self._sim.blockSize is None:
            self._sim.SetBlockSize(sys.maxint)
        # Handle output
        self._suppress_all_output = True
        self._restart = False # internal restart toggle

    def _get_system(self):
        return System(self._sim.sample)

    def _set_system(self, value): 
        self._sim.sample = value.sample

    system = property(_get_system, _set_system, 'System')

    @property
    def initial_state(self):
        return System(self._initial_sample)

    def __str__(self):
        return 'RUMD v%s' % GetVersion()

    @property
    def rmsd(self):
        """ Compute the mean square displacement between actual sample and the reference sample """
        # TODO: not sure it is the backend responsibility
        if self._sim.sample is self._initial_sample:
            raise Exception('rmsd between two references of the same system does not make sense (use deepecopy?)')
        ndim = 3 # hard coded
        N = self._sim.sample.GetNumberOfParticles()
        L = [self._sim.sample.GetSimulationBox().GetLength(i) for i in range(ndim)]
        # Unfold positions using periodic image information
        ref = self._initial_sample.GetPositions() + self._initial_sample.GetImages() * L
        unf = self._sim.sample.GetPositions() + self._sim.sample.GetImages() * L
        return (sum(sum((unf-ref)**2)) / N)**0.5

    def write_checkpoint(self):        
        if self.output_path is None:
            log.warning('output_path is not set so we cannot write checkpoint  %d' % self.steps)
            return
        f = os.path.join(self.output_path, 'trajectory.chk')
        with Trajectory(f, 'w') as t:
            t.write(self.system, None)
        with open(f + '.step', 'w') as fh:
            fh.write('%d' % self.steps)

    def read_checkpoint(self):
        f = os.path.join(self.output_path, 'trajectory.chk')
        log.debug('reading own restart file %s' % f)
        self._sim.sample.ReadConf(f)
        with open(f + '.step') as fh:
            self.steps = int(fh.read())
        log.info('backend rumd restarting from %d' % self.steps)

    def run_pre(self, restart):
        # Copy of initial state. This way even upon repeated calls to
        # run() we still get the right rmsd Note restart is handled
        # after this, so sample here is really the initial one.
        self._initial_sample = self._sim.sample.Copy()
        self._restart = restart
        self._rumd_block_index = None        

        if self.output_path is not None:
            self._sim.sample.SetOutputDirectory(self.output_path + '/rumd')

        if not self._restart:
            # We initialize RUMD writers state. Since we make 0 steps,
            # this won't create the output dir Note: RUMD complains if
            # we attempt to run some steps and we dont suppress output
            # without initializing the writers.
            self._sim.Run(0, suppressAllOutput=True, initializeOutput=True)
        else:
            log.debug('restart attempt')
            if self.output_path is None:
                log.warn('it does not make sense to restart when writing is disabled')
                return
            if os.path.exists(os.path.join(self.output_path, 'trajectory.chk')):
                # Use our own checkpoint file. We ignore RUMD checkpoint
                self.read_checkpoint()
            elif os.path.exists(self.output_path + '/rumd/LastComplete_restart.txt'):
                # Use RUMD checkpoint
                # @thomas unfortunately RUMD does not seem to write the last restart 
                # when the simulation is over therefore the last block is always rerun
                # RUMD restart file contains the block index (block_index) 
                # and the step within the block (nstep) of the most recent backup
                log.debug('reading rumd restart %s' % (self.output_path + '/rumd/LastComplete_restart.txt'))
                with open(self.output_path + '/rumd/LastComplete_restart.txt') as fh:
                    ibl, nstep = fh.readline().split()
                self._rumd_block_index = int(ibl)
                self.steps = int(nstep) * int(ibl)
                # Cleanup: delete old RUMD restart files
                import glob
                restart_files = glob.glob(self.output_path + '/rumd/restart*')
                restart_files.sort()
                for f in restart_files[:-1]:
                    log.debug('removing restart file %s' % f)
                    os.remove(f)
            else:
                log.warn('restart requested but no checkpoint is found')


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
                
        # If we use our own restart we must tell RUMD
        # to run only the difference n-self.steps. However, if we use
        # the native restart, we must keep n as is.
        log.debug('RUMD running from %d to %d steps' % (self.steps, n))
        if self._rumd_block_index is not None:
            # Restart from RUMD checkpoint
            self._sim.Run(n, restartBlock=self._rumd_block_index)
        else:
            # Not restarting or restart from our checkpoint.
            if self._restart:
                # We must toggle it here to prevent future calls to pre to restart. TODO: why??
                self._restart = False 
            self._sim.Run(n-self.steps, suppressAllOutput=self._suppress_all_output, initializeOutput=False)
            self.steps = n

import numpy
from atooms.system.particle import Particle
from atooms.system.cell import Cell


class Thermostat(object):

    """Wrap a RUMD integrator as a thermostat."""

    # TODO: problem with this approach is that we rely on RUMD keeping the same order in future versions. We should unit test it.
    # Info string looks like IntegratorNVT,0.004,0.3602,0.2,-0.7223

    def __init__(self, integrator):
        self._integrator = integrator

    def reset(self):
        info = self._integrator.GetInfoString(18).split(',')
        info[4] = '1.0'
        info = ','.join(info)
        self._integrator.InitializeFromInfoString(info)
        
    def _get_temperature(self):
        info = self._integrator.GetInfoString(18).split(',')
        return float(info[2])

    def _set_temperature(self, value):
        info = self._integrator.GetInfoString(18).split(',')
        info[2] = '%g' % value
        info = ','.join(info)
        self._integrator.InitializeFromInfoString(info)

    temperature = property(_get_temperature, _set_temperature, 'Temperature')
    

class MolecularDynamics(object):

    """Wrap an integrator as MolecularDynamics."""

    def __init__(self, integrator):
        self._integrator = integrator

    def _get_timestep(self):
        return self._integrator.GetTimeStep()

    def _set_timestep(self, value):
        self._integrator.SetTimeStep(value)

    timestep = property(_get_timestep, _set_timestep, 'Timestep')


class System(object):
    
    def __init__(self, sample):
        self.sample = sample
        self.thermostat = Thermostat(self.sample.GetIntegrator())
        self.dynamics = MolecularDynamics(self.sample.GetIntegrator())

    def __copy__(self):
        # This is not really needed, it's just there for reference
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        result.sample = self.sample.Copy()
        return result

    def __deepcopy__(self, memo):
        # TODO: @nick ask to implement
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
        # TODO: use double IntegratorNVT::GetKineticEnergy(bool copy) const{
        from atooms.system.particle import total_kinetic_energy
        return total_kinetic_energy(self.particle)

    def __get_mass(self):
        # TODO: cache it (but what if masses change?)
        npart = self.sample.GetNumberOfParticles()
        nsp = self.sample.GetNumberOfTypes()
        mass = numpy.ndarray(npart, dtype=float)
        ii = 0
        for i in range(nsp):
            ni = self.sample.GetNumberThisType(i)
            try:
                # This will work with rumd <= 2.0.1 I think
                # meta = self.sample.GetTrajectoryConfMetaData()
                # then get meta.GetMassOfType(i)
                mi = self.sample.GetMass(i)
            except:
                log.warning('cannot get mass from RUMD interface, setting to 1.0')
                mi = 1.0
            mass[ii:ii+ni] = mi
            ii += ni
        return mass

    def temperature(self):
        npart = self.sample.GetNumberOfParticles()
        ndof = self.sample.GetNumberOfDOFs()
        vel = self.sample.GetVelocities()
        mass = self.__get_mass()
        return 2 * numpy.sum(mass * numpy.sum(vel**2.0, 1))/ndof

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

    @property
    def cell(self):
        box = self.sample.GetSimulationBox()
        L = [box.GetLength(i) for i in range(3)]
        return Cell(L)

    @property
    def particle(self):
        nmap = ['A', 'B', 'C', 'D']
        n = self.sample.GetNumberOfParticles()
        pos = self.sample.GetPositions()
        vel = self.sample.GetVelocities()
        nsp = self.sample.GetNumberOfTypes()
        ima = self.sample.GetImages()
        mass = self.__get_mass()
        spe = numpy.ndarray(n, dtype=int)
        name = numpy.ndarray(n, dtype='|S1')
        ii = 0
        for i in range(nsp):
            ni = self.sample.GetNumberThisType(i)
            spe[ii:ii+ni] = i+1
            name[ii:ii+ni] = nmap[i]
            ii += ni
        p = [Particle(s, n, m, p, v) for s, n, m, p, v in zip(spe, name, mass, pos, vel)]
        for pi, i in zip(p, ima):
            pi.periodic_image = i
        return p


class Trajectory(object):

    suffix = 'xyz'

    def __init__(self, filename, mode='w'):
        self.filename = filename
        self.mode = mode

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def exclude(self, p):
        pass

    def include(self, p):
        pass

    def write(self, system, step):
        if step is None:
            f = self.filename # + '.' + self.suffix
        else:
            tag = '%011d.' % step
            if os.path.isdir(self.filename):
                f = os.path.join(self.filename, tag + self.suffix)
            else:
                f = self.filename + '.' + tag + '.' + self.suffix
        log.debug('writing config via backend to %s at step %s, %s' % (f, step, self.mode))
        system.sample.WriteConf(f, self.mode)

    def close(self):
        # This only unzips files with no step info
        if os.path.exists(self.filename + '.gz'):
            os.system("gunzip -f %s.gz" % self.filename)


def single(sim_input, potential=None, T=None, dt=0.001, interval_energy=None, interval_config=None):
    from rumd import IntegratorNVT
    from rumdSimulation import rumdSimulation

    if type(sim_input) is str:
        sim = rumdSimulation(sim_input)
        for pot in potential():
            sim.AddPotential(pot)
    else:
        sim = sim_input

    itg = IntegratorNVT(targetTemperature=T, timeStep=dt)
    sim.SetIntegrator(itg)
    sim.SetMomentumResetInterval(10000)
    sim.SetOutputScheduling("energies","none")
    sim.SetOutputScheduling("trajectory","none")
    if interval_energy is not None:
        sim.SetOutputScheduling("energies","linear",interval=interval_energy)
    if interval_config is not None:
        sim.SetOutputScheduling("trajectory","linear",interval=interval_config)

    return sim

def multi(input_file, potential, T, dt):
    from atooms.utils import size, rank, barrier
    
    # Create simulation and integrators
    for i in range(size):
        if i == rank:
            sim = [rumdSimulation(f) for f in input_file]
        barrier()
    igt = [IntegratorNVT(targetTemperature=Ti, timeStep=dti) for Ti, dti in zip(T, dt)]

    # Add potentials
    for s in sim:
        s.SetOutputScheduling("energies", "none")
        s.SetOutputScheduling("trajectory", "none")
        for p in potential():
            s.AddPotential(p)

    for s, i in zip(sim, igt):
        s.SetMomentumResetInterval(0)
        s.SetIntegrator(i)

    return sim

# Forcefields

def kalj():
    pot = rumd.Pot_LJ_12_6(cutoff_method = rumd.ShiftedPotential)
    pot.SetParams(i=0, j=0, Epsilon=1.0, Sigma=1.0, Rcut=2.5)
    pot.SetParams(i=1, j=0, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
    pot.SetParams(i=0, j=1, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
    pot.SetParams(i=1, j=1, Epsilon=0.5, Sigma=0.88, Rcut=2.5)
    return [pot]
