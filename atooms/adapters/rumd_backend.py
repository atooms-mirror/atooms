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


class WriterCheckpoint(object):

    def __call__(self, e):
        log.debug('write checkpoint %d' % e.steps)
        e.write_checkpoint()

class WriterConfig(object):

    def __call__(self, e):
        log.debug('write config %d' % e.steps)
        with Trajectory(e.output_file, 'a') as t:
            t.write(e.system, e.steps)

class WriterThermo(object):

    def __call__(self, e):
        log.debug('write thermo %d' % e.steps)
        with open(e.output_file + '.thermo', 'a') as fh:
            fh.write('%d %g %g\n' % (e.steps, e.system.potential_energy(), e.rmsd))

# TODO: can we have backend not inherit from simulation base class? 
# Strategy is better than inheritance. What is the minimal subset of methods / attributes that
# need to be exposed to Simulation and its subclasses, if we only were to implement it as a backend?

class Simulation(simulation.Simulation):

    _WRITER_CHECKPOINT = WriterCheckpoint
    _WRITER_CONFIG = WriterConfig
    _WRITER_THERMO = WriterThermo
    STORAGE = 'directory'

    def __init__(self, sim, *args, **kwargs):
        # System is an instance in base class, but this adapter redefines it as a property
        # So perhaps we might pass None and make an explicit copy of initial_state here
        self._sim = sim
        simulation.Simulation.__init__(self, self.system, *args, **kwargs)
        self._sim.sample.EnableBackup(False)
        self._sim.sample.SetVerbose(False)
        self._sim.SetVerbose(False)
        # TODO: need some switch to use or not RUMD checkpoint. If checkpoint_interval is set e.g. we suppress RUMD's one
        # Copy of initial state
        self.initial_state = self._sim.sample.Copy()

    def _get_system(self):
        return System(self._sim)

    def _set_system(self, value): 
        self._sim.sample = value.sample

    system = property(_get_system, _set_system, 'System')

    def __str__(self):
        return 'RUMD %s' % GetVersion()        

    @property
    def rmsd(self):
        """ Compute the mean square displacement between actual sample and the reference sample """
        if self._sim.sample is self.initial_state:
            raise Exception('rmsd between two references of the same system does not make sense (use deepecopy?)')
        ndim = 3 # hard coded
        N = self._sim.sample.GetNumberOfParticles()
        L = [self._sim.sample.GetSimulationBox().GetLength(i) for i in range(ndim)]

        # Unfold positions using periodic image information
        ref = self.initial_state.GetPositions() + self.initial_state.GetImages() * L
        unf = self._sim.sample.GetPositions() + self._sim.sample.GetImages() * L
                
        return (sum(sum((unf-ref)**2)) / N)**0.5

    def write_checkpoint(self, filename=None):
        if filename is None:
            filename = self.trajectory.filename + '.chk'
        with Trajectory(filename, 'w') as t:
            t.write(self.system, None)
        with open(filename + '.step', 'w') as fh:
            fh.write('%d' % self.steps)

    def read_checkpoint(self, filename=None):
        if filename is None:
            filename = self.trajectory.filename + '.chk'
        log.debug('reading own restart file %s' % filename)
        self._sim.sample.ReadConf(filename)
        with open(filename + '.step') as fh:
            self.steps = int(fh.read())
        log.info('rumd restarting from %d' % self.steps)

    def __check_restart(self):
        self._ibl = None
        if self.output_path is None:
            return
        # @thomas unfortunately RUMD does not seem to write the last restart 
        # when the simulation is over therefore the last block is always rerun
        # To work around it we check is final.xyz.gz exists (we write it in run_end)        
        # TODO: This may also a problem for self restarting blocks...??
        if self.restart:
            log.debug('restart attempt')
            # TODO: we assume this is the name of the checkpoint! but we had introduced 
            # an optional filename, which will be skipped here. One soultion is to be rigid
            # but then it means that we swap trajectories in PT (need refactoring)
            # Seriously, what does all the above mean?? 20.02.2016
            if os.path.exists(self.trajectory.filename + '.chk'):
                # If we find our own checkpoint file, we ignore RUMD checkpoint
                self.read_checkpoint()
            else:
                # Use RUMD checkpoint
                if os.path.exists(self.output_path + '/final.xyz.gz'):
                    if os.path.getmtime(self.output_path + '/final.xyz.gz') > \
                            os.path.getmtime(self.output_path + '/LastComplete_restart.txt'):
                        # Update the sample from final configuartion
                        self._sim.sample.ReadConf(self.output_path + '/final.xyz.gz')
                elif os.path.exists(self.output_path + '/LastComplete_restart.txt'):
                    log.debug('reading rumd restart file %s' % (self.output_path + '/LastComplete_restart.txt'))
                    with open(self.output_path + '/LastComplete_restart.txt') as fh:
                        ibl, nstep = fh.readline().split()
                    self._ibl = int(ibl)
                    self.steps = int(nstep) * int(ibl)

                    # Delete old RUMD restart files
                    import glob
                    restart_files = glob.glob(self.output_path + '/restart*')
                    restart_files.sort()
                    for f in restart_files[:-1]:
                        log.debug('removing restart file %s' % f)
                        os.remove(f)
        
    def run_pre(self):
        super(Simulation, self).run_pre()
        # If by now we havent set output_path it means we wont write anything to disk
        # This will avoid RUMD restart files to pollute the folder, however, it must be kept False the first time RUMD is run to at least get rid of the old directory. Or we do it manually
        self._suppress_all_output = True
        if self.output_path is not None:
            self._suppress_all_output = False # we let it clean the dir upon entrance, then we suppress
            self._sim.sample.SetOutputDirectory(self.output_path)
            # TODO: output_file can be dropped as instance variable I guess: we cacn check output_path is not None, grab the basename of trajectory. Or perhaps this should be defined by the writer, no?
            self.output_file = os.path.join(self.output_path, 'config.xyz')
            # We need to keep a reference to the trajectory backend used here
            # TODO: do we really need to create a file for that? Cant we we just inspect the backend? Can we just keep a reference of the class?
            self.trajectory = Trajectory(self.output_file, 'w')
            self.trajectory.close()

        if self._sim.blockSize is None:
            # We set RUMD block to infinity unless the user set it already
            self._sim.SetBlockSize(sys.maxint)

        log.debug('RUMD block is %d' % self._sim.blockSize)

        self.__check_restart()
        # Every time we call run() and we are not restarting, we clear output_path
        # TODO: check the use case of repeated run() without restart
        # Initialize output in RUMD is checked upon calling Run()
        # and takes care of creating backup if requested. If no backup is done and this is True
        # the directory is deleted. That's why this variable must be set to false
        # after the first call to Run()
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
                
        # If we use our own restart we must tell RUMD
        # to run only the difference n-self.steps. However, if we use
        # the native restart, we must keep n as is.
        # TODO: what is ibl? Use expressive names!
        if self._ibl is not None:
            # We are restarting from RUMD checkpoint. We don't need
            # to worry about initializeOutput.
            # TODO: in this shouldnt we keep all output?
            self._sim.Run(n, restartBlock=self._ibl,
                          suppressAllOutput=self._suppress_all_output)
        else:
            # Either we are not restarting or we restart
            # from our checkpoint. In the first case we make sure
            # the directory is cleared the first time this is called
            # TODO: Actually, we could avoid the inconvenience and do it on our
            # own (in this case set initializeOutput fo False)            
            if self.restart:
                self._suppress_all_output = True
                self._initialize_output = False
                # We must toggle it here to prevent future calls to pre to restart
                # TODO: improve restart here (we need an internal variable)
                self.restart = False 
            self._sim.Run(n-self.steps, suppressAllOutput=self._suppress_all_output,
                          initializeOutput = self._initialize_output)
            # We let it clean the dir upon entrance, then we suppress this
            # Note this will leave a restart file.
            # TODO: perhaps it is better to clean up the dir on our own.
            self._suppress_all_output = True
            # After calling run_until once, we prevent RUMD from
            # clearing the directory (note we have disabled backups but
            # still it clears the output_path)
            self._initialize_output = False
            self.steps = n

    def run_end(self):
        # Make sure we write final.xyz.gz, this way we can avoid
        # restarting from the but to last block
        if self.output_path is not None:
            self._sim.WriteConf(self.output_path + '/final.xyz.gz')

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
        return info[2]

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
    
    def __init__(self, sim):
        self._sim = sim
        self.sample = sim.sample
        # TODO: system may not have one right now, what will this give? A None?
        self.thermostat = Thermostat(self.sample.GetIntegrator())
        self.dynamics = MolecularDynamics(self.sample.GetIntegrator())

    def __copy__(self):
        # This is not really needed, it's just there for reference
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        result.sample = self.sample.Copy()
        return result

    # TODO: @nick ask to implement
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
        return p

class Trajectory(object):

    def __init__(self, filename, mode='w'):
        self.filename = filename
        self.mode = mode

        # Remove the .gz extension from path. It will be added by rumd anyway
        base, ext = os.path.splitext(self.filename)
        if ext == '.gz':
            self.filename = base

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def exclude(self, p):
        pass

    def include(self, p):
        pass

    def write(self, system, step):
        f = self.filename 
        if step:
            base, ext = os.path.splitext(f)
            f = base + '_%011d' % step + ext
        system.sample.WriteConf(f, self.mode)

    def close(self):
        # Assuming something has been written, unzip the trajectory file
        # TODO: this won't unzip config_*gz files
        if os.path.exists(self.filename + '.gz'):
            os.system("gunzip -f %s.gz" % self.filename)
    
def single(input_file, potential, T, dt, interval_energy=None, interval_config=None):
    from rumd import IntegratorNVT
    from rumdSimulation import rumdSimulation

    sim = rumdSimulation(input_file)
    for pot in potential():
        sim.AddPotential(pot)

    itg = IntegratorNVT(targetTemperature=T, timeStep=dt)
    sim.SetIntegrator(itg)
    sim.SetMomentumResetInterval(10000)
    sim.SetOutputScheduling("energies","none")
    sim.SetOutputScheduling("trajectory","none")
    if interval_energy:
        sim.SetOutputScheduling("energies","linear",interval=interval_energy)
    if interval_config:
        sim.SetOutputScheduling("trajectory","linear",interval=interval_config)

    # Output dir??
    return sim
    #yield si

def multi(input_file, potential, T, dt):
    # from rumd import IntegratorNVT
    from atooms.utils import size, rank, barrier
    # from rumdSimulation import rumdSimulation
    
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
