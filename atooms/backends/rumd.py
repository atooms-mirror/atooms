# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""Simulation backend for RUMD (http://rumd.org/)."""

# This enables importing the top-level rumd package and still call
# this file rumd.py
from __future__ import absolute_import

import sys
import os
import numpy
import logging

import rumd  # pylint:disable=import-error
from rumd.Simulation import Simulation  # pylint:disable=import-error
from atooms.system import System as _System
from atooms.system import Particle, Cell

_log = logging.getLogger(__name__)
_version = rumd.GetVersion()

if int(_version.split('.')[0]) < 3:
    raise ImportError('RUMD version must be at least 3')


def unfold(system):
    # s = system
    # particle = system.particle
    # for i, p in enumerate(particle):
    #     p.position += p.periodic_image * s.cell.side
    # system.particle = particle
    # return s
    npart = system.sample.GetNumberOfParticles()
    pos = system.sample.GetPositions()
    nsp = system.sample.GetNumberOfTypes()
    ima = system.sample.GetImages()
    spe = numpy.ndarray(npart, dtype=int)
    ii = 0
    for i in range(nsp):
        ni = system.sample.GetNumberThisType(i)
        spe[ii: ii + ni] = i
        ii += ni
    particle = [Particle(species=spe_i, position=pos_i) for spe_i, pos_i in
                zip(spe, pos)]
    for p, i in zip(particle, ima):
        p.position += i * system.cell.side

    return _System(particle=particle, cell=system.cell)


class RUMD(object):

    """RUMD simulation backend."""

    version = _version

    def __init__(self, input_file_or_sim_or_system, potentials=None,
                 integrator=None, temperature=None, dt=0.001,
                 fixcm_interval=0, thermostat_relaxation_time=None):

        # Keep a reference of the Trajectory backend class
        self.trajectory_class = Trajectory
        self.timestep = dt

        # If we pass a System instance, we write a temporary RUMD
        # trajectory to initialize the simulation and reset the
        # input_file_or_sim_or_system variable
        if isinstance(input_file_or_sim_or_system, _System):
            import tempfile, os
            from atooms.trajectory import TrajectoryRUMD            
            tmp = os.path.join(tempfile.mkdtemp(), 'config.xyz')
            with TrajectoryRUMD(tmp, 'w') as th:
                th.write(input_file_or_sim_or_system)
            input_file_or_sim_or_system = tmp
        
        # Store internal rumd simulation instance.
        # It is exposed as self.rumd_simulation for further customization
        if isinstance(input_file_or_sim_or_system, Simulation):
            self.rumd_simulation = input_file_or_sim
            self._suppress_all_output = False
            self._initialize_output = True

        else:
            self.rumd_simulation = Simulation(input_file_or_sim_or_system, verbose=False)
            self.rumd_simulation.SetVerbose(False)
            self.rumd_simulation.sample.SetVerbose(False)
            self.rumd_simulation.sample.EnableBackup(False)
            self.rumd_simulation.SetMomentumResetInterval(fixcm_interval)
            self.rumd_simulation.SetBlockSize(sys.maxsize)
            self.rumd_simulation.write_timing_info = False

            # By default we mute RUMD output.
            self.rumd_simulation.SetOutputScheduling("energies", "none")
            self.rumd_simulation.SetOutputScheduling("trajectory", "none")
            self._suppress_all_output = True
            self._initialize_output = False

        # Set the forcefield
        if potentials is not None:
            # We are provided a list of rumd potentials
            for potential in potentials:
                self.rumd_simulation.AddPotential(potential)

        # Add a rumd integrator
        if temperature is not None:
            integrator = 'nvt'

        if integrator is not None:
            if integrator in ['nvt', 'NVT']:
                itg = rumd.IntegratorNVT(targetTemperature=temperature,
                                         timeStep=dt)
                if thermostat_relaxation_time is not None:
                    itg.SetRelaxationTime(float(thermostat_relaxation_time))

            elif integrator in ['nve', 'NVE']:
                itg = rumd.IntegratorNVE(timeStep=dt)
            self.rumd_simulation.SetIntegrator(itg)

        # Copy of initial state
        self._initial_sample = self.rumd_simulation.sample.Copy()

        # Hold a reference to the system
        # self.system = System(self.rumd_simulation.sample, self.rumd_simulation.potentialList)

    # Wrapping system is needed because rumd holds a reference to the
    # potentials in rumd_simulation and they are needed to create a
    # working sample from scratch
    #
    # It appears mostly necessary because we must operate on the underlying sample.
    #
    def _get_system(self):
        system = System(self.rumd_simulation.sample)
        # system.__itg_infoStr_start = self.rumd_simulation.itg.GetInfoString(8)
        return system

    def _set_system(self, value):
        self.rumd_simulation.sample.Assign(value.sample)
        # TODO: improve copying over of thermostat
        # self.rumd_simulation.itg.InitializeFromInfoString(value._itg_infoStr_start)
        # Setting sample this way is useless.
        #   self.rumd_simulation.sample = value.sample
        # Rumd actually sets samples via a file, there seems to be no other way.
        # TODO: to retain modifications to system, use atooms trajectory but at the moment we would loose info on thermostat
        # import tempfile
        # from atooms.core.utils import rmd
        # # Why should we set the output dir? It should not change
        # #tmp = value.sample.GetOutputDirectory()
        # dirout = tempfile.mkdtemp()
        # file_tmp = os.path.join(dirout, 'sample.xyz.gz')
        # value.sample.WriteConf(file_tmp)
        # self.rumd_simulation.sample.ReadConf(file_tmp, init_itg=False)
        # # Why should we set the output dir? It should not change
        # # value.sample.SetOutputDirectory(tmp)
        # # Clean up
        # rmd(dirout)

    system = property(_get_system, _set_system, 'System')

    def __str__(self):
        return 'RUMD'

    @property
    def rmsd(self):
        """
        Compute the mean square displacement between actual sample and the
        reference sample.
        """
        if self.rumd_simulation.sample is self._initial_sample:
            raise Exception('rmsd between two references of the same system does not make sense (use deepecopy?)')
        ndim = 3  # hard coded
        N = self.rumd_simulation.sample.GetNumberOfParticles()
        L = [self.rumd_simulation.sample.GetSimulationBox().GetLength(i) for i in range(ndim)]
        # Unfold positions using periodic image information
        ref = self._initial_sample.GetPositions() + self._initial_sample.GetImages() * L
        unf = self.rumd_simulation.sample.GetPositions() + self.rumd_simulation.sample.GetImages() * L
        return (sum(sum((unf - ref)**2)) / N)**0.5

    def write_checkpoint(self, output_path):
        with Trajectory(output_path + '.chk', 'w') as t:
            t.write(self.system, None)

    def read_checkpoint(self, output_path):
        if os.path.exists(output_path + '.chk'):
            self.rumd_simulation.sample.ReadConf(output_path + '.chk')
        else:
            _log.debug('could not find checkpoint')

    def run(self, steps):
        self.rumd_simulation.Run(steps,
                                 suppressAllOutput=self._suppress_all_output,
                                 initializeOutput=self._initialize_output)
        self._initialize_output = False

class Thermostat(object):

    """Wrap a RUMD integrator as a thermostat."""

    # TODO: problem, RUMD must keep the same order in future versions
    # We should unit test this
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


class System(object):

    """System wrapper for RUMD."""

    def __init__(self, sample):
        self.sample = sample
        self.thermostat = Thermostat(self.sample.GetIntegrator())
        self.barostat = None
        self.reservoir = None

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
        # result.__dict__.update(self.__dict__)  #these are refs...
        # Copy() does not copy the integrator
        result.sample = self.sample.Copy()
        # This way we set the integrator. Note that it is always the same object...
        # result.sample.SetIntegrator(self.sample.GetIntegrator())
        # result._itg_infoStr_start = self.thermostat._integrator.GetInfoString(18)
        # TODO: thermostat should be a property this way we would not need to update this
        result.thermostat = Thermostat(self.sample.GetIntegrator())
        return result

    def update(self, other):
        self.sample.Assign(other.sample)
        # self.sample.SetIntegrator(other.sample.GetIntegrator())  # maybe not necessary
        self.thermostat = Thermostat(self.sample.GetIntegrator())

    def potential_energy(self, per_particle=False, normed=False, cache=False):
        self.sample.CalcF()
        if normed or per_particle:
            return self.sample.GetPotentialEnergy() / len(self.particle)
        else:
            return self.sample.GetPotentialEnergy()

    def kinetic_energy(self, per_particle=False, normed=False):
        # TODO: use double IntegratorNVT::GetKineticEnergy(bool copy) const{
        ekin = sum([p.kinetic_energy for p in self.particle])
        if normed or per_particle:
            return ekin / len(self.particle)
        else:
            return ekin

    def total_energy(self, per_particle=False, normed=False, cache=False):
        return self.potential_energy(per_particle=per_particle, normed=normed, cache=cache) +\
            self.kinetic_energy(per_particle=per_particle, normed=normed)

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
                _log.warn('cannot get mass from RUMD interface, setting to 1.0')
                mi = 1.0
            mass[ii: ii + ni] = mi
            ii += ni
        return mass

    @property
    def temperature(self):
        ndof = self.sample.GetNumberOfDOFs()
        vel = self.sample.GetVelocities()
        mass = self.__get_mass()
        return numpy.sum(mass * numpy.sum(vel**2.0, 1)) / ndof

    # TODO: this should be the attribute setter
    def set_temperature(self, T):
        # Scale velocities from temperature Told to T
        # TODO: use maxwellian
        # TODO: subtract CM velocity
        Told = self.temperature
        velocity_factor = (T/Told)**0.5
        self.sample.ScaleVelocities(velocity_factor)

    def scale_velocities(self, factor):
        self.sample.ScaleVelocities(factor)

    @property
    def cell(self):
        box = self.sample.GetSimulationBox()
        L = [box.GetLength(i) for i in range(3)]
        return Cell(L)

    @property
    def particle(self):
        # Warning n.1 : this is read only. If you change the particles, the
        # modification won't be propoagated to the RUMD objects.
        # One would have to create a new system.
        #
        # Warning n.2 : it ia assumed that particles are sorted by species.
        # since RUMD does not have accessors to particle types (why??)
        # and we can only access the number of particles of a given type.
        npart = self.sample.GetNumberOfParticles()
        pos = self.sample.GetPositions()
        vel = self.sample.GetVelocities()
        nsp = self.sample.GetNumberOfTypes()
        ima = self.sample.GetImages()
        mass = self.__get_mass()
        spe = numpy.ndarray(npart, dtype=int)
        ii = 0
        for i in range(nsp):
            ni = self.sample.GetNumberThisType(i)
            spe[ii: ii + ni] = i
            ii += ni
        p = [Particle(species=spe_i, mass=mass_i, position=pos_i,
                      velocity=vel_i) for spe_i, mass_i, pos_i, vel_i
             in zip(spe, mass, pos, vel)]
        for pi, i in zip(p, ima):
            pi.periodic_image = i
        return p

    def dump(self, what):
        import atooms.system
        system = atooms.system.System(self.particle, self.cell)
        return system.dump(what)


class Trajectory(object):

    suffix = 'xyz'

    def __init__(self, filename, mode='w'):
        self.filename = filename
        self.mode = mode

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def write(self, system, step):
        """
        If step is not None, output will follow a folder-based logic and
        filename will be considered as the root folder
        """
        if step is None:
            f = self.filename
        else:
            fbase = '%011d.%s' % (step, self.suffix)
            f = os.path.join(self.filename, fbase)
            if not os.path.exists(self.filename):
                os.makedirs(self.filename)
        _log.debug('writing config via backend to %s at step %s, %s', f, step, self.mode)
        system.sample.WriteConf(f, self.mode)

    def close(self):
        # This only unzips files with no step info
        if os.path.exists(self.filename + '.gz'):
            os.system("gunzip -f %s.gz" % self.filename)
