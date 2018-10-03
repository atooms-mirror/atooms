# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""HDF5 trajectory format."""

import numpy
import h5py
import logging
import copy

from .base import TrajectoryBase
from atooms.core import ndim
from atooms.system import System
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.interaction import Interaction
from atooms.interaction.potential import PairPotential
from atooms.interaction.cutoff import CutOff

_log = logging.getLogger(__name__)


#  Helper functions and classes

class _SafeFile(h5py.File):
    # TODO: decorate hdf5 class so that error messages contain the path of the offending file
    def create_group_safe(self, group):
        # TODO: recursively create all h5 groups that do not exist?
        # TODO: redefine create_group unless this is a serious performace issue?
        if group not in self:
            self.create_group(group)

def _get_cached_list_h5(fh, h5g, data):
    """ Replace |data| with read data in group |h5g| only if data is None """
    if data is None:
        try:
            data = [entry[0] for entry in fh[h5g].values()]
        except:
            data = []
    return data

def _write_datasets(fh, group, datasets):
    """Write several data sets stored in datasets dict in group of fh"""
    for name, dataset in datasets.items():
        fh[group + name] = dataset

def add_interaction_hdf5(finp, ff):
    """Add interaction to hdf5 file"""
    import os
    import glob

    pid = os.getpid()
    f_ref = '/tmp/cnv_%s.h5' % pid
    # TODO: we can cache a ref file if ff is the same
    os.system('system.x -n 2 -f %s %s 1>/dev/null 2>/dev/null' % (ff, f_ref))
    ref = h5py.File(f_ref, 'r')
    fout = finp + '.bak'

    # Add interaction
    os.system('/bin/cp %s %s' % (finp, fout))
    h5 = h5py.File(fout, 'r+')
    # Make sure interaction does not exist
    try:
        del h5['initialstate/interaction']
    except:
        pass
    h5.copy(ref['initialstate/interaction'], 'initialstate/interaction')

    # Final cleanup
    h5.close()
    os.remove(finp)
    os.rename(fout, finp)
    ref.close()
    for f in glob.glob(f_ref + '*'):
        os.remove(f)


class TrajectoryHDF5(TrajectoryBase):

    """Trajectory layout based on HDF5 library. """

    suffix = 'h5'

    def __init__(self, filename, mode='r+'):
        TrajectoryBase.__init__(self, filename, mode)

        self.general_info = {}
        self._grandcanonical = False
        self._system = None
        self.fields = ['position', 'velocity', 'cell']

        if self.mode == 'r' or self.mode == 'r+':
            self.trajectory = h5py.File(self.filename, mode)
            # gather general info on file
            for entry in self.trajectory['/']:
                if type(self.trajectory[entry]) == h5py.highlevel.Dataset:
                    self.general_info[entry] = self.trajectory[entry]

        elif self.mode == 'w' or self.mode == 'r+' or self.mode == "w-":
            self.trajectory = _SafeFile(self.filename, self.mode)

        else:
            raise ValueError('Specify mode (r/w) for file %s (invalid: %s)' % (self.filename, self.mode))

    def read_steps(self):
        return [d[0] for d in self.trajectory['trajectory/realtime/stepindex'].values()]

    def read_len(self):
        return len(self.trajectory['trajectory/realtime/stepindex'])

    def close(self):
        try:
            self.trajectory.close()
        except ValueError:
            _log.error('file %s already closed', self.filename)
            raise

    def read_timestep(self):
        try:
            return self.trajectory['trajectory/realtime/timestep'][0]
        except:
            _log.warning('no time step in file %s, set to 1', self.filename)
            return 1.0

    def write_timestep(self, value):
        self.trajectory.create_group_safe('/trajectory')
        self.trajectory.create_group_safe('/trajectory/realtime')
        self.trajectory['trajectory/realtime/timestep'] = [value]

    def read_block_size(self):
        try:
            return self.trajectory['trajectory/realtime/block_period'][0]
        except:
            return None

    def write_block_size(self, value):
        self.trajectory.create_group_safe('/trajectory')
        self.trajectory.create_group_safe('/trajectory/realtime')
        self.trajectory['trajectory/realtime/block_period'] = [value]

    def write_init(self, system):
        from atooms.system.particle import distinct_species
        self.trajectory.create_group_safe('/initialstate')
        self.trajectory['DIMENSIONS'] = [3]
        self.trajectory['NAME_SYS'] = ['Unknown']
        self.trajectory['VERSION_TRJ'] = ['1.2']
        self.trajectory['VERSION_MD'] = ['X.X.X']

        # Particles
        group = '/initialstate/particle/'
        if system.particle is not None:
            self.trajectory.create_group_safe(group)
            particle = system.particle
            species = distinct_species(particle)
            particle_h5 = {'number_of_species': [len(species)],
                           'number_of_particles': [len(particle)],
                           'identity': [species.index(p.species)+1 for p in particle],
                           'element': ['%3s' % p.species for p in particle],
                           'mass': [p.mass for p in particle],
                           'radius': [p.radius for p in particle],
                           'position': [p.position for p in particle],
                           'velocity': [p.velocity for p in particle],
                           }
            _write_datasets(self.trajectory, group, particle_h5)

        # Matrix
        group = '/initialstate/matrix/'
        if system.matrix is not None:
            self.trajectory.create_group_safe(group)
            matrix = system.matrix
            species = distinct_species(matrix)
            matrix_h5 = {'type': [''],
                         'id': [0],
                         'number_of_species': [len(species)],
                         'number_of_particles': [len(matrix)],
                         'identity': [species.index(p.species)+1 for p in matrix],
                         'element': ['%3s' % p.species for p in matrix],
                         'mass': [p.mass for p in matrix],
                         'position': [p.position for p in matrix],
                         }
            _write_datasets(self.trajectory, group, matrix_h5)

        # Cell
        group = '/initialstate/cell/'
        if system.cell is not None:
            self.trajectory.create_group_safe(group)
            self.trajectory[group + 'sidebox'] = system.cell.side

        # Thermostat
        group = '/initialstate/thermostat/'
        if system.thermostat is not None:
            self.trajectory.create_group_safe(group)
            self.trajectory[group + 'temperature'] = system.thermostat.temperature
            self.trajectory[group + 'type'] = system.thermostat.name
            self.trajectory[group + 'mass'] = system.thermostat.mass
            self.trajectory[group + 'collision_period'] = system.thermostat.collision_period

        # Interaction
        if system.interaction is not None:
            if type(system.interaction) is list:
                raise TypeError('interaction must be list')
            self.trajectory.copy(system.interaction, '/initialstate/interaction/')
            # self.write_interaction([system.interaction])

    def write_interaction(self, interaction):
        rgr = '/initialstate/interaction/'
        self.trajectory.create_group_safe('/initialstate/')
        # If the group exisist we delete it. This does not actual clear space in h5 file.
        # We could do it on a dataset basis via require_dataset, or visit the group and delete everything.
        try:
            del self.trajectory['/initialstate/interaction/']
        except:
            pass
        self.trajectory.create_group_safe('/initialstate/interaction/')
        self.trajectory[rgr + 'number_of_interactions'] = [len(interaction)]
        for i, term in enumerate(interaction):
            igr = rgr + '/interaction_%d/' % (i+1)
            # Numbering is fortran style for backward compatibility
            self.trajectory.create_group_safe(igr)
            self.trajectory[igr + 'interaction_type'] = [term.name]
            self.trajectory[igr + 'number_of_potentials'] = [len(term.potential)]
            for j, phi in enumerate(term.potential):
                pgr = igr + '/potential_%d/' % (j+1)
                self.trajectory.create_group_safe(pgr)
                self.trajectory[pgr + 'potential'] = [str(phi)]
                self.trajectory[pgr + 'interacting_bodies'] = [phi.interacting_bodies]
                self.trajectory[pgr + 'interacting_species'] = phi.species
                self.trajectory[pgr + 'parameters_number'] = [len(phi.params)]
                self.trajectory[pgr + 'parameters_name'] = sorted(phi.params.keys())
                self.trajectory[pgr + 'parameters'] = [phi.params[k] for k in sorted(phi.params.keys())]
                self.trajectory[pgr + 'cutoff_scheme'] = [phi.cutoff.scheme]
                self.trajectory[pgr + 'cutoff_radius'] = [phi.cutoff.radius]
                self.trajectory[pgr + 'lookup_points'] = [phi.npoints]

    def write_sample(self, system, step):
        self.trajectory.create_group_safe('/trajectory')
        self.trajectory.create_group_safe('/trajectory/realtime')
        self.trajectory.create_group_safe('/trajectory/realtime/stepindex')
        self.trajectory.create_group_safe('/trajectory/realtime/sampleindex')

        frame = len(self.steps) + 1
        csample = '/sample_%7.7i' % frame

        try:
            self.trajectory['/trajectory/realtime/stepindex' + csample] = [step]
            self.trajectory['/trajectory/realtime/sampleindex' + csample] = [frame]
        except RuntimeError:
            _log.error('error when writing step %s sample %s to file %s', step, frame, self.filename)
            raise

        if system.particle is not None:
            self.trajectory.create_group_safe('/trajectory/particle')
            if 'position' in self.fields:
                self.trajectory.create_group_safe('/trajectory/particle/position')
                self.trajectory['/trajectory/particle/position' + csample] = [p.position for p in system.particle]
            if 'velocity' in self.fields:
                self.trajectory['/trajectory/particle/velocity' + csample] = [p.velocity for p in system.particle]
                self.trajectory.create_group_safe('/trajectory/particle/velocity')
            if 'radius' in self.fields:
                self.trajectory.create_group_safe('/trajectory/particle/radius')
                self.trajectory['/trajectory/particle/radius' + csample] = [p.radius for p in system.particle]
            if 'species' in self.fields:
                self.trajectory.create_group_safe('/trajectory/particle/species')
                self.trajectory['/trajectory/particle/species' + csample] = ['%-3s' % p.species for p in system.particle]
                from atooms.system.particle import distinct_species
                ids = distinct_species(system.particle)
                self.trajectory['/trajectory/particle/ids' + csample] = [1+ids.index(p.species) for p in system.particle]

        if system.cell is not None:
            self.trajectory.create_group_safe('/trajectory/cell')
            if 'cell' in self.fields:
                self.trajectory.create_group_safe('/trajectory/cell/sidebox')
                self.trajectory['/trajectory/cell/sidebox' + csample] = system.cell.side

    def read_init(self):
        # read particles
        group = self.trajectory['/initialstate/particle']
        n = self.trajectory['/initialstate/particle/number_of_particles'].value[0]
        rad = None
        for entry in group:
            # TODO: refactor this
            if entry == 'element':
                spe = group[entry][:]
            if entry == 'mass':
                mas = group[entry][:]
            if entry == 'position':
                pos = group[entry][:]
            if entry == 'velocity':
                vel = group[entry][:]
            if entry == 'radius':
                rad = group[entry][:]
        if rad is not None:
            particle = [Particle(species=spe[i].strip(), mass=mas[i],
                                 position=pos[i, :], velocity=vel[i, :],
                                 radius=rad[i]) for i in range(n)]
        else:
            particle = [Particle(species=spe[i].strip(), mass=mas[i],
                                 position=pos[i, :], velocity=vel[i, :])
                        for i in range(n)]

        # read cell
        group = self.trajectory['/initialstate/cell']
        for entry in group:
            if entry == 'sidebox':
                sidebox = group[entry][:]
        cell = Cell(sidebox)

        # read interaction as h5 group
        try:
            interaction = self.trajectory['/initialstate/interaction']
        except:
            interaction = None

        # build system
        self._system = System(particle, cell, interaction)

        # read matrix
        if 'matrix' in self.trajectory['/initialstate']:
            group = self.trajectory['/initialstate/matrix']
            for entry in group:
                if entry == 'element':
                    spe = group[entry][:]
                if entry == 'mass':
                    mas = group[entry][:]
                if entry == 'position':
                    pos = group[entry][:]
            matrix = [Particle(species=spe[i].strip(), mass=mas[i],
                               position=pos[i, :])
                      for i in range(len(spe))]
            self._system.matrix = copy.deepcopy(matrix)

        return self._system

    def read_interaction(self):
        # read interaction terms
        if 'interaction' not in self.trajectory['/initialstate']:
            return None

        n = self.trajectory['/initialstate/interaction/number_of_interactions'][0]
        interactions = []
        for i in range(n):
            g = '/initialstate/interaction/interaction_%d/' % (i+1)
            np = self.trajectory[g + 'number_of_potentials'][0]
            name = self.trajectory[g + 'interaction_type'][0]
            potentials = []
            for j in range(np):
                sg = self.trajectory[g + 'potential_%d/' % (j+1)]
                # params = {k:v for k, v in zip(sg['parameters_name'][:], sg['parameters'][:])}
                # make it compatible with 2.6
                params = {}
                for k, v in zip(sg['parameters_name'][:], sg['parameters'][:]):
                    params[k] = v
                p = PairPotential(sg['potential'][0], params, sg['interacting_species'][:],
                                  CutOff(sg['cutoff_scheme'][0], sg['cutoff_radius'][0]),
                                  sg['lookup_points'][0])

                potentials.append(p)
            interactions.append(Interaction(potentials, name))
        return interactions

    def read_sample(self, frame, unfolded=False):
        # TODO: due to unfolded argument this differs from the base class method Can we drop this?
        # We must increase frame by 1 if we iterate over frames with len().
        # This is some convention to be fixed once and for all
        # TODO: read cell on the fly NPT
        # TODO: are keys cached?
        csample = '/' + self.trajectory['/trajectory/realtime/stepindex'].keys()[frame]

        # read particles
        group = self.trajectory['/trajectory/particle']
        if unfolded:
            if 'position_unfolded' not in group:
                raise NotImplementedError('cannot unfold like this, use decorator instead')
            else:
                # fix for unfolded positions that were not written at the first step
                # should be fixed once and for all in md.x
                if frame == 0:
                    pos = self.trajectory['/initialstate/particle/position'][:]
                else:
                    pos = group['position_unfolded' + csample][:]
        else:
            pos = group['position' + csample][:]

        try:
            vel = group['velocity' + csample][:]
        except:
            vel = numpy.zeros([len(pos), ndim])

        # Dynamic properties
        p = []
        for r, v in zip(pos, vel):
            p.append(Particle(position=r, velocity=v))

        # Static properties
        # TODO: optimize, this takes quite some additional time, almost x2
        for pi, r in zip(p, self._system.particle):
            # TODO: if id changes dynamically (like in swap) we will miss it. We should update them after this loop!
            pi.mass = r.mass
            pi.species = r.species
            pi.radius = r.radius

        # Try update radii. This must be done after setting defaults.
        try:
            r = group['radius' + csample][:]
            for i, pi in enumerate(p):
                pi.radius = r[i]
            if 'radius' not in self.fields:
                self.fields.append('radius')
        except KeyError:
            if 'radius' in self.fields:
                self.fields.remove('radius')

        # Read cell
        group = self.trajectory['/trajectory/cell']
        side = group['sidebox' + csample][:]
        # This fixes an issue with some hdf5 trajectories that stored
        # cell as (1,3) array
        if len(side.shape) == 2:
            side = side[0]
        self._system.cell.side = side

        # Read also interaction.
        has_int = True
        try:
            group = self.trajectory['/trajectory/interaction']
        except:
            has_int = False

        if has_int:
            self._system.interaction.total_energy = group['energy' + csample][0]
            self._system.interaction.total_virial = group['virial' + csample][0]
            self._system.interaction.total_stress = group['stress' + csample][:]

        return System(p, self._system.cell, self._system.interaction)
