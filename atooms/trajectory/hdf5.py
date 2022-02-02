# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""HDF5 trajectory format."""

import numpy
import h5py
import logging
import warnings

from .base import TrajectoryBase
from atooms.core import ndim
from atooms.system import System
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system.interaction import Interaction

_log = logging.getLogger(__name__)


class _CutOff(object):

    def __init__(self, scheme, radius):
        self.scheme = scheme
        self.radius = radius
        self.rcutsq = radius**2
        self.radius_mid = radius
        self.radius_mid_sq = radius**2


class _PairPotential(object):

    interacting_bodies = 2

    def __init__(self, func, params, species, cutoff=None,
                 npoints=20000):
        self.func = func
        self.params = params
        self.species = species
        self.cutoff = cutoff
        self.npoints = npoints


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

    """In-house trajectory layout in HDF5 format."""

    suffix = 'h5'

    def __init__(self, filename, mode='r+'):
        super(TrajectoryHDF5, self).__init__(filename, mode)

        self.general_info = {}
        self._grandcanonical = False
        self._system = None
        # TODO: move this to write mode
        self.variables = ['particle.position', 'particle.velocity']
        if self.mode == 'r' or self.mode == 'r+':
            self._file = h5py.File(self.filename, mode)
            # gather general info on file
            for entry in self._file['/']:
                if type(self._file[entry]) == h5py.Dataset:
                    self.general_info[entry] = self._file[entry]

        elif self.mode == 'w' or self.mode == 'r+' or self.mode == "w-":
            self._file = _SafeFile(self.filename, self.mode)

        else:
            raise ValueError('Specify mode (r/w) for file %s (invalid: %s)' % (self.filename, self.mode))

    def read_steps(self):
        return [d[0] for d in self._file['trajectory/realtime/stepindex'].values()]

    def read_len(self):
        return len(self._file['trajectory/realtime/stepindex'])

    def close(self):
        try:
            self._file.close()
        except ValueError:
            _log.error('file %s already closed', self.filename)
            raise

    def read_timestep(self):
        try:
            return self._file['trajectory/realtime/timestep'][0]
        except:
            _log.warning('no time step in file %s, set to 1', self.filename)
            return 1.0

    def write_timestep(self, value):
        self._file.create_group_safe('/trajectory')
        self._file.create_group_safe('/trajectory/realtime')
        self._file['trajectory/realtime/timestep'] = [value]

    def read_block_size(self):
        try:
            return self._file['trajectory/realtime/block_period'][0]
        except:
            return None

    def write_block_size(self, value):
        self._file.create_group_safe('/trajectory')
        self._file.create_group_safe('/trajectory/realtime')
        self._file['trajectory/realtime/block_period'] = [value]

    def write_init(self, system):
        from atooms.system.particle import distinct_species
        self._file.create_group_safe('/initialstate')
        self._file['DIMENSIONS'] = [3]
        self._file['NAME_SYS'] = [b'Unknown']
        self._file['VERSION_TRJ'] = [b'1.2']
        self._file['VERSION_MD'] = [b'X.X.X']

        # Particles
        group = '/initialstate/particle/'
        if system.particle is not None:
            self._file.create_group_safe(group)
            particle = system.particle
            species = distinct_species(particle)
            particle_h5 = {'number_of_species': [len(species)],
                           'number_of_particles': [len(particle)],
                           'identity': [species.index(p.species)+1 for p in particle],
                           'element': [str(p.species) for p in particle],
                           'mass': [p.mass for p in particle],
                           'radius': [p.radius for p in particle],
                           'position': [p.position for p in particle],
                           'velocity': [p.velocity for p in particle],
                           }
            _write_datasets(self._file, group, particle_h5)

        # Cell
        group = '/initialstate/cell/'
        if system.cell is not None:
            self._file.create_group_safe(group)
            self._file[group + 'sidebox'] = system.cell.side

        # Thermostat
        group = '/initialstate/thermostat/'
        if system.thermostat is not None:
            self._file.create_group_safe(group)
            self._file[group + 'temperature'] = system.thermostat.temperature
            self._file[group + 'type'] = system.thermostat.name
            self._file[group + 'mass'] = system.thermostat.mass
            self._file[group + 'collision_period'] = system.thermostat.collision_period

        # Interaction
        if system.interaction is not None:
            if type(system.interaction) is list:
                raise TypeError('cannot handle more than one interaction')
            self.write_interaction([system.interaction])

    def write_interaction(self, interaction):
        rgr = '/initialstate/interaction/'
        self._file.create_group_safe('/initialstate/')
        # If the group exisist we delete it. This does not actual clear space in h5 file.
        # We could do it on a dataset basis via require_dataset, or visit the group and delete everything.
        try:
            del self._file['/initialstate/interaction/']
        except:
            pass
        self._file.create_group_safe('/initialstate/interaction/')
        self._file[rgr + 'number_of_interactions'] = [len(interaction)]
        for i, term in enumerate(interaction):
            igr = rgr + '/interaction_%d/' % (i+1)
            # Numbering is fortran style for backward compatibility
            self._file.create_group_safe(igr)
            self._file[igr + 'interaction_type'] = [term.name.encode()]
            self._file[igr + 'number_of_potentials'] = [len(term.potential)]
            for j, phi in enumerate(term.potential):
                pgr = igr + '/potential_%d/' % (j+1)
                self._file.create_group_safe(pgr)
                self._file[pgr + 'potential'] = [str(phi).encode()]
                self._file[pgr + 'interacting_bodies'] = [phi.interacting_bodies]
                self._file[pgr + 'interacting_species'] = phi.species
                self._file[pgr + 'parameters_number'] = [len(phi.params)]
                self._file[pgr + 'parameters_name'] = [_ for _ in sorted(phi.params.keys())]
                self._file[pgr + 'parameters'] = [phi.params[k] for k in sorted(phi.params.keys())]
                self._file[pgr + 'cutoff_scheme'] = [phi.cutoff.scheme.encode()]
                self._file[pgr + 'cutoff_radius'] = [phi.cutoff.radius]
                self._file[pgr + 'lookup_points'] = [phi.npoints]

    def write_system(self, system, step):
        variables = self.variables

        self._file.create_group_safe('/trajectory')
        self._file.create_group_safe('/trajectory/realtime')
        self._file.create_group_safe('/trajectory/realtime/stepindex')
        self._file.create_group_safe('/trajectory/realtime/sampleindex')

        frame = len(self.steps) + 1
        csample = '/sample_%7.7i' % frame

        try:
            self._file['/trajectory/realtime/stepindex' + csample] = [step]
            self._file['/trajectory/realtime/sampleindex' + csample] = [frame]
        except RuntimeError:
            _log.error('error when writing step %s sample %s to file %s', step, frame, self.filename)
            raise

        if system.particle is not None:
            self._file.create_group_safe('/trajectory/particle')
            if 'particle.position' in variables:
                self._file.create_group_safe('/trajectory/particle/position')
                self._file['/trajectory/particle/position' + csample] = [p.position for p in system.particle]
            if 'particle.velocity' in variables:
                self._file['/trajectory/particle/velocity' + csample] = [p.velocity for p in system.particle]
                self._file.create_group_safe('/trajectory/particle/velocity')
            if 'particle.radius' in variables:
                self._file.create_group_safe('/trajectory/particle/radius')
                self._file['/trajectory/particle/radius' + csample] = [p.radius for p in system.particle]
            if 'particle.species' in variables:
                self._file.create_group_safe('/trajectory/particle/species')
                data = ['%-3s' % p.species for p in system.particle]
                self._file['/trajectory/particle/species' + csample] = [_.encode() for _ in data]
                from atooms.system.particle import distinct_species
                ids = distinct_species(system.particle)
                self._file['/trajectory/particle/ids' + csample] = [1+ids.index(p.species) for p in system.particle]

        if system.cell is not None:
            self._file.create_group_safe('/trajectory/cell')
            self._file.create_group_safe('/trajectory/cell/sidebox')
            self._file['/trajectory/cell/sidebox' + csample] = system.cell.side

    def read_init(self):
        # read particles
        group = self._file['/initialstate/particle']
        n = self._file['/initialstate/particle/number_of_particles'][0]
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
            particle = [Particle(species=spe[i].decode().strip(), mass=mas[i],
                                 position=pos[i, :], velocity=vel[i, :],
                                 radius=rad[i]) for i in range(n)]
        else:
            particle = [Particle(species=spe[i].decode().strip(), mass=mas[i],
                                 position=pos[i, :], velocity=vel[i, :])
                        for i in range(n)]

        # read cell
        group = self._file['/initialstate/cell']
        for entry in group:
            if entry == 'sidebox':
                sidebox = group[entry][:]
        cell = Cell(sidebox)

        # read interaction
        interaction = self.read_interaction()

        # build system
        self._system = System(particle, cell, interaction)

        return self._system

    def read_interaction(self):
        # read interaction terms
        if 'interaction' not in self._file['/initialstate']:
            return None

        n = self._file['/initialstate/interaction/number_of_interactions'][0]
        if n > 1:
            warnings.warn('can only read one interaction term')

        i = 0
        g = '/initialstate/interaction/interaction_%d/' % (i+1)
        np = self._file[g + 'number_of_potentials'][0]
        name = self._file[g + 'interaction_type'][0].decode()
        potentials = []
        for j in range(np):
            sg = self._file[g + 'potential_%d/' % (j+1)]
            # params = {k:v for k, v in zip(sg['parameters_name'][:], sg['parameters'][:])}
            # make it compatible with 2.6
            params = {}
            for k, v in zip(sg['parameters_name'][:], sg['parameters'][:]):
                params[k] = v
            p = _PairPotential(sg['potential'][0].decode(), params,
                               sg['interacting_species'][:],
                               _CutOff(sg['cutoff_scheme'][0].decode(),
                                       sg['cutoff_radius'][0]))
            if 'lookup_points' in sg:
                p.npoints = sg['lookup_points']

            potentials.append(p)
        interaction = Interaction()
        interaction.potential = potentials
        interaction.name = name
        return interaction

    def read_system(self, frame):
        # TODO: refactor reading particle variables
        # TODO: removing variables in read_system() is pointless, it should be done before
        # We must increase frame by 1 if we iterate over frames with len().
        # This is some convention to be fixed once and for all
        # TODO: read cell on the fly NPT
        keys = list(self._file['/trajectory/realtime/stepindex'].keys())
        csample = '/' + keys[frame]
        # read particles
        group = self._file['/trajectory/particle']
        pos = group['position' + csample][:]

        if 'position_unfolded' in group:
            # fix for unfolded positions that were not written at the first step
            # should be fixed once and for all in md.x
            if frame == 0:
                pos_unf = self._file['/initialstate/particle/position'][:]
            else:
                pos_unf = group['position_unfolded' + csample][:]
            if 'particle.position_unfolded' not in self.variables:
                self.variables.append('particle.position_unfolded')

        if 'velocity' in group and len(group['velocity']) > 0:
            vel = group['velocity' + csample][:]
        else:
            vel = numpy.zeros((len(pos), ndim))
            if 'particle.velocity' not in self.variables:
                self.variables.remove('particle.velocity')

        # Dynamic properties
        p = []
        for r, v in zip(pos, vel):
            p.append(Particle(position=numpy.array(r, dtype='float64'),
                              velocity=numpy.array(v, dtype='float64')))

        # Static properties
        # TODO: optimize, this takes quite some additional time, almost x2
        for pi, po in zip(p, self._system.particle):
            pi.mass = po.mass
            pi.species = po.species
            pi.radius = po.radius

        # Try update radii. This must be done after setting defaults.
        try:
            r = group['radius' + csample][:]
            for i, pi in enumerate(p):
                pi.radius = r[i]
            if 'particle.radius' not in self.variables:
                self.variables.append('particle.radius')
        except KeyError:
            if 'particle.radius' in self.variables:
                self.variables.remove('particle.radius')

        # Try update species. This must be done after setting defaults.
        try:
            spe = group['species' + csample][:]
            for i, pi in enumerate(p):
                pi.species = spe[i].decode().strip()
            if 'particle.species' not in self.variables:
                self.variables.append('particle.species')
        except KeyError:
            if 'particle.species' in self.variables:
                self.variables.remove('particle.species')

        # Read cell
        group = self._file['/trajectory/cell']
        side = group['sidebox' + csample][:]
        # This fixes an issue with some hdf5 trajectories that stored
        # cell as (1,3) array
        if len(side.shape) == 2:
            side = side[0]
        self._system.cell.side = side

        # Read also interaction.
        has_int = True
        try:
            group = self._file['/trajectory/interaction']
        except:
            has_int = False

        if has_int:
            self._system.interaction.total_energy = group['energy' + csample][0]
            self._system.interaction.total_virial = group['virial' + csample][0]
            self._system.interaction.total_stress = group['stress' + csample][:]

        return System(p, self._system.cell, self._system.interaction)
