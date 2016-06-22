# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import numpy
import h5py
from base import TrajectoryBase
from atooms import ndim
from atooms.system import System
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.interaction.interaction import Interaction
from atooms.potential.potential import PairPotential
from atooms.potential.cutoff import CutOff

class _SafeFile(h5py.File):
    # TODO: decorate hdf5 class so that error messages contain the path of the offending file
    def create_group_safe(self, group):
        # TODO: recursively create all h5 groups that do not exist?
        # TODO: actually redefine create_group unless this is a serious performace issue?
        if not group in self:
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

           
class TrajectoryHDF5(TrajectoryBase):

    """ Trajectory layout based on HDF5 libraries """

    suffix = 'h5'

    def __init__(self, filename, mode='r+'):
        TrajectoryBase.__init__(self, filename, mode)

        self.general_info = {}
        self._grandcanonical = False
        self._system = None
        self.fmt = ['position', 'velocity', 'cell']

        if self.mode == 'r' or self.mode == 'r+':
            self.trajectory = h5py.File(self.filename, mode)
            # gather general info on file
            for entry in self.trajectory['/']:
                if type(self.trajectory[entry]) == h5py.highlevel.Dataset:
                    self.general_info[entry] = self.trajectory[entry]
            try:
                # get steps list (could be cached and put in init_read())           
                self.steps = [d[0] for d in self.trajectory['trajectory/realtime/stepindex'].values()]
                # private list of samples. This solves the problem that samples may start from 0
                # or 1 depending on the code that initially produced the data
                # TODO: can we drop this for performance?
                self._samples = [d[0] for d in self.trajectory['trajectory/realtime/sampleindex'].values()]
            except KeyError:
                self.steps = []
                self._samples = []

            # Block period
            self._block_period = self.read_blockperiod()
            
        elif self.mode == 'w' or self.mode == 'r+' or self.mode == "w-":
            self.trajectory = _SafeFile(self.filename, self.mode)

        else:
            raise ValueError('Specify mode (r/w) for file %s (invalid: %s)' % (self.filename, self.mode))

    def close(self):
        self.trajectory.close()

    def read_timestep(self):
        try:
            return self.trajectory['trajectory/realtime/timestep'][0]
        except:
            print 'Warning: could not find dt in hdf5 file %s' % self.filename
            return 1.0

    def write_timestep(self, value):
        self.trajectory.create_group_safe('/trajectory')
        self.trajectory.create_group_safe('/trajectory/realtime')
        self.trajectory['trajectory/realtime/timestep'] = [value]

    def read_blockperiod(self):
        try:
            return self.trajectory['trajectory/realtime/block_period'][0]
        except:
            return None

    def write_blockperiod(self, value):
        self.trajectory.create_group_safe('/trajectory')
        self.trajectory.create_group_safe('/trajectory/realtime')
        self.trajectory['trajectory/realtime/block_period'] = [value]

    def write_init(self, system):
        self.trajectory.create_group_safe('/initialstate')
        self.trajectory['DIMENSIONS'] = [3]
        self.trajectory['NAME_SYS'] = ['Unknown'] #system.name
        self.trajectory['VERSION_TRJ'] = ['1.2']
        self.trajectory['VERSION_MD'] = ['X.X.X']

        # Particles
        group = '/initialstate/particle/'
        if system.particle is not None:
            self.trajectory.create_group_safe(group)
            particle = system.particle

            # Check that species id are ok (problems might arise when converting from RUMD
            if min([p.id for p in particle]) < 1:
                raise ValueError('Particles ids are smaller than 1. Use NormalizeId decorator to fix this.')
 
            particle_h5 = {'number_of_species': [len(list(set([p.id for p in particle])))], #particle.numberSpecies(),
                           'number_of_particles': [len(particle)],
                           'identity': [p.id for p in particle],
                           'element': ['%3s' % p.name for p in particle],
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
            matrix_h5 = {'type' : [''],
                         'id' : [0],
                         'number_of_species': [len(list(set([p.id for p in matrix])))],
                         'number_of_particles': [len(matrix)],
                         'identity': [p.id for p in matrix],
                         'element': ['%3s' % p.name for p in matrix],
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
            #self.write_interaction([system.interaction])
        
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
            self.trajectory.create_group_safe(igr) # numbering is from one to comply with atooms
            self.trajectory[igr + 'interaction_type'] = [term.name]
            self.trajectory[igr + 'number_of_potentials'] = [len(term.potential)]
            for j, phi in enumerate(term.potential):
                pgr = igr + '/potential_%d/' % (j+1)
                self.trajectory.create_group_safe(pgr)
                self.trajectory[pgr + 'potential'] = [phi.name]
                self.trajectory[pgr + 'interacting_bodies'] = [phi.interacting_bodies]
                self.trajectory[pgr + 'interacting_species'] = phi.species
                self.trajectory[pgr + 'parameters_number'] = [len(phi.params)]
                self.trajectory[pgr + 'parameters_name'] = sorted(phi.params.keys())
                self.trajectory[pgr + 'parameters'] = [phi.params[k] for k in sorted(phi.params.keys())]
                self.trajectory[pgr + 'cutoff_scheme'] = [phi.cutoff.name]
                self.trajectory[pgr + 'cutoff_radius'] = [phi.cutoff.formal_radius]
                self.trajectory[pgr + 'lookup_points'] = [phi.npoints]

    def write_sample(self, system, step):
        self.trajectory.create_group_safe('/trajectory')
        self.trajectory.create_group_safe('/trajectory/realtime')
        self.trajectory.create_group_safe('/trajectory/realtime/stepindex')
        self.trajectory.create_group_safe('/trajectory/realtime/sampleindex')

        sample = len(self.steps) + 1
        csample = '/sample_%7.7i' % sample

        self.trajectory['/trajectory/realtime/stepindex' + csample] = [step]
        self.trajectory['/trajectory/realtime/sampleindex' + csample] = [sample]

        if system.particle != None:
            self.trajectory.create_group_safe('/trajectory/particle')
            if 'position' in self.fmt:
                self.trajectory.create_group_safe('/trajectory/particle/position')
                self.trajectory['/trajectory/particle/position' + csample] = [p.position for p in system.particle]
            if 'velocity' in self.fmt:
                self.trajectory['/trajectory/particle/velocity' + csample] = [p.velocity for p in system.particle]
                self.trajectory.create_group_safe('/trajectory/particle/velocity')
            if 'radius' in self.fmt:
                self.trajectory.create_group_safe('/trajectory/particle/radius')
                self.trajectory['/trajectory/particle/radius' + csample] = [p.radius for p in system.particle]

        if system.cell != None:
            self.trajectory.create_group_safe('/trajectory/cell')
            if 'cell' in self.fmt:
                self.trajectory.create_group_safe('/trajectory/cell/sidebox')
                self.trajectory['/trajectory/cell/sidebox' + csample] = [system.cell.side]

    def read_init(self):
        # read particles
        group = self.trajectory['/initialstate/particle']
        n = self.trajectory['/initialstate/particle/number_of_particles'].value[0]
        rad = None
        for entry in group:
            if entry == 'identity': spe = group[entry][:]
            if entry == 'element' : ele = group[entry][:]
            if entry == 'mass'    : mas = group[entry][:]
            if entry == 'position': pos = group[entry][:]
            if entry == 'velocity': vel = group[entry][:]
            if entry == 'radius'  : rad = group[entry][:]
        if rad is not None:
            particle = [Particle(spe[i],ele[i],mas[i],pos[i,:],vel[i,:],rad[i]) for i in range(n)]
        else:
            particle = [Particle(spe[i],ele[i],mas[i],pos[i,:],vel[i,:]) for i in range(n)]

        # read cell
        group = self.trajectory['/initialstate/cell']
        for entry in group:
            if entry == 'sidebox' : sidebox = group[entry][:]
        cell = Cell(sidebox)

        # read interaction as h5 group
        try:
            interaction = self.trajectory['/initialstate/interaction']
        except:
            interaction = None

        # build system
        self._system = System(particle, cell, interaction)

        # read matrix
        if 'matrix' in  self.trajectory['/initialstate']:
            group = self.trajectory['/initialstate/matrix']
            for entry in group:
                if entry == 'identity': spe = group[entry][:]
                if entry == 'element' : ele = group[entry][:]
                if entry == 'mass'    : mas = group[entry][:]
                if entry == 'position': pos = group[entry][:]
            matrix = [Particle(spe[i],ele[i],mas[i],pos[i,:]) for i in range(len(spe))] 
            self._system.add_porous_matrix(matrix)

        return self._system

    def read_interaction(self):
        # read interaction terms
        if not 'interaction' in self.trajectory['/initialstate']:
            return None

        group = self.trajectory['/initialstate/interaction']
        n = self.trajectory['/initialstate/interaction/number_of_interactions'][0]
        interactions = []
        for i in range(n):
            g = '/initialstate/interaction/interaction_%d/' % (i+1)
            np = self.trajectory[g + 'number_of_potentials'][0]
            name = self.trajectory[g + 'interaction_type'][0]
            potentials = []
            for j in range(np):
                sg = self.trajectory[g + 'potential_%d/' % (j+1)]
                #params = {k:v for k, v in zip(sg['parameters_name'][:], sg['parameters'][:])}
                # make it compatible with 2.6
                params = {}
                for k, v in zip(sg['parameters_name'][:], sg['parameters'][:]):
                    params[k] = v
                p = PairPotential(sg['potential'][0], params, sg['interacting_species'][:], 
                                  CutOff(sg['cutoff_scheme'][0], sg['cutoff_radius'][0]),
                                  sg['lookup_points'][0])

                potentials.append(p)
            interactions.append(Interaction(name, potentials))
        return interactions

    def read_sample(self, sample, unfolded=False):
        # We must increase sample by 1 if we iterate over samples with len().
        # This is some convention to be fixed once and for all
        isample = self._samples[sample]
        csample = '/sample_%7.7i' % isample
        # read particles
        group = self.trajectory['/trajectory/particle']
        if unfolded:
            if not 'position_unfolded' in group:
                raise NotImplementedError('cannot unfold like this, use decorator instead')
            else:
                # fix for unfolded positions that were not written at the first step
                # should be fixed once and for all in md.x
                if sample == 0:
                    pos = self.trajectory['/initialstate/particle/position'][:]
                else:
                    pos = group['position_unfolded' + csample][:]
        else:
            pos = group['position' + csample][:]

        try:
            vel = group['velocity' + csample][:]
        except:
            vel = numpy.zeros([len(pos),ndim])

        # Dynamic properties
        p = []
        for r, v in zip(pos, vel):
            p.append(Particle(position=r, velocity=v))

        # Static properties
        # TODO: optimize, this takes quite some additional time, almost x2
        for pi, r in zip(p, self._system.particle):
            # TODO: if id changes dynamically (like in swap) we will miss it. We should update them after this loop
            pi.id = r.id
            pi.mass = r.mass
            pi.name = r.name
            pi.radius = r.radius

        # Try update radii. This must be done after setting defaults
        try:
            r = group['radius' + csample][:]           
            for i, pi in enumerate(p):
                pi.radius = r[i]
            self.include(['radius'])
        except:
            self.exclude(['radius'])

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

    def add_interaction(self, ff):

        """Add interaction from a fortran forcefield file"""

        import os

        pid = os.getpid()
        f_ref = '/tmp/cnv_%s.h5' % pid
        # TODO: we can cache a ref file if ff is the same
        if not os.path.exists(ff):
            raise IOError('forcefield file does not exist %s' % ff)
        os.system('system.x -n 2 -f %s %s 1>/dev/null 2>/dev/null' % (ff, f_ref))

        # Add interaction
        ref = h5py.File(f_ref, 'r')
        # Make sure interaction does not exist
        try:
            del self.trajectory['initialstate/interaction']
        except:
            pass
        self.trajectory.copy(ref['initialstate/interaction'], 'initialstate/interaction')

        # Cleanup
        ref.close()
        os.remove(f_ref)

if __name__ == '__main__':
    add_interaction_hdf5('/tmp/1', '/home/coslo/codes/atooms/forcefields/LJ.ff')
